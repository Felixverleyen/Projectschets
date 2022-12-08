#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "nbody.hpp"
using namespace std;


double G=1;

Vec a(int i, nbody sim) {
    
    Vec ri= sim.r(i);
    Vec ar={0,0,0};

    for (int j =0; j<sim.bodies(); ++j){
        Vec rj = sim.r(j);
        double mj= sim.m(j);
        double mi= sim.m(i);
        Vec afst= ri-rj;

        if (i!=j){
            ar-= G*mj*afst/afst.norm3();
            } 
        }
    return ar;}

double Energy(nbody sim){
    double E=0;

    for (int i=0; i<sim.bodies();++i){
        Vec ri = sim.r(i);
        double mi=sim.m(i);
        Vec vi= sim.v(i);
        E+=mi*vi.norm2()/2;

        for  (int j =0; j<sim.bodies(); ++j){
            Vec rj = sim.r(j);
            double mj= sim.m(j);
            Vec afst= ri-rj;

            if (i!=j){
                E-= 1/2*G*mi*mj/afst.norm();
                }
            }
    }

     return E;}

 // integratoren:

 nbody verlet_halfstep(double h, nbody sim_half){
    for (int i=0; i<sim_half.bodies(); i++){
        sim_half.swap_v(i, sim_half.v(i) + 0.5*h*a(i,sim_half));
    }
    return sim_half;
 }

 nbody verlet_step(double h, nbody sim, nbody sim_half){
    for (int i=0; i<sim.bodies(); i++) {
        sim.swap_r(i, sim.r(i) + h*sim_half.v(i));
        sim.swap_v(i, sim_half.v(i) + 0.5*h*a(i,sim));
    }
    return sim;
 }

// velocity verlet
void verlet(double h, double time, nbody sim) {
    ofstream outfile1("verlet.txt");
    ofstream outfile2("VerletE.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15); 
    int steps = int(time / h);
    
    for (double t=0; t<steps; t++) {
        nbody sim_half = verlet_halfstep(h, sim);
        sim = verlet_step(h, sim, sim_half);
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        double E = Energy(sim);
        outfile2 << E << '\n';
    }

    outfile1.close();
    outfile2.close();
}

int main(){
    string initial_i;
    nbody sim;
    int N = 0;
    int l = 0;
    fstream initialNfile("Initial_cond_Zonnestelsel.txt");

    while (getline(initialNfile, initial_i)){
        if(l>4){
            double m = stod(initial_i.substr(2,4));

            double rx = stod(initial_i.substr(7, 4));
            double ry = stod(initial_i.substr(12, 4));
            double rz = stod(initial_i.substr(17, 4));
            double vx = stod(initial_i.substr(22, 4));
            double vy = stod(initial_i.substr(27, 4));
            double vz = stod(initial_i.substr(32, 4));

                
            Vec pos{rx, ry, rz};
            Vec vel{vx,vy,vz};

            sim.add_mass(m);
            sim.add_pos(pos);
            sim.add_vel(vel);

            ++N;
        }

        ++l;
        };

    sim.set_N(N);
    initialNfile.close();
    double h = 1e-5;
    double time = 1;
    verlet(h, time, sim);
}



#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <string.h>
#include "nbody.hpp"
using namespace std;


double G=8;
double mu=0.01;
double h=0.0001;
int const N=4;

   
void print(Vec a)
{ cout << a.x() << ' ' << a.y() << ' '<< a.z() << endl; }

Vec a(int i, nbody sim) {
    
    Vec ri= sim.r(i);
    Vec ar={0,0,0};

    for (int j =0; j<N; ++j){
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

    for (int i=0; i<N;++i){
        Vec ri = sim.r(i);
        double mi=sim.m(i);
        Vec vi= sim.v(i);
        E+=mi*vi.norm2()/2;

        for  (int j =0; j<N; ++j){
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

    //Runge-Kutta
nbody RK4N_step(double h, nbody sim){
    nbody sim1 = sim;
    nbody sim2 = sim;
    nbody sim3 = sim;

    for (int i=0; i < sim.bodies(); i++){
        
        Vec kx1 = h * sim.v(i);
        Vec kv1 = h * a(i, sim);

        Vec lx1 = sim.r(i) + 0.5 * kx1;
        Vec lv1 = sim.v(i) + 0.5 * kv1;
        sim1.swap_r(i, lx1);
        sim1.swap_v(i, lv1);
        
        
        sim.swap_r(i, sim.r(i) + 1./6. * kx1);
        sim.swap_v(i, sim.v(i) + 1./6. * kv1);
        }
    

    for (int i=0; i < sim.bodies(); i++){
        Vec kx2 = h* sim1.v(i);
        Vec kv2 = h * a(i, sim1);

        Vec lx2 = sim1.r(i) + 0.5 * kx2;
        Vec lv2 = sim.v(i) + 0.5 * kv2;
        sim2.swap_r(i, lx2);
        sim2.swap_v(i, lv2);
        
        sim.swap_r(i, sim.r(i) + 1./6. * 2 *  kx2);
        sim.swap_v(i, sim.v(i) + 1./6. * 2 * kv2);
        }
    
    for (int i=0; i < sim.bodies(); i++){

        Vec kx3 = h * sim2.v(i);
        Vec kv3 = h * a(i, sim2);

        Vec lx3 = sim2.r(i) + kx3;
        Vec lv3 = sim2.v(i) + kv3;
        sim3.swap_r(i, lx3);
        sim3.swap_v(i, lv3);
        
        
        sim.swap_r(i, sim.r(i) + 1./6. * 2 * kx3);
        sim.swap_v(i, sim.v(i) + 1./6. * 2 * kv3);
        }
    
    for (int i=0; i < sim.bodies(); i++){

        Vec kx4 = h * sim3.v(i);
        Vec kv4 = h * a(i, sim3);

        sim.swap_r(i, sim.r(i) + 1./6. * kx4);
        sim.swap_v(i, sim.v(i) + 1./6. * kv4);
        }
    
    return sim;
    }


int main(){
    string initial_i;
    nbody sim(2);
    int l = 0;
    fstream initialNfile("Initial_cond.txt");

    while (getline(initialNfile, initial_i)){
        
           

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
        
        };
    initialNfile.close();
    int time = 1;
    int steps = time / h;
    ofstream outfile("RK4N.txt");
    for (int t=0; t < steps; t++){
        
        sim = RK4N_step(h, sim);
        outfile << sim.r(0).x() << ' ' << sim.r(0).y() << '\n';
    }
    
    };
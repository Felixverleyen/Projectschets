#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "nbody2.hpp"
using namespace std;


double G=1;
//double mu=0.01;
double h=1e-3;
//int const N=4;
const double t0 = 0.;
const double t_end = 10.;

  
void print(Vec a)
{ cout << a.x() << ' ' << a.y() << ' '<< a.z() << endl; }

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

// velocity verlet
void verlet(nbody sim) {

    // file
    ofstream outfile("verlet2.txt");
    outfile << setprecision(15); 

    //loop
    int N = sim.bodies();
    Vec vhalf = {0,0,0};
    Vec r = {0,0,0};
    Vec vfull = {0,0,0};
    //vector<Vec> all_r;
    //vector<Vec> all_v; 

    for (double t=t0; t<t_end+h/2.; t+=h) {

        for (int i=0; i<N; i++) {
            vhalf = sim.v(i) + 0.5*h*a(i,sim);
            r = sim.r(i) + h*vhalf;
            vfull = vhalf + 0.5*h*a(i,sim);

            //print
            //outfile << r.x() << ' ' << r.y() << ' ' << r.z() << ' ' ;
           
            //swap waarden in lijst
            sim.swap_r(i,r);
            sim.swap_v(i,vfull);

            outfile << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' ';

        }

        outfile << '\n';

    }

    outfile.close();

}

int main(){

    string initial_i;
    nbody sim;
    int N = 0;
    int l = 0;
    fstream initialNfile("Initial_cond.txt");

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

    
    verlet(sim);

}


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

void RKN4(double h, double time, nbody sim){
    ofstream outfile1("RK4N.txt");
    ofstream outfile2("RK4NLE.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15);
    int steps = int(time / h);

    for (int t=0; t < steps; t++){
        sim = RK4N_step(h, sim);
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        double E = Energy(sim);
        outfile2 << E << '\n';
    }
}


int main(){
    string file = "Initial_cond.txt";
    nbody sim = init_sim(file);
    double h= 1e-4;
    int time = 10;
    RKN4(h, time, sim);
};
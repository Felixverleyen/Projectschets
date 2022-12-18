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
    ofstream outfile2("VerletLE.txt");
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
    string file = "Initial_cond.txt";
    nbody sim = init_sim(file);
    double h = 1e-5;
    double time = 10;
    verlet(h, time, sim);
}



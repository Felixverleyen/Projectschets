#define _USE_MATH_DEFINES
#include "nbody.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
using namespace std;


nbody RK4_step(double h, nbody sim){
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

nbody AB_step(double h, nbody sim_AB, nbody sim0, nbody sim1, nbody sim2, nbody sim3){
    for(int i = 0; i < sim_AB.bodies(); i++){
        sim_AB.swap_v(i, sim3.v(i) + h/24. * (55*a(i, sim3) - 59*a(i, sim2) + 37*a(i, sim1) - 9*a(i, sim0)));
        sim_AB.swap_r(i, sim3.r(i) + h/24. * (55*sim3.v(i) - 59*sim2.v(i) + 37*sim1.v(i) - 9*sim0.v(i)));
    }
    return sim_AB;
}

nbody AM_step(double h, nbody sim, nbody sim1, nbody sim2, nbody sim3){
    for(int i = 0; i < sim.bodies(); i++){
        sim.swap_v(i, sim3.v(i) + h/24. * (9*a(i, sim) + 19*a(i, sim3) - 5*a(i, sim2) + 1*a(i, sim1)));
        sim.swap_r(i, sim3.r(i) + h/24. * (9*sim.v(i) + 19*sim3.v(i) - 5*sim2.v(i) + 1*sim1.v(i)));
    }
    return sim;
}

void AM(double h, double time, nbody sim){
    ofstream outfile1("Adams-Moulton.txt");
    ofstream outfile2("Adams-MoultonLE.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15);
    int steps = int(time / h);

    //We determine the first 4 positions with the Runge_Kutta 4 method
    array<nbody,4> sim_list;
    sim_list[0] = sim;
    for(int i=0; i<sim.bodies(); i++){
        outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
    }
    outfile1 << '\n';
    double E = Energy(sim);
    outfile2 << E << '\n';
    for(int i = 1; i < 4; i++){
        sim = RK4_step(h, sim);
        sim_list[i] = sim;
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        E = Energy(sim);
        outfile2 << E << '\n';
    }

    //The Adams-Moulton method is implicit, we first have to make a prediction of the position with the Adams-Bashford method
    //We then use the Adams-Moulton method to correct this prediction 
    for (int i=4; i<steps; i++){
        sim = AB_step(h, sim, sim_list[0], sim_list[1], sim_list[2], sim_list[3]);
        sim = AM_step(h, sim, sim_list[1], sim_list[2], sim_list[3]);
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        E = Energy(sim);
        outfile2 << E << '\n';
        sim_list[0]=sim_list[1];
        sim_list[1]=sim_list[2];
        sim_list[2]=sim_list[3];
        sim_list[3]=sim;
    }
    outfile1.close();
    outfile2.close();
}


int main(){
    string file = "Initial_cond.txt";
    nbody sim = init_sim(file);
    double h = 1e-4;
    double time = 10;
    AM(h, time, sim);
}

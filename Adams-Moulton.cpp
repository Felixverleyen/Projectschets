#define _USE_MATH_DEFINES
#include "nbody.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
using namespace std;

double G=8;

void print(Vec a){ 
    cout << a.x() << ' ' << a.y() << ' '<< a.z() << endl; }


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


nbody RK4_step(double h, nbody sim){
    
    return sim;
}

nbody AB_step(double h, nbody sim_AB, nbody sim0, nbody sim1, nbody sim2, nbody sim3){
    for(int i = 0; i < sim_AB.bodies(); i++){
        sim_AB.v(i) = sim3.v(i) + h/24. * (55*a(i, sim3) - 59*a(i, sim2) + 37*a(i, sim1) - 9*a(i, sim0));
        sim_AB.r(i) = sim3.r(i) + h/24. * (55*sim3.v(i) - 59*sim2.v(i) + 37*sim1.v(i) - 9*sim0.v(i));
    }
    return sim_AB;
}

nbody AM_step(double h, nbody sim_AM, nbody sim_AB, nbody sim1, nbody sim2, nbody sim3){
    for(int i = 0; i < sim_AB.bodies(); i++){
        sim_AM.v(i) = sim3.v(i) + h/24. * (9*a(i, sim_AM) + 19*a(i, sim3) - 5*a(i, sim2) + a(i, sim1));
        sim_AM.r(i) = sim3.r(i) + h/24. * (9*sim_AM.v(i) + 19*sim3.v(i) - 5*sim2.v(i) + sim1.v(i));
    }
    return sim_AM;
}

void AM(double h, double time, nbody sim){
    ofstream outfile("Adams-Moulton.txt");
    outfile << setprecision(15);
    int steps = int(time / h);

    //We determine the first 4 positions with the Runge_Kutta 4 method
    array<nbody,4> sim_list;
    sim_list[0] = sim;
    for(int i=0; i<sim.bodies(); i++){
        outfile << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
    }
    outfile << '\n';
    for (int i = 1; i < 4; i++){
        sim = RK4_step(h, sim);
        sim_list[i] = sim;
        for(int i=0; i<sim.bodies(); i++){
        outfile << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile << '\n';
    }

    //The Adams-Moulton method is implicit, we first have to make a prediction of the position with the Adams-Bashford method
    //We then use the Adams-Moulton method to correct this prediction 
    for (int i=4; i<steps; i++){
        nbody sim_AB = AB_step(h, sim, sim_list[0], sim_list[1], sim_list[2], sim_list[3]);
        nbody sim_AM = AM_step(h, sim, sim_AB, sim_list[1], sim_list[2], sim_list[3]);
       for(int i=0; i<sim.bodies(); i++){
            outfile << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile << '\n';
        sim_list[0]=sim_list[1];
        sim_list[1]=sim_list[2];
        sim_list[2]=sim_list[3];
        sim_list[3]=sim_AM;
    }
    outfile.close();
}


int main(){
    double h=0.01;
    double time = 500;
    string initial_i;
    nbody sim;
    int N = 0;
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
            
            ++N;
        };

    sim.set_N(N);
    initialNfile.close();
    AM(h, time, sim);
}

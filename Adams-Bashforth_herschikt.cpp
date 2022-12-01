#define _USE_MATH_DEFINES
#include "nbody.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
using namespace std;

double G = 1;

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
    return ar;
}

double Energy(nbody sim){
    double E=0;

    for (int i=0; i<sim.bodies();++i){
        Vec ri = sim.r(i);
        double mi=sim.m(i);
        Vec vi= sim.v(i);
        E+=mi*vi.norm2()/2;

        for(int j =0; j<sim.bodies(); ++j){
            Vec rj = sim.r(j);
            double mj= sim.m(j);
            Vec afst= ri-rj;

            if (i!=j){
                E-= 1/2*G*mi*mj/afst.norm();
            }
        }
    }
    return E;
}

nbody AB_1step(double h, nbody n, nbody n0){
    for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i,n0.r(i)+h*n0.v(i));
        n.swap_v(i, n0.v(i)+h*a(i,n0));
    }
    return n;
}

nbody AB_2step(double h, nbody n, nbody n0, nbody n1){
    for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i, n1.r(i) +h*(3./2.*n1.v(i)-0.5*n0.v(i)));
        n.swap_v(i, n1.v(i) +h*(3./2.*a(i,n1)-0.5*a(i,n0)));
    }
    return n;
}

nbody AB_3step(double h, nbody n, nbody n0, nbody n1, nbody n2){
     for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i,n2.r(i)+h/12.*(23*n2.v(i)-16*n1.v(i)+5*n0.v(i)));
        n.swap_v(i,n2.v(i)+h/12.*(23*a(i,n2)-16*a(i,n1)+5*a(i,n0)));
    }
    return n;
}

nbody AB_4step(double h, nbody n, nbody n0, nbody n1, nbody n2, nbody n3){
    for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i,n3.r(i) +h/24.*(55*n3.v(i)-59*n2.v(i)+37*n1.v(i)-9*n0.v(i)));
        n.swap_v(i,n3.v(i) +h/24.*(55*a(i,n3)-59*a(i,n2)+37*a(i,n1)-9*a(i,n0)));
    }
    return n;
}

void AB(double h, double time, nbody sim, int inputorder){
    ofstream outfile1("Adams-Bashforth.txt");
    ofstream outfile2("Adams-BashforthE.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15);
    int steps = int(time / h);
    double E = Energy(sim);

    array<nbody,4> sim_list;
    sim_list[0] = sim;

    sim = AB_1step(h, sim, sim_list[0]);
    sim_list[1] = sim;
    
    sim = AB_2step(h, sim, sim_list[0], sim_list[1]);
    sim_list[2] = sim;

    sim = AB_3step(h, sim, sim_list[0], sim_list[1], sim_list[2]);
    sim_list[3] = sim;


    if(inputorder==1){
        for (int i=1; i<steps; i++){
            sim = AB_1step(h, sim, sim_list[0]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            sim_list[0] = sim;
        }
    }

    if(inputorder==2){
        for (int i=1; i<steps; i++){
            sim = AB_2step(h, sim, sim_list[0], sim_list[1]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            sim_list[0] = sim_list[1];
            sim_list[1] = sim;
        }
    }

    if(inputorder==3){
        for (int i=1; i<steps; i++){
            sim = AB_3step(h, sim, sim_list[0], sim_list[1], sim_list[2]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            sim_list[0] = sim_list[1];
            sim_list[1] = sim_list[2];
            sim_list[2] = sim;
        }
    }

    if(inputorder==4){
        for (int i=1; i<steps; i++){
            sim = AB_4step(h, sim, sim_list[0], sim_list[1], sim_list[2], sim_list[3]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            sim_list[0] = sim_list[1];
            sim_list[1] = sim_list[2];
            sim_list[2] = sim_list[3];
            sim_list[3] = sim;
        }
    }
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
    double h = 1e-7;
    double time = 0.02;
    int inputorder = 4;
    AB(h, time, sim, inputorder);
}

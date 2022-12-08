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

double theta= 1.351207;

nbody FR_step(double h, nbody n0){
    for(int i = 0; i < n0.bodies(); i++){
        n0.swap_r(i, n0.r(i)+h*theta*n0.v(i)*0.5);
        n0.swap_v(i, n0.v(i)+theta*h*a(i,n0));
        n0.swap_r(i, n0.r(i)+h*(1-theta)*n0.v(i)*0.5);
        n0.swap_v(i, n0.v(i)+h*(1-2*theta)*a(i,n0));
        n0.swap_r(i, n0.r(i)+ h*(1-theta)*n0.v(i)*0.5); 
        n0.swap_v(i, n0.v(i)+theta*h*a(i,n0));
        n0.swap_r(i, n0.r(i)+h*theta*n0.v(i)*0.5);
    }
    return n0;
}

void FR(double h, double time, nbody sim){
    ofstream outfile1("FR.txt");
    ofstream outfile2("FRLE.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15); 
    int steps = int(time / h);
    
    for (double t=0; t<steps; t++){
        sim = FR_step(h, sim);
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
    FR(h, time, sim);
}

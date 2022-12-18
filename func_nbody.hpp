#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include "class_nbody.hpp"
using namespace std;

Vec a(int i, nbody sim) {
    drivercount.countplus();
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
                E-= (1./2.*G*mi*mj/afst.norm());
                }
            }
    }

     return E;}

nbody init_sim(string file){
    string initial_i;
    nbody sim;
    int N = 0;
    int l = 0;
    fstream initialNfile(file);

    while (getline(initialNfile, initial_i)){
        if(l>7){
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
    return sim;
}

// we choose as an adaptive time scheme a function that changes h according to the minimum distance between 2 bodies
// s = r / scale(for example = 1 AU) => h'=f(s)*h with f(s)= s^power (=> f(1)=1), we can vary the power to give different adaptive schemes
// As a power function rises fast for f(s) > 1 we choose the maximum of f(s)=10 and constant for any value of s that gives s^power > 10

double time_step_scale(nbody sim,  double scale, double power){
    //first the minimum distance between any 2 bodies is calculated

    vector<double> diff;
    //put all distances between 2 bodies in a vector
    for (int i = 0; i < sim.bodies(); i++){
        for (int j = 0; j < sim.bodies(); j++){

            if (j != i){
                Vec d = sim.r(i) - sim.r(j);
                diff.push_back(d.norm());
            }

            else { continue;}
        
        }
    }
    
    //find the minimum in the vector
    double min = diff[0];
    
    for (int i = 0; i < diff.size(); i++){

        if(diff[i] < min){
            min = diff[i];
        }
    }
    
    //last step returns f(s)
    double value = pow(min / scale, power);
    if(value < 10) {return value;}
    else {return 10;};
}





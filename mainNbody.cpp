#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <string.h>
#include "int_nbody.hpp"
using namespace std;



int main(){
    string file = "Initial_cond.txt";
    nbody sim = init_sim(file);
    double h = 1e-4;
    double time = 0.01;
    bool adapt = false;
    //see func_nbody.hpp for the source-code of the integrators

    //RK4integrator(h, time, sim, adapt, 1., 2.);
    //verlet(h, time, sim);
    //FR(h, time, sim);
    //AM(h, time, sim);
    //AB(h, time, sim, 2);

    bool hlooping = true;
    bool vary_adapt_looping = false;

    if (hlooping){
        hloop(h, time, sim);    //see bottom of int_nbody.hpp
    };
    if (vary_adapt_looping){
        adapt_loop(h, time, 1., sim);   //see bottom of int_nbody.hpp
    };
    }
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <string.h>
#include "loops_nbody.hpp"
using namespace std;



int main(){
    string file = "Initial_cond.txt";
    nbody sim = init_sim(file);         //reads "Initial_cond.txt" and initializes sim accordingly (see func_nbody.hpp)
                                        //for class nbody see class_nbody.hpp
                                        
    double h = 1e-4;                    //timestep
    double time = 0.1;                  //integration time

    bool adapt = false;     //=true if you want an adaptive time step scheme

    //see func_nbody.hpp for the source-code of the integrators
    //remove the slashes of the integrators you want to test:

    //RK4integrator(h, time, sim, adapt, 1., 2.5);      //last 2 numbers are parameters of the adaptive scheme (see func_nbody line 89)
    //verlet(h, time, sim);
    //FR(h, time, sim);
    //AM(h, time, sim);
    //AB(h, time, sim, 2);


    bool hlooping = false;          //=true if you want to generate output for comparing the rel. E. error for different h
    bool vary_adapt_looping = false; //=true if you want to generate output for comparing the rel. E. error for different parameters...
                                     //... of the adaptive time step scheme

    if (hlooping){
        hloop(h, time, sim);    //see loops_nbody.hpp
    };
    if (vary_adapt_looping){
        adapt_loop(h, time, 1., sim);   //see loops_nbody.hpp
    };
    
    }
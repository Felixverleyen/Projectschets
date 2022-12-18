#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <string.h>
#include <chrono>
#include "int_nbody.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///special functions for looping over certain parameters:

// function for looping over different h and writing the final relative energy error to a file

void hloop(double h , double time, nbody sim){
    ofstream outfile("Adams-BashforthE(h).txt");    //writingfile has name [integrator]E(h).txt
    outfile << setprecision(15);
    ofstream outfile2("Adams-MoultonE(h).txt");
    outfile2 << setprecision(15);
    ofstream outfile3("Foresth-RuthE(h).txt");
    outfile3 << setprecision(15);
    ofstream outfile4("RK4NE(h).txt");
    outfile4 << setprecision(15);
    ofstream outfile5("VerletE(h).txt");
    outfile5 << setprecision(15);

    //We will loop over i with h changing as h_i = 10e-i * h
    for (int i=0; i<2; i++){
        
        double u = pow(10, -i) * h; // h_i = 10e-i * h
        int steps = time / u;

        //simulate for h_i
        RK4integrator(u, time, sim, false, 1., 1.);
        verlet(u, time, sim);
        FR(u, time, sim);
        AM(u, time, sim);
        AB(u, time, sim, 4);

        //now for h_i the energies been written to the [integrator]E.txt files, we need the first (E_0) and last (E_f) to calculate...
        //... the relative energy error. this is done for each integrator in the same way:

        string En;
        ifstream E("Adams-BashforthE.txt");
    
        int l = 0;
        double E0;
        double Ef;
        while (getline(E, En)){
            if (l==0){
                E0 = stod(En);
            };

            if (l==steps-1){
                Ef = stod(En);
            };
            ++l;
            
        }  
        outfile << abs((Ef-E0)/E0) << ' ' << u << '\n';
        E.close();

        ifstream E2("Adams-MoultonE.txt");

        int l2 = 0;
        double E02;
        double Ef2;
        while (getline(E2, En)){
            if (l2==0){
                E02 = stod(En);
            };

            if (l2==steps-1){
                Ef2 = stod(En);
            };
            ++l2;
            
        }  
        outfile2 << abs((Ef2-E02)/E02) << ' ' << u << '\n';

        E2.close();

        ifstream E3("FRE.txt");

        int l3 = 0;
        double E03;
        double Ef3;
        while (getline(E3, En)){
            if (l3==0){
                E03 = stod(En);
            };

            if (l3==steps-2){
                Ef3 = stod(En);
            };
            ++l3;
            
        }  
        outfile3 << abs((Ef3-E03)/E03) << ' ' << u << '\n';
        E3.close();

        ifstream E4("RK4NE.txt");
        
        int l4 = 0;
        double E04;
        double Ef4;
        while (getline(E4, En)){
            if (l4==0){
                E04 = stod(En);
            };

            if (l4==steps-1){
                Ef4 = stod(En);
            };
            ++l4;
            
        }  
        outfile4 << abs((Ef4-E04)/E04) << ' ' << u << '\n';
        E4.close();

        ifstream E5("VerletE.txt");
        int l5 = 0;
        double E05;
        double Ef5;
        while (getline(E5, En)){
            if (l5==0){
                E05 = stod(En);
            };

            if (l5==steps-1){
                Ef5 = stod(En);
            };
            ++l5;
            
        }  
        outfile5 << abs((Ef5-E05)/E05) << ' ' << u << '\n';
        E5.close();
    }}


// function for looping over different parameters of the adaptive timestep scheme and writing the final relative energy error to a file
// We will change the power of the adaptive scheme from 1 to 3 in steps of 0.25

void adapt_loop(double h, double time, double scale, nbody sim){
    //initialize file; output is 2 columns, column 1 = rel. energy error at indicated time, column 2 = power in the time adaptive scheme
    ofstream outfile("RK4N_adapt_loop_E.txt");
    outfile << setprecision(15);

    //Different powers in the adaptive scheme causes that the simulations will reach the same physical time in a different amount... 
    //... of steps. As we want to compare energies at the same physical time we first look how far the simulation with the highest...
    //... power reaches in pysical time (higher power means more smaller steps). We then find the energy of the other simulations...
    //... at that time.
    //note that the argument (higher power means more smaller steps) only counts if the scale is the maximum distance that is reached...
    //... during the simulation. We run this loop for the initial conditions of the earth but with less orbital speed (more eccentric...
    //... orbit: yvel = 400. instead of = 581.

    // for power = 3:
    RK4integrator(h, time, sim, true, 1., 3.);
    double maxt;
    int w = 0;

    string filename = "time.txt";
    ifstream fin;
    fin.open(filename);
    string tn;

    while(getline(fin, tn)){        //this loop makes the last value in the file equal to maxt and gives the length of the file w
        maxt = stod(tn);
        ++w;
    }

    fin.close();
      
    // we calculate the relative energy error at maxt 
    string En;
    ifstream E("RK4NE.txt");
    
    int l = 0;
    double E0;
    double Ef;
    while (getline(E, En)){
        if (l==0){
            E0 = stod(En);
        };

        if (l==w-2){
            Ef = stod(En);
        };
            
        ++l;
            
    }  
    outfile << abs((Ef-E0)/E0) << ' ' << 3. << ' ' << maxt << '\n'; //in the first line of the file we also specify the physical time...
                                                                    //... we consider when comparing the energies
    E.close();         
    
    // Now for other powers
    double min = maxt - 1e-5;       //we must make an interval because the steps in physical time don't correspond perfectly for...
    double max = maxt + 1e-5;       //... different simulations with different powers in the adaptive scheme (see further)

    for (int i = 1; i < 9; i++){
        double power = (12-i) / 4.;
        RK4integrator(h, time, sim, true, 1., power);

        string tn_;
        ifstream t2("time.txt");
        int u = 0;

        while (getline(t2, tn_)){   //this loop finds how many steps u the simulation needs to reach maxt
            ++u;
            double tn_d = stod(tn_);

            if (min < tn_d &&  tn_d < max ){ //the considered interval
                break;
            };
        };

        //the energy at step u will be at physical time t
        string En;
        ifstream E("RK4NE.txt"); 
    
        int l = 0;
        double E0;
        double Ef;
        //we calculate the relative energy error
        while (getline(E, En)){     
            if (l==0){
                E0 = stod(En);
            };

            if (l==u){
                Ef = stod(En);
            };
            
            ++l;     
        }  
        outfile << abs((Ef-E0)/E0) << ' ' << power << '\n';
        E.close(); 
        
    }
}

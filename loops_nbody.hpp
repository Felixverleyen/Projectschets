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
    ofstream outfile("Adams-BashforthE(h).txt");    //file has name [integrator]E(h).txt
    outfile << setprecision(15);
    ofstream outfile2("Adams-MoultonE(h).txt");
    outfile2 << setprecision(15);
    ofstream outfile3("Foresth-RuthE(h).txt");
    outfile3 << setprecision(15);
    ofstream outfile4("RK4NE(h).txt");
    outfile4 << setprecision(15);
    ofstream outfile5("VerletE(h).txt");
    outfile5 << setprecision(15);

    //We will loop over i with h changing to h_i=10e-i * h
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
        //string En3;
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
        //string En4;
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

void adapt_loop(double h, double time, double scale, nbody sim){
    ofstream outfile("RK4N_adapt_loop_E.txt");
    outfile << setprecision(15);

    RK4integrator(h, time, sim, true, 1., 3.);
    double maxt;
    int w = 0;

    string filename = "time.txt";
    ifstream fin;
    fin.open(filename);
    string tn;

    while(getline(fin, tn)){
        maxt = stod(tn);
        ++w;
    }

    fin.close();
      
        
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
    outfile << abs((Ef-E0)/E0) << ' ' << 3. << '\n';
    E.close();         
    
    double min = maxt - 1e-5;
    double max = maxt + 1e-5;

    for (int i = 1; i < 9; i++){
        double power = (12-i) / 4.;
        RK4integrator(h, time, sim, true, 1., power);
        string tn_;
        ifstream t2("time.txt");
        int u = 0;
        while (getline(t2, tn_)){
            ++u;
            double tn_d = stod(tn_);
            if (min < tn_d &&  tn_d < max ){
                cout << tn_d;
                break;
            };
        };

        string En;
        ifstream E("RK4NE.txt");
    
        int l = 0;
        double E0;
        double Ef;
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
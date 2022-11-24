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


double G=8;
double mu=0.01;
double h=0.01;
int const N=4;


int main(){

};



   
void print(Vec a)
{ cout << a.x() << ' ' << a.y() << ' '<< a.z() << endl; }

Vec a(int i, nbody sim) {
    
    Vec ri= sim.p(i);
    Vec ar={0,0,0};

    for (int j =0; j<N; ++j){
        Vec rj = sim.p(j);
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

    for (int i=0; i<N;++i){
        Vec ri = sim.p(i);
        double mi=sim.m(i);
        Vec vi= sim.v(i);
        E+=mi*vi.norm2()/2;

        for  (int j =0; j<N; ++j){
            Vec rj = sim.p(j);
            double mj= sim.m(j);
            Vec afst= ri-rj;

            if (i!=j){
                E-= 1/2*G*mi*mj/afst.norm();
                }
            }
    }

     return E;}

 // integratoren:

    //Runge-Kutta

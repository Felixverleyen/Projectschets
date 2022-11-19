#define _USE_MATH_DEFINES
#include "3Vec.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
using namespace std;


Vec a(Vec r, double m){
    return -1 * m * (r / r.norm3());
}

array<Vec,2> RK4_step(double h, Vec ri, Vec vi, double mu){
    Vec kv1 = h *a(ri, mu);
    Vec kx1 = h * vi;

    Vec lv = ri + 0.5 * kv1;
    Vec lx = vi + 0.5 * kx1;

    Vec kv2 = h *a(lv, mu);
    Vec kx2 = h * lx;

    lv = ri + 0.5 * kv2;
    lx = vi + 0.5 * kx2;

    Vec kv3 = h *a(lv, mu);
    Vec kx3 = h * lx;

    lv = ri + kv3;
    lx = vi + kx3;

    Vec kv4 = h *a(lv, mu);
    Vec kx4 = h * lx;

    Vec vo = vi + (1./6.) * (kv1 + 2 * kv2 + 2 * kv3 + kv4);
    Vec ro = ri + (1./6.) * (kx1 + 2 * kx2 + 2 * kx3 + kx4);
    array<Vec,2> opl = {ro, vo};

    return opl;
}

array<Vec,2> AB_step(double h, array<Vec,2> opl0, array<Vec,2> opl1, array<Vec,2> opl2, array<Vec,2> opl3, double m){
    Vec v_AB = opl3[1] + h/24. * (55*a(opl3[0], m) - 59*a(opl2[0], m) + 37*a(opl1[0], m) - 9*a(opl0[0], m));
    Vec r_AB = opl3[0] + h/24. * (55*opl3[1] - 59*opl2[1] + 37*opl1[1] - 9*opl0[1]);
    array<Vec,2> opl_AB = {r_AB, v_AB};
    return opl_AB;
}

array<Vec,2> AM_step(double h, array<Vec,2> opl1, array<Vec,2> opl2, array<Vec,2> opl3, array<Vec,2> opl_AB, double m){
    Vec v_AM = opl3[1] + h/24. * (9*a(opl_AB[0], m) + 19*a(opl3[0], m) - 5*a(opl2[0], m) + a(opl1[0], m));
    Vec r_AM = opl3[0] + h/24. * (9*opl3[1] + 19*opl2[1] - 5*opl1[1] + opl1[1]);
    array<Vec,2> opl_AM = {r_AM, v_AM};
    return opl_AM;
}

void AM(double h, double mu, double time, Vec r0, Vec v0){
    ofstream outfile("Adams-Moulton.txt");
    outfile << setprecision(15);
    int steps = int(time / h);
    array<Vec,2> opli[steps];

    //We determin the first 4 positions with the Runge_Kutta 4 method
    array<Vec,2>  opl = RK4_step(h, r0, v0, mu);
    Vec r = opl[0];
    Vec v = opl[1];
    opli[0] = opl;
    outfile << r.x() << ' ' << r.y() << ' ' << r.z() << '\n'; 
    outfile << setprecision(15);
    for (int i = 1; i < 4; i++){
        array<Vec,2>  opl = RK4_step(h, r, v, mu);
        r = opl[0];
        v = opl[1];
        opli[i] = opl;
        outfile << r.x() << ' ' << r.y() << ' ' << r.z() << '\n'; 
    }

    //The Adams-Moulton method is implicit, we first have to make a prediction of the position with the Adams-Bashford method
    //We then use the Adams-Moulton method to correct this prediction 
    for (int i=4; i<steps; i++){
        array<Vec,2>  opl_AB = AB_step(h, opli[0], opli[1], opli[2], opli[3], mu);
        array<Vec,2>  opl_AM = AM_step(h, opli[1], opli[2], opli[3], opl_AB, mu);
        outfile << opl_AM[0].x() << ' ' << opl_AM[0].y() << ' ' << opl_AM[0].z() << '\n'; 
        opli[0]=opli[1];
        opli[1]=opli[2];
        opli[2]=opli[3];
        opli[3]=opl_AM;
    }
}


int main(){
    Vec r0(1.2, 0.5, 0.1);
    Vec v0(0., 0.5, 0.);
    
    double mu = 1;
    double h = 0.001;
    double t = 50;
    AM(h, mu, t, r0, v0);
}

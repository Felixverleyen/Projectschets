#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <string.h>
#include "nbodyf.hpp"
using namespace std;


double G=1;
counter drivercount;

///// integratoren:

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Runge-Kutta4

nbody RK4N_step(double h, nbody sim){

    nbody sim1 = sim;
    nbody sim2 = sim;
    nbody sim3 = sim;

    for (int i=0; i < sim.bodies(); i++){
        

        Vec kx1 = h * sim.v(i);
        Vec kv1 = h * a(i, sim);

        Vec lx1 = sim.r(i) + 0.5 * kx1;
        Vec lv1 = sim.v(i) + 0.5 * kv1;
        sim1.swap_r(i, lx1);
        sim1.swap_v(i, lv1);
        
        
        sim.swap_r(i, sim.r(i) + 1./6. * kx1);
        sim.swap_v(i, sim.v(i) + 1./6. * kv1);
        }
    

    for (int i=0; i < sim.bodies(); i++){
        

        Vec kx2 = h * sim1.v(i);
        Vec kv2 = h * a(i, sim1);

        Vec lx2 = sim1.r(i) + 0.5 * kx2;
        Vec lv2 = sim.v(i) + 0.5 * kv2;
        sim2.swap_r(i, lx2);
        sim2.swap_v(i, lv2);
        
        sim.swap_r(i, sim.r(i) + 1./6. * 2 *  kx2);
        sim.swap_v(i, sim.v(i) + 1./6. * 2 * kv2);
        }
    
    for (int i=0; i < sim.bodies(); i++){
        

        Vec kx3 = h * sim2.v(i);
        Vec kv3 = h * a(i, sim2);

        Vec lx3 = sim2.r(i) + kx3;
        Vec lv3 = sim2.v(i) + kv3;
        sim3.swap_r(i, lx3);
        sim3.swap_v(i, lv3);
        
        
        sim.swap_r(i, sim.r(i) + 1./6. * 2 * kx3);
        sim.swap_v(i, sim.v(i) + 1./6. * 2 * kv3);
        }
    
    for (int i=0; i < sim.bodies(); i++){


        Vec kx4 = h * sim3.v(i);
        Vec kv4 = h * a(i, sim3);

        sim.swap_r(i, sim.r(i) + 1./6. * kx4);
        sim.swap_v(i, sim.v(i) + 1./6. * kv4);
        }
    
    return sim;
    };

void RK4integrator(double u, double time, nbody sim){
    
    ofstream outfile1("RK4N.txt");
    ofstream outfile2("RK4N_E.txt");
    bool stop = false;
    int steps = int(time/u);

    for (int t=0; t < steps; t++){

        double scale = time_step_scale(sim, 0, 1);
        double h = scale * u;
        sim = RK4N_step(h, sim);       
        

        for(int i=0; i<sim.bodies(); i++){
            if (sim.r(i).norm() > 5){
                stop = true;
            }
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            if (stop){
                break;
            }
        
        double E = Energy(sim);
        outfile2 << E << '\n';
        
        }
    outfile1.close();
    outfile2.close();
    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << averaged_count << '\n';
    drivercount.reset();}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//velocity Verlet

 nbody verlet_halfstep(double h, nbody sim_half){
    for (int i=0; i<sim_half.bodies(); i++){
        sim_half.swap_v(i, sim_half.v(i) + 0.5*h*a(i,sim_half));
    }
    return sim_half;
 }

 nbody verlet_step(double h, nbody sim, nbody sim_half){
    for (int i=0; i<sim.bodies(); i++) {
        sim.swap_r(i, sim.r(i) + h*sim_half.v(i));
        sim.swap_v(i, sim_half.v(i) + 0.5*h*a(i,sim));
    }
    return sim;
 }


void verlet(double h, double time, nbody sim) {
    ofstream outfile1("verlet.txt");
    ofstream outfile2("VerletE.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15); 
    int steps = int(time / h);
    
    for (double t=0; t<steps; t++) {
        nbody sim_half = verlet_halfstep(h, sim);
        sim = verlet_step(h, sim, sim_half);
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        double E = Energy(sim);
        outfile2 << E << '\n';
    }

    outfile1.close();
    outfile2.close();

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << averaged_count << '\n';
    drivercount.reset();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Forest-Ruth

double theta= 1.351207;

nbody FR_step(double h, nbody sim){
    for(int i = 0; i < n0.bodies(); i++){
        sim.swap_r(i, sim.r(i)+h*theta*sim.v(i)*0.5);
        sim.swap_v(i, sim.v(i)+theta*h*a(i,sim));
        sim.swap_r(i, sim.r(i)+h*(1-theta)*sim.v(i)*0.5);
        sim.swap_v(i, sim.v(i)+h*(1-2*theta)*a(i,sim));
        sim.swap_r(i, sim.r(i)+ h*(1-theta)*sim.v(i)*0.5); 
        sim.swap_v(i, sim.v(i)+theta*h*a(i,sim));
        sim.swap_r(i, sim.r(i)+h*theta*sim.v(i)*0.5);
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

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << averaged_count << '\n';
    drivercount.reset();
}
        
    


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Adams-Bashfort

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
    outfile1.close();
    outfile2.close();

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << averaged_count << '\n';
    drivercount.reset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Adams-Moulton

nbody AM_step(double h, nbody sim, nbody sim1, nbody sim2, nbody sim3){
    
    for(int i = 0; i < sim.bodies(); i++){
        sim.swap_v(i, sim3.v(i) + h/24. * (9*a(i, sim) + 19*a(i, sim3) - 5*a(i, sim2) + 1*a(i, sim1)));
        sim.swap_r(i, sim3.r(i) + h/24. * (9*sim.v(i) + 19*sim3.v(i) - 5*sim2.v(i) + 1*sim1.v(i)));
    }
    return sim;
}

void AM(double h, double time, nbody sim){
    ofstream outfile1("Adams-Moulton.txt");
    ofstream outfile2("Adams-MoultonE.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15);
    int steps = int(time / h);

    //We determine the first 4 positions with the Runge_Kutta 4 method
    array<nbody,4> sim_list;
    sim_list[0] = sim;
    for(int i=0; i<sim.bodies(); i++){
        outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
    }
    outfile1 << '\n';
    double E = Energy(sim);
    outfile2 << E << '\n';
    for(int i = 1; i < 4; i++){
        //double h = time_step_scale(sim, 0, 0.5) * u;
        sim = RK4N_step(h, sim);
        sim_list[i] = sim;
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        E = Energy(sim);
        outfile2 << E << '\n';
    }

    //The Adams-Moulton method is implicit, we first have to make a prediction of the position with the Adams-Bashford method
    //We then use the Adams-Moulton method to correct this prediction 
    for (int i=4; i<steps; i++){
        sim = AB_4step(h, sim, sim_list[0], sim_list[1], sim_list[2], sim_list[3]);
        sim = AM_step(h, sim, sim_list[1], sim_list[2], sim_list[3]);
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        E = Energy(sim);
        outfile2 << E << '\n';
        sim_list[0]=sim_list[1];
        sim_list[1]=sim_list[2];
        sim_list[2]=sim_list[3];
        sim_list[3]=sim;
    }
    outfile1.close();
    outfile2.close();

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << averaged_count << '\n';
    drivercount.reset();
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//main:

int main(){
    string file = "Initial_cond.txt";
    nbody sim = init_sim(file);
    double h = 1e-5;
    double time = 10;
    
    //RK4integrator(h, time, sim);
    //verlet(h, time, sim);
    //FR(h, time, sim);
    //foresth-ruth(h)
    //AM(h, time, sim);
    //AB(h, time, sim, 3);
    
};
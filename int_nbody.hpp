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
#include "func_nbody.hpp"
using namespace std;
using namespace std::chrono;

///// integratoren:

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Runge-Kutta
nbody RK4N_step(double h, nbody sim){

    nbody sim1 = sim;  //used as input for driverfunction in eq. 3.36 of the syllabus
    nbody sim2 = sim;  //used as input for driverfunction in eq. 3.37 of the syllabus
    nbody sim3 = sim;  //used as input for driverfunction in eq. 3.38 of the syllabus
    nbody sim4 = sim;  //after each calculation of k_i a corresponding term of eq. 3.39 is added. at the end it is sim(t_(n+1))

    for (int i=0; i < sim.bodies(); i++){
        

        Vec kx1 = h * sim.v(i);
        Vec kv1 = h * a(i, sim);

        Vec lx1 = sim.r(i) + 0.5 * kx1;
        Vec lv1 = sim.v(i) + 0.5 * kv1;
        sim1.swap_r(i, lx1);
        sim1.swap_v(i, lv1);
              
        sim4.swap_r(i, sim.r(i) + 1./6. * kx1);
        sim4.swap_v(i, sim.v(i) + 1./6. * kv1);
        }
    

    for (int i=0; i < sim.bodies(); i++){
        

        Vec kx2 = h * sim1.v(i);
        Vec kv2 = h * a(i, sim1);

        Vec lx2 = sim.r(i) + 0.5 * kx2;
        Vec lv2 = sim.v(i) + 0.5 * kv2;
        sim2.swap_r(i, lx2);
        sim2.swap_v(i, lv2);

        sim4.swap_r(i, sim4.r(i) + 1./6. * 2 * kx2);
        sim4.swap_v(i, sim4.v(i) + 1./6. * 2 * kv2);
        }
    
    for (int i=0; i < sim.bodies(); i++){
        

        Vec kx3 = h * sim2.v(i);
        Vec kv3 = h * a(i, sim2);

        Vec lx3 = sim.r(i) + kx3;
        Vec lv3 = sim.v(i) + kv3;
        sim3.swap_r(i, lx3);
        sim3.swap_v(i, lv3);
              
        sim4.swap_r(i, sim4.r(i) + 1./6. * 2 * kx3);
        sim4.swap_v(i, sim4.v(i) + 1./6. * 2 * kv3);
        }
    
    for (int i=0; i < sim.bodies(); i++){


        Vec kx4 = h * sim3.v(i);
        Vec kv4 = h * a(i, sim3);

        sim4.swap_r(i, sim4.r(i) + 1./6. * kx4);
        sim4.swap_v(i, sim4.v(i) + 1./6. * kv4);
        }

    sim = sim4;
    return sim;
    };

void RK4integrator(double u, double time, nbody sim, bool adapt, double scale, double power){
    ofstream outfile1("RK4N.txt");
    ofstream outfile2("RK4NE.txt");
    ofstream outfile3("time.txt");      //needed in adapt_loop function see loops_nbody.hpp
    outfile1 << setprecision(15);
    outfile2 << setprecision(15);
    outfile3 << setprecision(15);
    
    double phystime = 0.;
    int steps = int(time/u);
    double tot = 0.;


    for (int t=0; t < steps; t++){
        
        auto start = high_resolution_clock::now(); // take time at start of integration time step


        if (adapt){ // change time step if adaptive scheme is used
            double scale_ = time_step_scale(sim, scale, power);
            double h = scale_ * u;
            sim = RK4N_step(h, sim);
            phystime += h;          //keeps track of physical time
            outfile3 << phystime << '\n';
        }
        else {sim = RK4N_step(u, sim);}     
        

        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; }
            
        outfile1 << '\n';
            
        
        double E = Energy(sim);
        outfile2 << E << '\n';

        auto stop = high_resolution_clock::now();  // take time at end of integration time step

        auto duration = duration_cast<microseconds>(stop - start); 
        tot += duration.count(); // count total time
        
        }
    outfile1.close();
    outfile2.close();
    
    double meanduration = tot / steps; // mean execution time per integration time step

    cout << "Runge-Kutta 4: " << endl;
    cout << "For N =" << ' ' << sim.bodies() << ' ' << "bodies: The mean execution time per integration time step is" << ' ' << meanduration << ' ' << "microseconds." << endl;
    if (adapt) {cout << "with adaptive time step." << endl;}
    else {cout << "without adaptive time step." << endl;}

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << "The number of driver function evaluations is: " << averaged_count << '\n' << endl;
    drivercount.reset();}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Verlet:
 nbody verlet_halfstep(double h, nbody sim_half){
    nbody sim0 = sim_half;
    for (int i=0; i<sim_half.bodies(); i++){
        sim_half.swap_v(i, sim0.v(i) + 0.5*h*a(i,sim0));
    }
    return sim_half;
 }

 nbody verlet_step(double h, nbody sim, nbody sim_half){
    nbody sim0 = sim;
    for (int i=0; i<sim.bodies(); i++) {
        sim.swap_r(i, sim0.r(i) + h*sim_half.v(i));
        
    }
    sim0 = sim;
    for (int i=0; i<sim.bodies(); i++) {
        sim.swap_v(i, sim_half.v(i) + 0.5*h*a(i,sim0));
    }
    return sim;
 }

// velocity verlet
void verlet(double h, double time, nbody sim) {
    ofstream outfile1("verlet.txt");
    ofstream outfile2("VerletE.txt");
    outfile1 << setprecision(20);
    outfile2 << setprecision(20); 
    int steps = int(time / h);
    double tot =0.;
    
    for (double t=0; t<steps; t++) {
        
        auto start = high_resolution_clock::now();
        
        nbody sim_half = verlet_halfstep(h, sim);
        sim = verlet_step(h, sim, sim_half);
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        
        outfile1 << '\n';
        double E = Energy(sim);
        outfile2 << E << '\n';
        
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(stop - start); 
        tot += duration.count();
    }

    outfile1.close();
    outfile2.close();
    
    double meanduration = tot / steps;
    cout << "Verlet: " << endl;
    cout << "For N =" << ' ' << sim.bodies() << ' ' << "bodies: The mean execution time per integration time step is" << ' ' << meanduration << ' ' << "microseconds." << endl;

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << "The number of driver function evaluations is: " << averaged_count << '\n' << endl;
    drivercount.reset();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Forest-Ruth
double theta= 1.351207;

nbody FR_step(double h, nbody sim){
    nbody sim2 = sim;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_r(i, sim.r(i)+h*theta*sim.v(i)*0.5);
    }
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_v(i, sim.v(i)+theta*h*a(i,sim));}
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_r(i, sim.r(i)+h*(1-theta)*sim.v(i)*0.5);}
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_v(i, sim.v(i)+h*(1-2*theta)*a(i,sim));}
    sim = sim2;
    
    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_r(i, sim.r(i)+ h*(1-theta)*sim.v(i)*0.5);}
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_v(i, sim.v(i)+theta*h*a(i,sim));}
    sim = sim2;
    
    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_r(i, sim.r(i)+h*theta*sim.v(i)*0.5);}
    sim = sim2;
    
    return sim;
}

void FR(double h, double time, nbody sim){
    ofstream outfile1("FR.txt");
    ofstream outfile2("FRE.txt");
    outfile1 << setprecision(20);
    outfile2 << setprecision(20); 
    int steps = int(time / h);
    double tot = 0.;
    
    for (double t=0; t<steps; t++){
        
        auto start = high_resolution_clock::now();
        
        sim = FR_step(h, sim);
        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
        }
        outfile1 << '\n';
        double E = Energy(sim);
        outfile2 << E << '\n';
        
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        tot += duration.count();
    }

    outfile1.close();
    outfile2.close();
    
    double meanduration = tot / steps;
    cout << "Forest-Ruth: " << endl;
    cout << "For N =" << ' ' << sim.bodies() << ' ' << "bodies: The mean execution time per integration time step is" << ' ' << meanduration << ' ' << "microseconds." << endl;

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << "The number of driver function evaluations is: " << averaged_count << '\n' << endl;
    drivercount.reset();
}
        
    


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Adams-Bashfort

nbody AB_1step(double h, nbody n, nbody n0){
    for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i,n0.r(i)+h*n0.v(i));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_v(i, n0.v(i)+h*n0.a(i));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_a(i, a(i, n));
    }
    return n;
}

nbody AB_2step(double h, nbody n, nbody n0, nbody n1){
    for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i, n1.r(i) +h*(3./2.*n1.v(i)-0.5*n0.v(i)));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_v(i, n1.v(i) +h*(3./2.*n1.a(i)-0.5*n0.a(i)));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_a(i, a(i, n));
    }
    return n;
}

nbody AB_3step(double h, nbody n, nbody n0, nbody n1, nbody n2){
    for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i,n2.r(i)+h/12.*(23*n2.v(i)-16*n1.v(i)+5*n0.v(i)));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_v(i,n2.v(i)+h/12.*(23*n2.a(i)-16*n1.a(i)+5*n0.a(i)));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_a(i, a(i, n));
    }
    return n;
}

nbody AB_4step(double h, nbody n, nbody n0, nbody n1, nbody n2, nbody n3){
    for(int i = 0; i < n.bodies(); i++){
        n.swap_r(i,n3.r(i) +h/24.*(55*n3.v(i)-59*n2.v(i)+37*n1.v(i)-9*n0.v(i)));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_v(i,n3.v(i) +h/24.*(55*n3.a(i)-59*n2.a(i)+37*n1.a(i)-9*n0.a(i)));
    }
    for(int i = 0; i < n.bodies(); i++){
        n.swap_a(i, a(i, n));
    }
    return n;
}

void AB(double h, double time, nbody sim, int inputorder){
    ofstream outfile1("Adams-Bashforth.txt");
    ofstream outfile2("Adams-BashforthE.txt");
    outfile1 << setprecision(20);
    outfile2 << setprecision(20);
    int steps = int(time / h);
    double E = Energy(sim);
    double tot = 0.;

    for (int i = 0; i < sim.bodies(); i++){     //here the accelerations have to be added to the nbody datatype
        sim.add_a(a(i, sim));
    }

    array<nbody,4> sim_list;

    sim_list[0] = sim;

    sim = AB_1step(h, sim, sim_list[0]);
    sim_list[1] = sim;
    
    sim = AB_2step(h, sim, sim_list[0], sim_list[1]);
    sim_list[2] = sim;

    sim = AB_3step(h, sim, sim_list[0], sim_list[1], sim_list[2]);
    sim_list[3] = sim;


    if(inputorder==1){
        for (int i=0; i<steps; i++){
            
            auto start = high_resolution_clock::now();
            
            sim = AB_1step(h, sim, sim_list[0]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            outfile2 << E << '\n';
            sim_list[0] = sim;
            
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            tot += duration.count();
        }
    }

    if(inputorder==2){
        for (int i=0; i<steps; i++){
            
            auto start = high_resolution_clock::now();
            
            sim = AB_2step(h, sim, sim_list[0], sim_list[1]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            outfile2 << E << '\n';
            sim_list[0] = sim_list[1];
            sim_list[1] = sim;
            
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            tot += duration.count();
        }
    }

    if(inputorder==3){
        for (int i=0; i<steps; i++){
            
            auto start = high_resolution_clock::now();

            sim = AB_3step(h, sim, sim_list[0], sim_list[1], sim_list[2]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            outfile2 << E << '\n';
            sim_list[0] = sim_list[1];
            sim_list[1] = sim_list[2];
            sim_list[2] = sim;
            
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            tot += duration.count();
        }
    }

    if(inputorder==4){
        for (int i=0; i<steps; i++){
            
            auto start = high_resolution_clock::now();
            
            sim = AB_4step(h, sim, sim_list[0], sim_list[1], sim_list[2], sim_list[3]);
            for(int i=0; i<sim.bodies(); i++){
                outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; 
            }
            outfile1 << '\n';
            E = Energy(sim);
            outfile2 << E << '\n';
            sim_list[0] = sim_list[1];
            sim_list[1] = sim_list[2];
            sim_list[2] = sim_list[3];
            sim_list[3] = sim;
            
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            tot += duration.count();
        }
    }
    outfile1.close();
    outfile2.close();
    
    double meanduration = tot / steps;
    cout << "Adams-Bashforth " << inputorder << ": " << endl;
    cout << "For N =" << ' ' << sim.bodies() << ' ' << "bodies: The mean execution time per integration time step is" << ' ' << meanduration << ' ' << "microseconds." << endl;

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << "The number of driver function evaluations is: " << averaged_count << '\n' << endl;
    drivercount.reset();
}

        
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Adams-Moulton:
nbody AM_step(double h, nbody sim, nbody sim1, nbody sim2, nbody sim3){
    nbody simAB = sim;
    for(int i = 0; i < sim.bodies(); i++){
        sim.swap_v(i, sim3.v(i) + h/24. * (9*simAB.a(i) + 19*sim3.a(i) - 5*sim2.a(i) + 1*sim1.a(i)));
    }
    for(int i = 0; i < sim.bodies(); i++){
        sim.swap_r(i, sim3.r(i) + h/24. * (9*simAB.v(i) + 19*sim3.v(i) - 5*sim2.v(i) + 1*sim1.v(i)));
    }
    for(int i = 0; i < sim.bodies(); i++){
        sim.swap_a(i, a(i, sim));
    }
    return sim;
}

void AM(double h, double time, nbody sim){
    ofstream outfile1("Adams-Moulton.txt");
    ofstream outfile2("Adams-MoultonE.txt");
    outfile1 << setprecision(20);
    outfile2 << setprecision(20);
    int steps = int(time / h);
    double tot = 0.;
    for (int i = 0; i < sim.bodies(); i++){
        sim.add_a(a(i, sim));
    }

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
        sim = RK4N_step(h, sim);
        for (int i = 0; i < sim.bodies(); i++){
            sim.swap_a(i, a(i, sim));
        }
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

    for (int i=0; i<steps; i++){
        
        auto start = high_resolution_clock::now();
        
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
        
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        tot += duration.count();
    }
    outfile1.close();
    outfile2.close();
    
    double meanduration = tot / (steps-4);
    cout << "Adams-Moulton 4: " << endl;
    cout << "For N =" << ' ' << sim.bodies() << ' ' << "bodies: The mean execution time per integration time step is" << ' ' << meanduration << ' ' << "microseconds." << endl;

    double averaged_count = drivercount.order(sim.bodies(), steps);
    cout << "The number of driver function evaluations is: " << averaged_count << '\n' << endl;
    drivercount.reset();
    
}


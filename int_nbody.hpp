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

    nbody sim1 = sim;
    nbody sim2 = sim;
    nbody sim3 = sim;
    nbody sim4 = sim;

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
        
        Vec ux2 = sim4.r(i) + 1./6. * 2 * kx2;
        Vec uv2 = sim4.v(i) + 1./6. * 2 * kv2;
        sim4.swap_r(i, ux2 );
        sim4.swap_v(i, uv2);
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
    ofstream outfile3("time.txt");
    outfile1 << setprecision(15);
    outfile2 << setprecision(15);
    outfile3 << setprecision(15);
    
    double phystime = 0.;
    int steps = int(time/u);
    double tot = 0.; // variable to take total execution time


    for (int t=0; t < steps; t++){
        
        auto start = high_resolution_clock::now(); // take time at start of integration time step


        if (adapt){
            double scale_ = time_step_scale(sim, scale, power);
            double h = scale_ * u;
            sim = RK4N_step(h, sim);
            phystime += h;
            outfile3 << phystime << '\n';
        }
        else {sim = RK4N_step(u, sim);}     
        

        for(int i=0; i<sim.bodies(); i++){
            outfile1 << sim.r(i).x() << ' ' << sim.r(i).y() << ' ' << sim.r(i).z() << ' '; }
            
        outfile1 << '\n';
            
        
        double E = Energy(sim);
        outfile2 << E << '\n';

        auto stop = high_resolution_clock::now(); // take time at end of integration time step

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

    for (int i = 0; i < sim.bodies(); i++){
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///special case functions:

// function for looping over different h and writing the final relative energy error to a file
void hloop(double h , double time, nbody sim){
    ofstream outfile("Adams-BashforthE(h).txt");
    outfile << setprecision(15);
    ofstream outfile2("Adams-MoultonE(h).txt");
    outfile2 << setprecision(15);
    ofstream outfile3("Foresth-RuthE(h).txt");
    outfile3 << setprecision(15);
    ofstream outfile4("RK4NE(h).txt");
    outfile4 << setprecision(15);
    ofstream outfile5("VerletE(h).txt");
    outfile5 << setprecision(15);

    for (int i=0; i<4; i++){
        
        double u = pow(10, -i) * h;
        int steps = time / u;

        RK4integrator(u, time, sim, false, 1., 1.);
        verlet(h, time, sim);
        FR(u, time, sim);
        AM(u, time, sim);
        AB(u, time, sim, 4);

        string En;
        ifstream E("Adams-BashforthE.txt");
    
        int l = 0;
        double E0;
        double Ef;
        while (getline(E, En)){
            if (l==0){
                E0 = stod(En);
            };

            if (l==steps-2){
                Ef = stod(En);
            };
            ++l;
            
        }  
        outfile << abs((Ef-E0)/E0) << ' ' << u << '\n';
        E.close();

        ifstream E2("Adams-MoultonE.txt");
        //string En2;
        int l2 = 0;
        double E02;
        double Ef2;
        while (getline(E2, En)){
            if (l2==0){
                E02 = stod(En);
            };

            if (l2==steps-2){
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

            if (l4==steps-2){
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

            if (l5==steps-2){
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

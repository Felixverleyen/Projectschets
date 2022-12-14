We have 5 files:

INT_Nbody.hpp: 
    This file contains the five integrator functions:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        RK4(h, time, sim): Runge-Kutta 4 integrator

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        Verlet(h, time, sim): Verlet integrator

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        FR(h, time, sim): Forest-Ruth integrator

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        AB(h, time, sim): Adams-Bashforth integrator
            Here we have implemented that one can choose the order of integration from 1 to 4

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        AM(h, time, sim): Adams-Moulton integrator
            The first 3 steps are done with the Adams-Bashforth orders 1 to 3.
            This method is implicit, we first make a prediction of the position with the Adams-Bashford method
            We then use the Adams-Moulton method to correct this prediction 
   
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    hloop: this function loops over different


CLASS_Nbody.hpp:
    This headerfile contains some assistance classes:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        Vec: a vector datatype to represent vectors in 3D-space

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        nbody: datatype that holds the masses, positions and velocities of all bodies at a specifiek time

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        counter: dataype that counts the number of driverfunction evaltuations for the integrater that is used


FUNC_Nbody.hpp:
    This headerfile contains some assistance functions:
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        a(i, sim): funtion to calculate the accelartion of one body due to the N-1 other bodies

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Energy(sim): function to calculate the energie of the system at a given time

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        time_step_scale(sim, i, schaal): adaptive time step

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        init_sim(file): reads initial conditions from the file and initiates an nbody-type  system


MAIN_Nbody.cpp:
    Here all the previous headerfiles are brought together.
    In the main-function, the initial conditions from the file "inintial_cond.txt" are read into an nbody datatype sim.
    Further one has to specify: 
        time: the Time over which one wants to integrate (a multiple of 2.9e9)
        h: the step size
        RK4/Verlet/FR/AB/AM: the integrator one wants to use

initial_cond.txt:
    File that contains initial conditions of N-body system 


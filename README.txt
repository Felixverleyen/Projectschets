We have 6 files: int_nbody.hpp, class_nbody.hpp, func_nbody.hpp, loops_nbody.hpp, MAIN_Nbody.cpp and initial_cond.txt


int_nbody.hpp: 
    This file contains the five integrator functions:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        RK4integrator(h, time, sim, adapt, scale, power): Runge-Kutta 4 integrator

            if adapt=true this intagrator will use an adaptive time step scheme with scale and power 2 parameters of the scheme.
            The explanation and source code of the scheme are given at func_nbody.hpp line 89

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        Verlet(h, time, sim): Verlet integrator

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        FR(h, time, sim): Forest-Ruth integrator

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        AB(h, time, sim): Adams-Bashforth integrator

            Here we have implemented that one can choose the order of integration from 1 to 4
            The first positions of the higher orders are determined with the lower order Adams-Bashforth method

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        AM(h, time, sim): Adams-Moulton integrator

            We determine the first 4 positions with the Runge-Kutta method
            This method is implicit, we first make a prediction of the position with the Adams-Bashford method
            We then use the Adams-Moulton method to correct this prediction 
   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////



class_nbody.hpp:
    This headerfile contains some assistance classes:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        Vec: a vector datatype to represent vectors in 3D-space

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        nbody: datatype that holds the masses, positions, velocities and (if needed) accelarations of all bodies
                at a specific time. Only the AB and AM methods need to hold the accelarations.

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        counter: dataype that can be used to count the number of driverfunction evaluations for the integrater 
                 that is used.

            Each time that the accelaration is calculated for a body (calling of function a(), see FUNC_Nbody) 
            the counter will go +1. The driver function is an N-dimensional array (N bodies) of 3-Vectors a(i)
            the counter goes up for every a(i) so we should divide by N to get the correct number of driver-
            function evaluations (it should be independent of #bodies=N). It should also be divided by the amount
            of time steps to give the average evaluations per unit of simulated time. The counter is reset after 
            integration with any of the intagrators.

    

loops_nbody.hpp:
    This headerfile is not needed for the operation of our simulations. We have included it here for completeness.
    It generates the proper output for when we want to see how the relative energy error evolves with varying 
    time step (for fixed h) and generating output for varying parameters of the adaptive scheme.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        hloop(h, time, sim): for varying the time step
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        adapt_loop(h, time, scale, sim): for varying parameters of the adaptive scheme



func_nbody.hpp:
    This headerfile contains some assistance functions:
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        a(i, sim): function to calculate the accelaration of one body due to the N-1 other bodies

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Energy(sim): function to calculate the energy of the system at a given time

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        time_step_scale(sim, i, schaal): adaptive time step (see source code for explanation of our method)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        init_sim(file): reads initial conditions from the file and initiates an nbody-type  system




MAIN_Nbody.cpp:
    Here all the previous headerfiles are brought together.
    In the main-function, the initial conditions from the file "inintial_cond.txt" are read into an nbody datatype sim.
    Further one has to specify: 
        time: the Time over which one wants to integrate (in units of 2.9722 s \approx 94 years)
        h: the step size
        RK4/Verlet/FR/AB/AM: the integrator one wants to use
    
    Here we can also initialize the loops that provide output for varying h and varying parameters of the adaptive time sheme

initial_cond.txt:
    File that contains initial conditions of N-body system. Also the used simulation units are explained


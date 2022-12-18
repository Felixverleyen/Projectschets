We have 6 files: INT_Nbody.hpp, CLASS_Nbody.hpp, FUNC_Nbody.hpp, LOOPS_Nbody.hpp, MAIN_Nbody.cpp and initial_cond.txt

INT_Nbody.hpp: 
    This file contains the five integrator functions:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        RK4(h, time, sim, adapt, scale, power): Runge-Kutta 4 integrator

            if adapt=true this intagrator will use an adaptive time step scheme with scale and power 2 parameters of the scheme.
            The explanation and source code of the scheme are given at func_nbody.hpp line 89

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
   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////


CLASS_Nbody.hpp:
    This headerfile contains some assistance classes:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        Vec: a vector datatype to represent vectors in 3D-space

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        nbody: datatype that holds the masses, positions, velocities and (if needed) accelarations of all bodies
               at a specific time.

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        counter: dataype that can be used to count the number of driverfunction evaltuations for the integrater 
                 that is used.

            Each time that the accelaration is calculated for a body (calling of function a(), see FUNC_Nbody) 
            the counter will go +1. The driver function is an N-dimensional array (N bodies) of 3-Vectors a(i)
            the counter goes up for every a(i) so we should divide by N to get the correct number of driver-
            function evaluations (independent of bodies). It should also be divided by the amount of time steps
            to give the average evaluations per unit of simulated time. The counter is reset after integration 
            with any of the intagrators.

    


LOOPS_Nbody.hpp:
    This headerfile is not needed for the operation of our simulations. We have included it here for completeness
    It generates the proper output for when we want to see how the relative energy error evolves with varying 
    time step (for fixed h) and generating output for varying parameters of the adaptive scheme.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        hloop(h, time, sim): for varying the time step
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        adapt_looph, time, scale, sim): for varying parameters of the adaptive scheme




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
        time: the Time over which one wants to integrate (in units of 2.9722 s \approx 94 years)
        h: the step size
        RK4/Verlet/FR/AB/AM: the integrator one wants to use
    
    Here we can also initialize the loops that provide output for varying h and varying parameters of the adaptivve time sheme

initial_cond.txt:
    File that contains initial conditions of N-body system. Also the used simulation units are explained


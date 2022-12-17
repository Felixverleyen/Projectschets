//Forest-Ruth
double theta= 1.351207;

nbody FR_step(double h, nbody sim){
    nbody sim2 = sim;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_r(i, sim.r(i)+h*theta*sim.v(i)*0.5);
    }
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_v(i, sim.v(i)+theta*h*a(i,sim)*(-1));}
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_r(i, sim.r(i)+h*(1-theta)*sim.v(i)*0.5);}
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_v(i, sim.v(i)+h*(1-2*theta)*a(i,sim)*(-1));}
    sim = sim2;
    
    for(int i = 0; i < sim.bodies(); i++){
        sim2.swap_r(i, sim.r(i)+ h*(1-theta)*sim.v(i)*0.5);}
    sim = sim2;

    for(int i = 0; i < sim.bodies(); i++){
        sim.swap_v(i, sim.v(i)+theta*h*a(i,sim)*(-1));}
    sim = sim2;
    
    for(int i = 0; i < sim.bodies(); i++){
        sim.swap_r(i, sim.r(i)+h*theta*sim.v(i)*0.5);}
    sim = sim2;
    
    return sim;
}
int main(){
    string initial_i;
    nbody sim;
    int l = 0;
    fstream initialNfile("Initial_cond.txt");

    while (getline(initialNfile, initial_i)){
        
           

            double m = stod(initial_i.substr(2,4));

            double rx = stod(initial_i.substr(7, 4));
            double ry = stod(initial_i.substr(12, 4));
            double rz = stod(initial_i.substr(17, 4));
            double vx = stod(initial_i.substr(22, 4));
            double vy = stod(initial_i.substr(27, 4));
            double vz = stod(initial_i.substr(32, 4));

            
            Vec pos{rx, ry, rz};
            Vec vel{vx,vy,vz};

            sim.add_mass(m);
            sim.add_pos(pos);
            sim.add_vel(vel);
        
        };
        
    
    initialNfile.close();
    };
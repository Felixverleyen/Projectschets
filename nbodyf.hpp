#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
using namespace std;

class Vec{
    double _x;
    double _y;
    double _z;

    public:
    Vec(){_x=0; _y=0;_z=0;}
    Vec(double x, double y, double z){_x=x; _y=y; _z=z;}

    double x() const{return _x;}
    double y() const{return _y;}
    double z() const{return _z;}

    double norm() const{return pow(_x*_x + _y*_y+ _z*_z, 0.5);}
    double norm2() const{return _x*_x + _y*_y+ _z*_z;}
    double norm3() const{double r=pow(_x*_x + _y*_y+_z*_z, 0.5); return r*r*r;}

    Vec& operator+=(Vec v) {
        _x += v._x;
        _y += v._y;
        _z += v._z;
        return *this;}

    Vec& operator-=(Vec v) {
        _x -= v._x;
        _y -= v._y;
        _z -= v._z;
        return *this;}

    Vec& operator*=(double s) {
        _x *= s;
        _y *= s;
        _z *= s;
        return *this;}

    Vec& operator/=(double s) {
        _x /= s;
        _y /= s;
        _z /= s;
        return *this;}


    };
Vec operator*(Vec a, double s) { return a *= s; }
Vec operator*(double s, Vec b) { return b *= s; }
Vec operator/(Vec a, double s) { return a /= s; }

Vec operator+(Vec a, Vec b) {return a+=b;}
Vec operator-(Vec a, Vec b) { return a -= b; }

class nbody{

    int N;
    vector<double> _mass;
    vector<Vec> _pos;
    vector<Vec> _vel;
    

    public:

    
    void set_N(int n) {N=n;};
    int bodies() const{return N;};

    //deze methods geven de waarde van een element in de lijst
    double m(int l) const{return _mass[l];};
    Vec r(int w) const{return _pos[w];};
    Vec v(int q) const{return _vel[q];};

    // methods om de lijsten op te vullen:
    void add_mass(double massa){
        _mass.push_back(massa);
    };

    void add_pos(Vec r){
        _pos.push_back(r);
    };

    void add_vel(Vec v){
        _vel.push_back(v);
    };

    vector<Vec> pos() const{return _pos;};
    //methods om een element in een lijst te vervangen, input: (positie in de lijst, nieuwe waarde)
    void swap_m(int k, double m){
        _mass[k] = m;
    }

    void swap_r(int h, Vec r){
        _pos[h] = r;
    }

    void swap_v(int f, Vec v){
        _vel[f] = v;
    }

};

class counter{
    int count = 0;

    public:
    void countplus() {count++;};
    int getcount() {return count;};
    void reset() {count = 0;};

    double order(int N, int steps) { double averaged_count = (1.* count) / (1. * steps * N);
    return averaged_count;};
};

double G=1;
counter drivercount;


// s = r / scale(normally =1AU) => h'=f(s)*h with f(s)= s^power , we vary the power to give different time_schemes
double time_step_scale(nbody sim,  double scale, double power){
    vector<double> diff;

    for (int i = 0; i < sim.bodies(); i++){
        for (int j = 0; j < sim.bodies(); j++){

            if (j != i){
                Vec d = sim.r(i) - sim.r(j);
                diff.push_back(d.norm());
            }

            else { continue;}
        
        }
    }
    
    double min = diff[0];
    
    for (int i = 0; i < diff.size(); i++){

        if(diff[i] < min){
            min = diff[i];
        }
    }
    
    return pow(min / scale, power);
}

void print(Vec a)
{ cout << a.x() << ' ' << a.y() << ' '<< a.z() << endl; }

Vec a(int i, nbody sim) {
    drivercount.countplus();
    Vec ri= sim.r(i);
    Vec ar={0,0,0};

    for (int j =0; j<sim.bodies(); ++j){
        Vec rj = sim.r(j);
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

    for (int i=0; i<sim.bodies();++i){
        Vec ri = sim.r(i);
        double mi=sim.m(i);
        Vec vi= sim.v(i);
        E+=mi*vi.norm2()/2;

        for  (int j =0; j<sim.bodies(); ++j){
            Vec rj = sim.r(j);
            double mj= sim.m(j);
            Vec afst= ri-rj;

            if (i!=j){
                E-= 1/2*G*mi*mj/afst.norm();
                }
            }
    }

     return E;}


nbody init_sim(string file){
    string initial_i;
    nbody sim;
    int N = 0;
    int l = 0;
    fstream initialNfile(file);

    while (getline(initialNfile, initial_i)){
        if(l>4){
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

            ++N;
        }

        ++l;
        };

    sim.set_N(N);
    initialNfile.close();
    return sim;
}
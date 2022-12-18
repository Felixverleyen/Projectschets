#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
using namespace std;

//3D-Vector class:
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
    double norm3() const{return pow(_x*_x + _y*_y+_z*_z, 1.5);}

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

//A class that holds all information of the simulation at time t_n (N masses, N positions, N velocities and if needed N accelarations)
class nbody{

    int N;
    vector<double> _mass;
    vector<Vec> _pos;
    vector<Vec> _vel;
    vector<Vec> _a;

    public:

    
    void set_N(int n) {N=n;};
    int bodies() const{return N;};

    //these methods give the values corresponding to specified bodies
    double m(int l) const{return _mass[l];};
    Vec r(int w) const{return _pos[w];};
    Vec v(int q) const{return _vel[q];};
    Vec a(int z) const{return _a[z];}

    //Methods to fill the lists when initializing the simulation (used in init_sim() function, see func_nbody)
    void add_mass(double massa){
        _mass.push_back(massa);
    };

    void add_pos(Vec r){
        _pos.push_back(r);
    };

    void add_vel(Vec v){
        _vel.push_back(v);
    };

    void add_a(Vec a){
        _a.push_back(a);
    }

    vector<Vec> pos() const{return _pos;};

    //methods to change a specific value in sim, input: (index of body, new value)
    void swap_m(int k, double m){
        _mass[k] = m;
    }

    void swap_r(int h, Vec r){
        _pos[h] = r;
    }

    void swap_v(int f, Vec v){
        _vel[f] = v;
    }

    void swap_a(int j, Vec a){
        _a[j] = a;
    }

};

//class used to count the driverfunction evaluations
class counter{
    int count = 0;

    public:
    void countplus() {count++;};
    int getcount() {return count;};
    void reset() {count = 0;};

    
    double order(int N, int steps) {            //gives the number of driver function evaluations of an integrator
    
    double averaged_count = (1.* count) / (1. * steps * N); 
    return averaged_count;};

}; 

double G=1; 
counter drivercount;
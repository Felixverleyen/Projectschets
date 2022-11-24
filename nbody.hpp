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

    int bodies() const{return N;};

    //deze methods geven de waarde van een element in de lijst
    double m(int l) const{return _mass[l];};
    Vec p(int w) const{return _pos[w];};
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

    //methods om een element in een lijst te vervangen, input: (positie in de lijst, nieuwe waarde)
    void swap_m(int k, double m){
        _mass[k] = m;
    }

    void swap_p(int h, Vec r){
        _pos[h] = r;
    }

    void swap_v(int f, Vec v){
        _vel[f] = v;
    }
};
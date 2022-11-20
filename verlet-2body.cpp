#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
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

    double norm() const{return sqrt(_x*_x + _y*_y+ _z*_z);}
    double norm2() const{return _x*_x + _y*_y+ _z*_z;}
    double norm3() const{double r=sqrt(_x*_x + _y*_y+_z*_z); return r*r*r;}

    void print(const vector<double>& v) {
    cout << v.size() << ": ";
    for (double elem : v) cout << elem << ' ';
    cout << endl;
    };
    void print(Vec a)
    { cout << a.x() << ' ' << a.y() << endl; }

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

Vec a(double mu, Vec r) {return -mu/r.norm3() * r; }

double energy(double mu, Vec r, Vec v) {return 0.5 * v.norm2() - mu/r.norm();}


int main() {
    
    // file
    ofstream outfile("verlet2b.txt");
    outfile << setprecision(15);

    //velocity verlet

    // initial conditions
    Vec r0(1,0,0.5);
    Vec v0(0.1,0.06,0.02);
    
    const double t0 = 0.;
    const double t_end = 1000.;
    const double mu = 0.01;
    const double h = 0.01; 

    // acceleration
    Vec a0 = -mu/r0.norm3() * r0;
        
    // energie
    double E0 = 0.5 * v0.norm2() - mu/r0.norm();
    
    outfile << t0 << ' ' << r0.x() << ' ' << r0.y() << ' ' << r0.z() << ' ' << v0.x() << ' ' << v0.y() << ' ' << v0.z() << ' ' << E0 << '\n';

    //first time step (n=0)
    Vec vhalf = v0 + 0.5*h*a0;
    Vec r = r0 + h*vhalf;
    Vec vfull = vhalf + 0.5*h*a(mu,r);

    double E = energy(mu,r,vfull);

    outfile << t0+h << ' ' << r.x() << ' ' << r.y() << ' ' << r.z() << ' ' << vfull.x() << ' ' << vfull.y() << ' ' << vfull.z() << ' ' << E << '\n';

    //loop (n>=1)
    for (double t=t0+2*h; t<t_end+h/2.; t+=h) {
        vhalf = vfull + 0.5*h*a(mu,r);
        r = r + h*vhalf;
        vfull = vhalf + 0.5*h*a(mu,r);
        E = energy(mu,r,vfull);
        outfile << t << ' ' << r.x() << ' ' << r.y() << ' ' << r.z() << ' ' << vfull.x() << ' ' << vfull.y() << ' ' << vfull.z() << ' ' << E << '\n';

    }

    outfile.close();

};
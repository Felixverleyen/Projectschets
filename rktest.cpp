//runge-kutta-4 integrator
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

double norm() const{return pow(_x*_x + _y*_y+ _z*_z, 0.5);}
double norm2() const{return _x*_x + _y*_y+ _z*_z;}
double norm3() const{double r=pow(_x*_x + _y*_y+_z*_z, 0.5); return r*r*r;}

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

Vec a(Vec r, double m){
    return -1 * m * (r / r.norm3());
}

array<Vec,2> RK4_step(double h, Vec ri, Vec vi, double mu){
    Vec kv1 = h *a(ri, mu);
    Vec kx1 = h * vi;

    Vec lv = ri + 0.5 * kv1;
    Vec lx = vi + 0.5 * kx1;

    Vec kv2 = h *a(lv, mu);
    Vec kx2 = h * lx;

    lv = ri + 0.5 * kv2;
    lx = vi + 0.5 * kx2;

    Vec kv3 = h *a(lv, mu);
    Vec kx3 = h * lx;

    lv = ri + kv3;
    lx = vi + kx3;

    Vec kv4 = h *a(lv, mu);
    Vec kx4 = h * lx;

    Vec vo = vi + (1./6.) * (kv1 + 2 * kv2 + 2 * kv3 + kv4);
    Vec ro = ri + (1./6.) * (kx1 + 2 * kx2 + 2 * kx3 + kx4);
    array<Vec,2> opl = {ro, vo};

    return opl;

}

void RK4(double h, double mu, double time, Vec r0, Vec v0){
        int steps = int(time / h);
        
        Vec r = RK4_step(h, r0, v0, mu)[0];
        Vec v = RK4_step(h, r0, v0, mu)[1];
        ofstream outfile("project.txt");
        outfile << setprecision(15);
        for (int i = 0; i < steps; i++){
            r = RK4_step(h, r, v, mu)[0];
            v = RK4_step(h, r, v, mu)[1];
            
            outfile << r.x() << ' ' << r.y() << ' ' << r.z() << '\n'; //<< ' ' << v.x() << ' ' << v.y() << ' ' << v.z() 
        }
}


int main(){
    ofstream outfile("project.txt");
    outfile << setprecision(15);
    
     
    Vec r0(1.2, 0, 0);
    Vec v0(0., 0.5, 0.);
    
    double mu = 1;
    double h = 0.001;
    double t = 50;
    RK4(h, mu, t, r0, v0);
    Vec p = a(r0, mu);
    
}


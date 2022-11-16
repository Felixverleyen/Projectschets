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

// telkens voor andere ri ar berekenen door te sommeren over elke andere rj
Vec a(int i){
    Vec ri= pos[0,i];
    Vec rj= pos[0,j];
    Vec ar;
for (int j =0; j<N, ++j){
    Vec rj = pos[j];
    double mj= mass[j];
    double mi=mass[i];
    Vec afst= ri-rj;
 if (i!=j){
    ar-= G*mj*afst/afst.norm3();
    return ar;

 }}} ;

Vec a(Vec r, double m){
    return -1 * m * (r / r.norm3())
}
//runge-kutta-4 integrator

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

    Vec v = vi + (1/6) * (kv1 + 2 * kv2 + 2 * kv3 + kv4);
    Vec r = ri + (1/6) * (kx1 + 2 * kx2 + 2 * kx3 + kx4);
    array<Vec,2> opl = {r, v};

    return opl

}

void RK4(double h, double time, Vec r0, Vec v0){
        int steps = int(time / h);
        for (int i = 0; i < steps; i++){
            Vec r = RK4_step(h, r)
        }
}


int main(){
    array<double> pos;
    array<double> mass;
    vector<double> list;
    ofstream outfile("project.txt");
    outfile << setprecision(8);
    
//begincondities (pos en speed)

    //outfile 

// initialiseren voor elke i (eerste stappen tot algoritme kan beginnen)

     

for (double t = 0.00; t <= 0.1; t+=h){
    list.pushback(pos);
    for(int i=0; i<N; i++){
        //programma
    Vec rn = pos[t,i];
    Vec vn = speed[t,i];


    



        
        //pos[t,i]= waarde van ri op tijd t;
        // speed[t,i]= ...


        // voor elke i telkens naar andere file schrijven

    };

    // na volledige tijdstap energie berekenen en uitschrijven

}
   
    outfile.close();
};

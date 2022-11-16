#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <string>
using namespace std;


double G=8;
double mu=0.01;
double h=0.01;
int const N=4;
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

    array<Vec,N> npos = {};
    array<Vec,N> opos = {Vec(2,2,1),Vec(2,1,-1),Vec(2,1,3),Vec(4,2,1)};
    array<double,N> mass = {4
    ,4,5,6};
    array<Vec,N> ospeed = {Vec(3,1,2),Vec(3,1,2),Vec(0,1,2),Vec(0,1,2)}; 
    array<Vec,N> nspeed = {};
    Vec list[N][100];
    vector<double> speedlist;
   
void print(Vec a)
{ cout << a.x() << ' ' << a.y() << ' '<< a.z() << endl; }

Vec a(int count ,int i) {
    // aanpassen naar verschillende tijdstippen met count
    Vec ri= opos[i];
    Vec ar={0,0,0};
for (int j =0; j<N; ++j){
    Vec rj = opos[j];
    double mj= mass[j];
    double mi=mass[i];
    Vec afst= ri-rj;
 if (i!=j){
    ar-= G*mj*afst/afst.norm3();
    
 } 
}
return ar;}
double E(){
    double E=0;
    for (int i=0; i<N;++i){
        Vec ri = opos[i];
        double mi=mass[i];
        Vec vi= ospeed[i];
        E+=mi*vi.norm2()/2;
        for  (int j =0; j<N; ++j){
            Vec rj = opos[j];
            double mj= mass[j];
            Vec afst= ri-rj;
            if (i!=j){
             E-= 1/2*G*mi*mj/afst.norm();
                }
    }} return E;}

 
int main(){
     ofstream outfile("projectpos.txt");
    outfile << setprecision(8);
    ofstream outfile2("projectspeed.txt");
    outfile2 << setprecision(8);
    ofstream outfile3("projectenergy.txt");
    outfile2 << setprecision(8);
    int count=0;
    ofstream time("projecttime.txt");
    time << setprecision(8);
    

//begincondities (pos en speed)

    //outfile 

// initialiseren voor elke i (eerste stappen tot algoritme kan beginnen)

     
    double theta= 1.351207;
    for (double t = 0.00; t <= 20; t+=h){
    
        
    time << t << '\n';
    outfile3 << E() << '\n';
    for(int i=0; i<N; i++){ 
    Vec rn = opos[i];
    Vec vn = ospeed[i];

    

    rn += h*theta*vn/2;
    vn += theta*h*a(count,i);
    rn += h*(1-theta)*vn/2;
    vn += theta*h*(1-2*theta)*a(count,i);
    rn += h*(1-theta)*vn/2;
    vn += theta*h*a(count,i);
    rn += h*theta*vn/2;

    npos[i]= rn;
    nspeed[i]= vn;
    
   outfile << rn.x() << ' ' <<rn.y()<<' '<< rn.z()<<' ' ;
    outfile2 << vn.x() << ' ' <<vn.y()<<' '<< vn.z()<<' ';
        
        //pos[t,i]= waarde van ri op tijd t;
        // speed[t,i]= ...


        // voor elke i telkens naar andere file schrijven

    };
    outfile  << '\n';
    outfile2 << '\n';
    int count += 1;
    ospeed = nspeed;
    opos = npos; 
    

}
  
};

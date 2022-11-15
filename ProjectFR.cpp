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
int const N=2;
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
    array<Vec,N> opos = {Vec(1,2,3),Vec(2,4,1)};
    array<double,N> mass = {4,10};
    array<Vec,N> ospeed = {Vec(1,5,3),Vec(-2,4,0)}; 
    array<Vec,N> nspeed = {};
    Vec list[N][100];
    vector<double> speedlist;
   
    


Vec a(int s,int i) {
    Vec ri= opos[i];
    Vec ar;
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



 
int main(){
     ofstream outfile("projectpos.txt");
    outfile << setprecision(8);
    ofstream outfile2("projectspeed.txt");
    outfile2 << setprecision(8);
    outfile << opos[0].x() << ' ' <<opos[0].y();
    int count=0;

//begincondities (pos en speed)

    //outfile 

// initialiseren voor elke i (eerste stappen tot algoritme kan beginnen)

     
    double theta= 1.351207;
    for (double t = 0.00; t <= 0.1; t+=h){
    
        

    
    for(int i=0; i<N; i++){ 
    Vec rn = opos[i];
    Vec vn = ospeed[i];
    

    rn += 1/2*h*theta*vn;
     outfile << rn.x() << ' ' <<rn.y()<<' '<< rn.z() ;
    vn += theta*h*a(count,i);
    outfile2 << vn.x() << ' ' <<vn.y()<<' '<< vn.z();
    rn += 1/2*h*(1-theta)*vn;
     outfile << rn.x() << ' ' <<rn.y()<<' '<< rn.z() ;
    vn += theta*h*(1-2*theta)*a(count,i);
    outfile2 << vn.x() << ' ' <<vn.y()<<' '<< vn.z();
    rn += 1/2*h*(1-theta)*vn;
     outfile << rn.x() << ' ' <<rn.y()<<' '<< rn.z() ;
    vn += theta*h*a(count,i);
    rn += 1/2*h*theta*vn;

    npos[i]= rn;
    nspeed[i]= vn;
    
   outfile << rn.x() << ' ' <<rn.y()<<' '<< rn.z()<<'\t' ;
    outfile2 << vn.x() << ' ' <<vn.y()<<' '<< vn.z()<<'\t';
outfile <<'\t';
outfile2 <<'\t';
        
        //pos[t,i]= waarde van ri op tijd t;
        // speed[t,i]= ...


        // voor elke i telkens naar andere file schrijven

    };
    outfile  << '\n';
    outfile2 << '\n';
    int count = count +1;
    ospeed = nspeed;
    opos = npos; 
    // na volledige tijdstap energie berekenen en uitschrijven

}
  
};
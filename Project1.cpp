(#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
using namespace std;
void print(const vector<double>& v) {
cout << v.size() << ": ";
for (double elem : v) cout << elem << ' ';
cout << endl;
};
double G=8;
double mu=0.01;
double h=0.01;
double N=8;

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
Vec a(int i) {
    Vec ri= pos[i]
    Vec ri= pos[j]

    Vec ar;
for (int j =0; j < N, ++j){
    Vec rj = pos[j];
    double mj= mass[j];
    double mi=mass[i];
    Vec afst= ri-rj;
 if (i!=j){
    ar-= G*mj*afst/afst.norm3();
    return ar;

 } ;
    



void print(Vec a)
{ cout << a.x() << ' ' << a.y() << endl; }



int main(){
    array<double> pos;
    array<double> mass;
    ofstream outfile("project.txt");
    outfile << setprecision(8);
    
    //begincondities (pos en speed)
    
    outfile << rn.x() << ' ' <<rn.y()<<' '<< vn.x() <<' '<< vn.y() << ' '
<< a(rn).x()<<' ' << a(rn).y() << '\n';

// initialiseren vooor elke i 

     

    outfile << rn.x() << ' ' <<rn.y()<<' '<< vn.x() <<' '<< vn.y() << ' '
<< a(rn).x()<<' ' << a(rn).y() << '\n';





for (int t = 0; t <= 0.1; t+=h){
    for(int i=0; i<N; i++){
        //programma


        //telkens naar andere file schrijven

    }

    // na volledige tijdstap energie berekenen en uitschrijven

}
   
    outfile.close();)
};
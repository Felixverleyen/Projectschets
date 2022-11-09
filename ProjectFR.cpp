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

// telkens voor andere ri ar berekenen door te sommeren over elke andere rj
Vec a(int i) {
    Vec ri= pos[0,i];
    Vec ar;
for (int j =0; j<N, ++j){
    Vec rj = pos[0,j];
    double mj= mass[j];
    double mi=mass[i];
    Vec afst= ri-rj;
 if (i!=j){
    ar-= G*mj*afst/afst.norm3();
    return ar;

 } ;





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

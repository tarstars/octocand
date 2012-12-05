//---------------------------------------------------------------------------


#pragma hdrstop

#include "matrix3.h"

//---------------------------------------------------------------------------

#include <iostream>
#include <sstream>

using namespace std;


ostream&
operator<<(ostream& os, const Matrix3 &r){
	for(int p = 0; p < 3; ++p){
		for(int q = 0; q < 3; ++q)
			os << r(p,q) << " ";
		os << endl;
	}
	return os;
}

istream&
operator>>(istream& is, Matrix3 &r){
 for (int p = 0; p < 3; ++p) {
	 for (int q = 0; q < 3; ++q)
		is >> r(p,q);
 }
	return is;
}

Matrix3
Matrix3::operator*(const Matrix3&r)const{

Matrix3 ret;
ret.dat[0][0]=dat[0][0]*r.dat[0][0]+dat[0][1]*r.dat[1][0]+dat[0][2]*r.dat[2][0];
ret.dat[1][0]=dat[1][0]*r.dat[0][0]+dat[1][1]*r.dat[1][0]+dat[1][2]*r.dat[2][0];
ret.dat[2][0]=dat[2][0]*r.dat[0][0]+dat[2][1]*r.dat[1][0]+dat[2][2]*r.dat[2][0];

ret.dat[0][1]=dat[0][0]*r.dat[0][1]+dat[0][1]*r.dat[1][1]+dat[0][2]*r.dat[2][1];
ret.dat[1][1]=dat[1][0]*r.dat[0][1]+dat[1][1]*r.dat[1][1]+dat[1][2]*r.dat[2][1];
ret.dat[2][1]=dat[2][0]*r.dat[0][1]+dat[2][1]*r.dat[1][1]+dat[2][2]*r.dat[2][1];

ret.dat[0][2]=dat[0][0]*r.dat[0][2]+dat[0][1]*r.dat[1][2]+dat[0][2]*r.dat[2][2];
ret.dat[1][2]=dat[1][0]*r.dat[0][2]+dat[1][1]*r.dat[1][2]+dat[1][2]*r.dat[2][2];
ret.dat[2][2]=dat[2][0]*r.dat[0][2]+dat[2][1]*r.dat[1][2]+dat[2][2]*r.dat[2][2];


return ret;
}

void test_mat(){
	stringstream sour_a ("11 12 13 21 22 23 31 32 33");

	Matrix3 a;
	sour_a >> a;

	stringstream sour_b ("7 13 22 63 12 8 1 3 2");

	Matrix3 b;
	sour_b >> b;

	cout << a << endl;
	cout << b << endl;

	cout << a*b << endl;
	}

	Poly
	Matrix3::characteristic()const{
	Poly ret(-1,
			dat[0][0]+dat[1][1]+dat[2][2],
			dat[0][1]*dat[1][0]-dat[0][0]*dat[1][1]-dat[0][0]*dat[2][2]+dat[0][2]*dat[2][0]-dat[1][1]*dat[2][2]+dat[1][2]*dat[2][1],
			dat[0][0]*dat[1][1]*dat[2][2]-dat[0][0]*dat[1][2]*dat[2][1]-dat[0][1]*dat[1][0]*dat[2][2]+dat[0][1]*dat[1][2]*dat[2][0]+dat[0][2]*dat[1][0]*dat[2][1]-dat[0][2]*dat[1][1]*dat[2][0]);
	return ret;
	}

#pragma package(smart_init)

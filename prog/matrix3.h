//---------------------------------------------------------------------------

#ifndef matrix3H
#define matrix3H
//---------------------------------------------------------------------------

#include <iostream>

#include "poly.h"

class Matrix3 {
	double dat[3][3];
public:
	double& operator()(int p, int q) {return dat[p][q];}
	const double& operator()(int p, int q)const {return dat[p][q];}

	Matrix3 operator*(const Matrix3& )const;

	Poly characteristic()const;

	std::ostream& operator>>(std::ostream& os)const;
	std::istream& operator<<(std::istream& is);
};

std::ostream& operator<<(std::ostream& os, const Matrix3& r);
std::istream& operator>>(std::istream& is, Matrix3& r);


#endif

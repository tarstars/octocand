//---------------------------------------------------------------------------

#ifndef Tensor3H
#define Tensor3H

#include <iostream>

#include "matrix3.h"
#include "vector3.h"

//---------------------------------------------------------------------------

class Tensor3 {
	double dat[3][3][3];
public:
	double& operator()(int, int, int);
	const double& operator()(int, int, int) const;
	Matrix3 crist(const Vector3&, const Matrix3& epsilon, const Tensor3& pc)const;

	Tensor3 operator*(const Matrix3&);

	std::ostream& operator >>(std::ostream& )const;
};

std::ostream& operator <<(std::ostream& ,const Tensor3&);



#endif

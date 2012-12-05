//---------------------------------------------------------------------------

#ifndef tensorH
#define tensorH
//---------------------------------------------------------------------------

#include <iostream>

#include "matrix3.h"
#include "vector3.h"

class Tensor3;

class Tensor {
	double dat[3][3][3][3];
public:
	double& operator()(int, int, int, int);
	const double& operator()(int, int, int, int) const;
	Matrix3 crist(const Vector3&)const;
    Matrix3 make_piezo_crist(const Vector3& n, const Matrix3& epsilon, const Tensor3& pc) const;

	Tensor operator*(const Matrix3&);

	std::ostream& operator >>(std::ostream& )const;
};

std::ostream& operator <<(std::ostream& ,const Tensor&);

#endif

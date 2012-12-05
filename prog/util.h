//---------------------------------------------------------------------------

#ifndef utilH
#define utilH

//---------------------------------------------------------------------------

#include <iostream>
#include <complex>
#include <vector>
#include "matrix3.h"
#include "vector3.h"
#include "tensor.h"
#include "tensor3.h"
#include "types.h"

Tensor
make_tetragonal_tensor(double,double,double,double,double,double);

Tensor
make_trigonal_tensor(double,double,double,double,double,double,double);

Tensor3
make_trigonal_piezo_tensor_3m(double,double,double,double);

Tensor3
make_trigonal_piezo_tensor_32(double,double);

Matrix3 rotx(double);
Matrix3 roty(double);
Matrix3 rotz(double);
Matrix3 euler(double, double, double);

Matrix3 eps2mat(double,double);

VD rho2v(const VCD& , double);

VD arr_re(const VCD&);

Vector3 pq2vec(int, int, int);

Vector3 get_polarization( const Matrix3&, double);

Vector3 ort(const Vector3& );

Vector3 slow_normal(const Vector3& n,
                    int ind,
                    const Matrix3& epsilon,
                    const Tensor& tens,
                    const Tensor3& piezotens);

std::ostream&
operator<<(std::ostream& os, const VCD& r);

std::ostream&
operator<<(std::ostream& os, const VD& r);

#endif

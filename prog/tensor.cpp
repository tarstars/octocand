//---------------------------------------------------------------------------


#pragma hdrstop

#include "tensor.h"
#include "tensor3.h"
#include <iostream>

using namespace std;

//---------------------------------------------------------------------------


Matrix3
Tensor::crist(const Vector3& v)const {
	Matrix3 ret;
	for (int i = 0; i < 3; i++) {
		for (int k = 0; k < 3; k++) {
			ret(i,k)=0;
			for (int j = 0; j < 3; j++) {
				for (int l = 0; l < 3; l++) {
					ret(i,k)+=dat[i][j][k][l]*v(j)*v(l);
				}
			}
		}
	}
	return ret;

}

Matrix3
Tensor::make_piezo_crist(const Vector3& n, const Matrix3& epsilon, const Tensor3& pc) const {
	Matrix3 ret;
	double eps=0;
	for (int p = 0; p < 3; p++) {
		for (int q = 0; q < 3; q++) {
			eps += epsilon(p,q)*n(p)*n(q);
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int k = 0; k < 3; k++) {
			ret(i,k)=0;
			double gamma_i=0, gamma_k=0;
			for (int j = 0; j < 3; j++) {
				for (int l = 0; l < 3; l++) {
					ret(i,k)+=dat[i][j][k][l]*n(j)*n(l);
					gamma_i+=pc(l,i,j)*n(j)*n(l);
					gamma_k+=pc(l,k,j)*n(j)*n(l);
				}
			}
			ret(i,k)+=gamma_i*gamma_k/eps;
		}
	}
return ret;
}

Tensor
Tensor::operator*(const Matrix3& mat){
Tensor ret;
for (int p=0; p < 3; p++) {
	for (int q = 0; q < 3; q++) {
		for (int r = 0; r < 3; r++) {
			for (int s = 0; s < 3; s++) {
				ret(p,q,r,s)=0;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						for (int k = 0; k < 3; k++) {
							for (int l = 0; l < 3; l++) {
								ret(p,q,r,s) += (*this)(i,j,k,l)*mat(p,i)*
																 mat(q,j)*
																 mat(r,k)*
																 mat(s,l);
							}
						}
					}
				}
			}
		}
	}
}

return ret;
}

double&
Tensor:: operator()(int i, int j, int k, int l) {
	return dat [i][j][k][l];
}

const double&
Tensor::operator()(int i, int j, int k, int l) const {
	return dat [i][j][k][l];
}

ostream&
Tensor::operator >>(ostream& os)const{
for (int q = 0; q < 3; q++) {
	for (int s = 0; s < 3; s++) {
		for (int p = 0; p < 3; p++) {
			for (int r = 0; r < 3; r++) {
				os << (*this)(p,q,r,s) << " ";
			}
			os << "\t\t";
		}
		os << endl;
	}
	os << endl;
}

return os;
}

ostream&
operator <<(ostream& os,const Tensor& r){
return r >> os;
}

#pragma package(smart_init)

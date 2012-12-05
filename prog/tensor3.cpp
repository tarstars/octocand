//---------------------------------------------------------------------------


#pragma hdrstop

#include "tensor3.h"
#include <iostream>

using namespace std;

//---------------------------------------------------------------------------


Tensor3
Tensor3::operator*(const Matrix3& mat){
Tensor3 ret;
for (int p=0; p < 3; p++) {
	for (int q = 0; q < 3; q++) {
		for (int r = 0; r < 3; r++) {
		ret(p,q,r)=0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						ret(p,q,r) += (*this)(i,j,k)*mat(p,i)*
													 mat(q,j)*
													 mat(r,k);
					}
				}
			}
		}
	}
}

return ret;
}

double&
Tensor3::operator()(int i, int j, int k){
	return dat [i][j][k];
}

const double&
Tensor3::operator()(int i, int j, int k) const {
	return dat [i][j][k];
}

ostream&
operator <<(ostream& os,const Tensor3& tens) {
for (int p = 0; p < 3; ++p) {
	for (int q = 0; q < 3; ++q) {
		for (int r = 0; r < 3; ++r) {
			 os << tens(p,q,r) << " ";
			}
			os << endl;
		}
		os << endl;
	}
	return os;
}

#pragma package(smart_init)

//---------------------------------------------------------------------------

#ifndef deriv_valsH
#define deriv_valsH

#include <iostream>
#include <complex>

//---------------------------------------------------------------------------

class DerivVals {
public:
	std::complex<double> d0;
	std::complex<double> d1;
	std::complex<double> d2;
	DerivVals (): d0(0),d1(0),d2(0){}
	DerivVals (std::complex<double> xd0, std::complex<double> xd1, std::complex<double> xd2): d0(xd0),d1(xd1),d2(xd2){}
	std::ostream& operator >>(std::ostream& os) const{
		return os<<d0<<" "<<d1<<" "<<d2;
	}
};

std::ostream& operator<<(std::ostream& os, const DerivVals& r);


#endif

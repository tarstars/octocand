//---------------------------------------------------------------------------

#ifndef polyH
#define polyH

#include <vector>
#include <complex>
#include <cmath>

#include "deriv_vals.h"
#include "types.h"

//---------------------------------------------------------------------------

class Poly {
	VCD dat;
public:
	Poly (){}
	Poly (double a1,double a0):dat(2) {
		dat[0]=a0;
		dat[1]=a1;
	}
	Poly (double a2,double a1,double a0):dat(3) {
		dat[0]=a0;
		dat[1]=a1;
		dat[2]=a2;
	}
	Poly (double a3,double a2,double a1,double a0):dat(4) {
		dat[0]=a0;
		dat[1]=a1;
		dat[2]=a2;
		dat[3]=a3;
	}

	Poly (std::istream&);

	std::ostream& operator >>(std::ostream& os) const{
		for (int i = (int)dat.size()-1; i >= 0 ; i--) {
			os<<dat[i]<<" ";
		}
		return os;
	}

	DerivVals vals_horn (CD q);

	VCD
	all_roots()const;


	/*DerivVals vals(double q) {
		DerivVals ret;
		for (int i = 0; i < (int)dat.size(); i++) {
			if (i>1) ret.d2 += i * (i - 1) * dat[i] * pow(q, i - 2);
			if (i>0) ret.d1 += i * dat[i] * pow(q, i - 1) ;

			ret.d0 += dat[i]*pow(q,i);
		}
		return ret;
	} */

	CD root();

	Poly bezu_div(const CD& r);

	int deg() const{return dat.size() -1;}
};

std::ostream& operator<<(std::ostream& os, const Poly& r);

#endif

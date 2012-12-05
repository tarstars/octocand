//---------------------------------------------------------------------------


//#pragma hdrstop

#include "poly.h"
#include "types.h"
#include "util.h"

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

//#include <vcl.h>
//#pragma hdrstop

#include <algorithm>
#include <vector>
#include <iostream>
//#include <tchar.h>
#include <string>

using namespace std;

//---------------------------------------------------------------------------

Poly::Poly (std::istream& is){
CD inp;
while (is >> inp)
	dat.push_back(inp);
reverse(dat.begin(), dat.end());
}

ostream& operator<<(ostream& os, const Poly& r){
	return r >> os;
}

	DerivVals Poly::vals_horn (CD q) {
		DerivVals ret;
		for (int i = (int)dat.size()-1; i >= 0; i--) {
			ret.d2=ret.d2*q+ret.d1;
			ret.d1=ret.d1*q+ret.d0;
			ret.d0=ret.d0*q+dat[i];
		}
		ret.d2 *= 2;
		return ret;
	}

	CD Poly::root() {
		CD ret(1e-19,0);
		int meter = 50;
		double eps=1e-15;
		int n=dat.size() - 1;
		while (meter > 0) {
			DerivVals pv=vals_horn(ret);
			if (std::abs(pv.d0)<eps) {
				return ret;
				}
			CD g=pv.d1/pv.d0;
			CD h=g*g-pv.d2/pv.d0;
			CD k=double(n-1)*(double(n)*h-g*g);
			k=sqrt(k);
			CD a;
			if (abs(g+k)>abs(g-k)) {
				a=double(n)/(g+k);
			} else a=double(n)/(g-k);

			if (abs(a)<eps) {
				return ret;
			}

			ret-=a;

			--meter;
		}
		return ret;
	}


void test_poly()
{
	Poly a(1, 2);
	Poly b(1, -8, 15);
	Poly c(1, 3, 3, 1);
	cout << a << endl;
	cout << b << endl;
	cout << c << endl;


	/*DerivVals test (1,2,3);
	cout << test << endl;*/

	DerivVals p=c.vals_horn(1);
	cout << p << endl;

	cout << c.root() << endl;

}

VCD
Poly::all_roots()const{
	VCD ret;
	Poly mute = *this;

	try{
	while (mute.deg()>0){
		CD root = mute.root();
		ret.push_back(root);
		mute = mute.bezu_div(root);
	}
    }catch(string msg){
	stringstream dest;
	dest << "Poly::all_roots() " << *this << endl << msg;
	throw dest.str();
	}


	return ret;
}

Poly
Poly::bezu_div(const CD& x0) {
	Poly ret;
	ret.dat.resize(dat.size()-1);
	for (int i = (int) dat.size()-2; i > -1; i--) {
		if (i==(int) dat.size()-2){
			ret.dat[i]=dat[i+1];
		}
		else {
			ret.dat[i]= ret.dat[i+1]*x0+dat[i+1];
		}
	}

	return ret;
}

//---------------------------------------------------------------------------

//#pragma package(smart_init)

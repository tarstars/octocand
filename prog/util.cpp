//---------------------------------------------------------------------------


//#pragma hdrstop

#include "util.h"

//---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>

using namespace std;

#include "matrix3.h"
#include "vector3.h"
#include "tensor.h"
#include "tensor3.h"
#include "util.h"


int jum(int p, int q) {
	int ind [3][3]={
		{0,5,4},
		{5,1,3},
		{4,3,2}
		};

	return ind [p][q];
}

void test_jum() {
	if (jum(0,0)==0) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(1,1)==1) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(2,2)==2) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(1,2)==3) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(0,2)==4) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(0,1)==5) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(2,1)==3) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(2,0)==4) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
	if (jum(1,0)==5) {cout << "OK" << endl;} else {cout << "!!!EPIC FAIL!!!" << endl;}
}

Tensor
make_tetragonal_tensor(	double c11,
						double c12,
						double c13,
						double c33,
						double c66,
						double c44) {
	Tensor ret;
	double el_c [6][6]={
		{c11, c12, c13,   0,   0,   0},
		{c12, c11, c13,   0,   0,   0},
		{c13, c13, c33,   0,   0,   0},
		{  0,   0,   0, c44,   0,   0},
		{  0,   0,   0,   0, c44,   0},
		{  0,   0,   0,   0,   0, c66}};

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					ret(i,j,k,l) = el_c [jum(i,j)][jum(k,l)];
				}
			}
		}
	}

	return ret;
}

Tensor
make_trigonal_tensor(	double c11,
						double c12,
						double c13,
						double c33,
						double c14,
						double c44,
						double c66) {
	Tensor ret;
	double el_c [6][6]={
		{c11, c12, c13, c14,   0,   0},
		{c12, c11, c13,-c14,   0,   0},
		{c13, c13, c33,   0,   0,   0},
		{c14,-c14,   0, c44,   0,   0},
		{  0,   0,   0,   0, c44, c14},
		{  0,   0,   0,   0, c14, c66}};

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					ret(i,j,k,l) = el_c [jum(i,j)][jum(k,l)];
				}
			}
		}
	}

	return ret;
}

Tensor3
make_trigonal_piezo_tensor_3m(	double e15,
								double e22,
								double e31,
								double e33) {
	Tensor3 ret;
	double el_p_c [3][6]={
		{   0,   0,   0,   0, e15, -e22},
		{-e22, e22,   0, e15,   0,    0},
		{ e31, e31, e33,   0,   0,    0}};

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				ret(i,j,k) = el_p_c [i][jum(j,k)];
			}
		}
	}

	return ret;
}

Tensor3
make_trigonal_piezo_tensor_32(	double e11,
								double e14) {
	Tensor3 ret;
	double el_p_c [3][6]={
		{e11,-e11, 0, e14,   0,   0},
		{  0,   0, 0,   0,-e14,-e11},
		{  0,   0, 0,   0,   0,   0}};

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				ret(i,j,k) = el_p_c [i][jum(j,k)];
			}
		}
	}

	return ret;
}

Matrix3 eps2mat (double e_xx, double e_zz) {
	Matrix3 ret;

	ret(0,0)=e_xx; ret(0,1)=0;	   ret(0,2)=0;
	ret(1,0)=0;    ret(1,1)= e_xx; ret(1,2)=0;
	ret(2,0)=0;	   ret(2,1)=0;	   ret(2,2)=e_zz;

	return ret;
}

Matrix3 rotz (double phi) {
	Matrix3 ret;

	ret(0,0)=cos(phi); ret(0,1)=-sin(phi); ret(0,2)=0;
	ret(1,0)=sin(phi); ret(1,1)= cos(phi); ret(1,2)=0;
	ret(2,0)=       0; ret(2,1)=        0; ret(2,2)=1;

	return ret;
}

Matrix3 rotx (double phi) {
	Matrix3 ret;

	ret(0,0)=1; ret(0,1)=  		0; ret(0,2)=  		0;
	ret(1,0)=0; ret(1,1)=cos(phi); ret(1,2)=-sin(phi);
	ret(2,0)=0; ret(2,1)=sin(phi); ret(2,2)= cos(phi);

	return ret;
}

Matrix3 roty (double phi) {
	Matrix3 ret;

	ret(0,0)=  cos(phi); ret(0,1)=0; ret(0,2)=sin(phi);
	ret(1,0)=	   	  0; ret(1,1)=1; ret(1,2)=		 0;
	ret(2,0)= -sin(phi); ret(2,1)=0; ret(2,2)=cos(phi);

	return ret;
}

Matrix3 turn_mat (double phi) {
    Matrix3 ret;

    ret(0,0)=  cos(phi); ret(0,1)= sin(phi); ret(0,2)= 0;
    ret(1,0)= -sin(phi); ret(1,1)= cos(phi); ret(1,2)= 0;
    ret(2,0)=         0; ret(2,1)=        0; ret(2,2)= 1;

    return ret;
}

Matrix3 euler(double alpha, double beta, double gamma){
return rotx(alpha)*roty(beta)*rotz(gamma);
}

VD rho2v(const VCD& arrg, double rho) {
	VD ret(arrg.size());

	for (int i = 0; i < (int) arrg.size(); i++) {
		ret[i]=sqrt(real(arrg[i])/rho);
	}

	return ret;
}

VD arr_re(const VCD& r) {
	VD ret(r.size());

	for (int i = 0; i < (int) r.size(); i++) {
		ret[i]=real(r[i]);
	}

	return ret;
}

Vector3 pq2vec(int p, int q, int n) {
Vector3 ret;
ret.z()=2.*p/n-1;
ret.y()=sqrt(1-ret.z()*ret.z())*sin(M_PI*q/n);
ret.x()=sqrt(1-ret.z()*ret.z())*cos(M_PI*q/n);
return ret;
}

Vector3
get_polarization(const Matrix3& g, double ev){
Matrix3 g1(g);
g1(0,0)-=ev;
g1(1,1)-=ev;
g1(2,2)-=ev;

Vector3 a(g1(0,0),g1(0,1),g1(0,2));
//cout << "a= " << a << endl;
Vector3 b(g1(1,0),g1(1,1),g1(1,2));
//cout << "b= " << b << endl;
Vector3 c(g1(2,0),g1(2,1),g1(2,2));
//cout << "c= " << c << endl;

Vector3 r1=a&b, r2=a&c, r3=b&c;
double m1=r1*r1, m2=r2*r2, m3=r3*r3;
if (m1 > m2 && m1 > m3) {
	return r1.normalized();
	}
if (m2 > m1 && m2 > m3) {
	return r2.normalized();
	}
return r3.normalized();
}


ostream&
operator<<(ostream& os, const VCD & r) {
	for(int i=0; i< (int) r.size(); i++) {
		os << r[i] << " ";
	}
	return os;
}

ostream&
operator<<(ostream& os, const VD & r) {
	for(int i=0; i< (int) r.size(); i++) {
		os << r[i] << " ";
	}
	return os;
}

Vector3
ort(const Vector3& n){
    Vector3 ret;
    Vector3 m;
    Vector3 n0=n.normalized();
    Vector3 n1(-n0.y(),n0.x(),n0.z());
    Vector3 n2(n0.x(),-n0.z(),n0.y());
    Vector3 n3(-n0.z(),n0.y(),n0.x());

    //cout << "n1=" << n1 << endl << "n2=" << n2 << endl << "n3=" << n3 << endl;

    double a=abs(n0*n1);
    double b=abs(n0*n2);
    double c=abs(n0*n3);

    //cout << "a=" << a << " b=" << b << " c=" << c << endl;

    if (b >= a && c >= a) {
        m = n1;
    } else if (a >= b && c >= b) {
        m = n2;
    } else {
        m = n3;
    }

    //cout << "m=" << m << endl;

    Vector3 k = m - n0 * (n0*m);

    ret=k.normalized();
    return ret;
}

double slowness(const Vector3& n,
                int ind,
                const Matrix3& epsilon,
                const Tensor& tens,
                const Tensor3& piezotens){

    Matrix3 tt = tens.make_piezo_crist(n, epsilon, piezotens);

    Poly pol = tt.characteristic();

    VCD comp_roots = pol.all_roots();
    vector<double> roots(comp_roots.size());
    for(int t = 0; t < int(comp_roots.size()); ++t)
            roots[t] = real(comp_roots[t]);
    sort(roots.begin(),roots.end());

    return sqrt(1/roots[ind]);
}

Vector3 slowness_vec(const Vector3& n,
                int ind,
                const Tensor& tens,
                double rho){


    Matrix3 cr=tens.crist(n);
    Poly pol= cr.characteristic();
    VCD roots = pol.all_roots();
    VD vels = rho2v(pol.all_roots(),rho);
    sort(vels.begin(),vels.end());
    return n*(1/vels[ind]);

}

Vector3 polaris(const Vector3& n,
                int ind,
                const Tensor& tens,
                double rho){

    Matrix3 cr=tens.crist(n);
    Poly pol = cr.characteristic();

    VCD comp_roots = pol.all_roots();
    vector<double> roots(comp_roots.size());
    for(int t = 0; t < int(comp_roots.size()); ++t)
            roots[t] = real(comp_roots[t]);
    sort(roots.begin(),roots.end());

    double gamma = roots[ind];
    return get_polarization( cr, gamma);
}

Vector3 slow_normal(const Vector3& n,
                    int ind,
                    const Matrix3& epsilon,
                    const Tensor& tens,
                    const Tensor3& piezotens){

    double delta=10E-4;

    Vector3 n0=n.normalized();
    Vector3 dn1=ort(n0)*delta;
    Vector3 dn2=n0&dn1;

    Vector3 n1=n0+dn1;
    Vector3 n2=n0+dn2;

    double s0 = slowness(n0,ind,epsilon,tens,piezotens);
    double s1 = slowness(n1,ind,epsilon,tens,piezotens);
    double s2 = slowness(n2,ind,epsilon,tens,piezotens);

    Vector3 sn0=n0*s0;
    Vector3 sn1=n1*s1;
    Vector3 sn2=n2*s2;

    Vector3 dsn1=sn1-sn0;
    Vector3 dsn2=sn2-sn0;

    return (dsn1&dsn2).normalized();
}

Matrix3 make_strain(const Vector3& q, const Vector3& s) {

    Matrix3 ret;

    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            ret(i,j)=(q(i)*s(j)+q(j)*s(i))/2;
        }
    }

    return ret;
}

//#pragma package(smart_init)

//---------------------------------------------------------------------------

#ifndef vector3H
#define vector3H
//---------------------------------------------------------------------------

#include <iostream>

class Vector3 {
	double dat[3];
public:
	Vector3(){
		dat[0]=0;
		dat[1]=0;
		dat[2]=0;
	}
	Vector3(double x, double y, double z){
		dat[0]=x;
		dat[1]=y;
		dat[2]=z;
	}
	double& operator()(int p) {return dat[p];}

	const double& operator()(int p) const {return dat[p];}

	std::ostream& operator>>(std::ostream& os)const;

	Vector3 operator&(const Vector3& ) const;

	double operator*(const Vector3& ) const;

	double& x() {return dat[0];}
	double& y() {return dat[1];}
	double& z() {return dat[2];}

	const double& x() const {return dat[0];}
	const double& y() const {return dat[1];}
	const double& z() const {return dat[2];}

	Vector3 normalized() const;

	Vector3& operator*=(double r) {
		dat[0] *= r;
		dat[1] *= r;
		dat[2] *= r;
		return *this;
	}

    Vector3 operator*(double r) const {
		Vector3 ret(*this);
		ret *= r;
		return ret;
	}

	Vector3 operator+(const Vector3& r) {
		return Vector3(x()+r.x(),y()+r.y(),z()+r.z());
	}


    Vector3 operator-(const Vector3& r) const;

};

std::ostream& operator<<(std::ostream& os, const Vector3& r);


#endif

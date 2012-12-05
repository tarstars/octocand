//---------------------------------------------------------------------------

//#pragma hdrstop

#include "vector3.h"

//---------------------------------------------------------------------------

#include <iostream>

using namespace std;

#include "matrix3.h"
#include "util.h"


  Vector3
  Vector3::operator&(const Vector3& r) const {
	return Vector3(	dat[1] * r.dat[2] - dat[2]* r.dat[1],
					dat[2] * r.dat[0] - dat[0]* r.dat[2],
					dat[0] * r.dat[1] - dat[1]* r.dat[0]);
  }


  double Vector3::operator*(const Vector3& r) const {
	return x()*r.x()+y()*r.y()+z()*r.z();
  }

  Vector3 Vector3::operator-(const Vector3& r) const {
      return Vector3 (x()-r.x(),y()-r.y(),z()-r.z());
  }

ostream&
Vector3::operator>>(ostream& os)const{
  return os << dat[0] << " " << dat[1] << " " << dat[2];
}

Vector3
Vector3::normalized() const{
//cout << "normalized: " << x() << " " << y() << " " << z() << endl;
double r=sqrt(x()*x()+y()*y()+z()*z());
return Vector3(x()/r,y()/r,z()/r);
}

ostream&
operator<<(ostream& os, const Vector3& r){
  return r >> os;
}



//#pragma package(smart_init)

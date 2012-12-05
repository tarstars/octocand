//---------------------------------------------------------------------------


//#pragma hdrstop

#include <cmath>
#include "pov_ray.h"
#include "vector3.h"

using namespace std;



//---------------------------------------------------------------------------

PovRay::PovRay(std::string flnm):dest(flnm.c_str()) {
dest << "\
background{color <1,1,1>}\n\
camera {\n\
  location <7.5,5.5,5>\n\
  look_at <0,0,0>\n\
  angle 20\n\
  }\n\
light_source {<1,1,10> color <0.8,0.8,0.8>} \n\
light_source {<10,10,-10> color <0.8,0.8,0.8>} \n\
light_source {<-10,10,-10> color <0.8,0.8,0.8>} \n\
light_source {<-10,10,10> color <0.8,0.8,0.8>} \n\
\n\
cylinder{ <0,0,0>,<1.5,0,0>, 0.01  pigment {color <1,0,0>}} \n\
cylinder{ <0,0,0>,<0,0,1.5>, 0.01  pigment {color <0,1,0>}} \n\
cylinder{ <0,0,0>,<0,1.5,0>, 0.01  pigment {color <0,0,1>}} \n\
\n\
cylinder{ <1.45,0.15,0>,<1.35,0.05,0>, 0.01  pigment {color <1,0,0>}} \n\
cylinder{ <1.45,0.05,0>,<1.35,0.15,0>, 0.01  pigment {color <1,0,0>}} \n\
\n\
cylinder{ <0,0.1,1.4>,<0,0.15,1.45>, 0.01  pigment {color <0,1,0>}} \n\
cylinder{ <0,0.15,1.35>,<0,0.1,1.4>, 0.01  pigment {color <0,1,0>}} \n\
cylinder{ <0,0.1,1.4>,<0,0.05,1.4>, 0.01  pigment {color <0,1,0>}} \n\
\n\
cylinder{ <0.15,1.45,0>,<0.05,1.45,0>, 0.01  pigment {color <0,0,1>}} \n\
cylinder{ <0.15,1.35,0>,<0.05,1.45,0>, 0.01  pigment {color <0,0,1>}} \n\
cylinder{ <0.15,1.35,0>,<0.05,1.35,0>, 0.01  pigment {color <0,0,1>}} \n\
\n";

/*
cylinder{ <0,0,0>,<1.5,0,0>, 0.01  pigment {color <0.9,0.2,0>}} \n\
cylinder{ <0,0,0>,<0,0,1.5>, 0.01  pigment {color <0.9,0.2,0>}} \n\
cylinder{ <0,0,0>,<0,1.5,0>, 0.01  pigment {color <0.9,0.2,0>}} \n\
\n\
cylinder{ <1.45,0.15,0>,<1.35,0.05,0>, 0.01  pigment {color <0.9,0.2,0>}} \n\
cylinder{ <1.45,0.05,0>,<1.35,0.15,0>, 0.01  pigment {color <0.9,0.2,0>}} \n\
\n\
cylinder{ <0,0.1,1.4>,<0,0.15,1.45>, 0.01  pigment {color <0.9,0.2,0>}} \n\
cylinder{ <0,0.15,1.35>,<0,0.1,1.4>, 0.01  pigment {color <0.9,0.2,0>}} \n\
cylinder{ <0,0.1,1.4>,<0,0.05,1.4>, 0.01  pigment {color <0.9,0.2,0>}} \n\
\n\
cylinder{ <0.15,1.45,0>,<0.05,1.45,0>, 0.01  pigment {color <0.9,0.2,0>}} \n\
cylinder{ <0.15,1.35,0>,<0.05,1.45,0>, 0.01  pigment {color <0.9,0.2,0>}} \n\
cylinder{ <0.15,1.35,0>,<0.05,1.35,0>, 0.01  pigment {color <0.9,0.2,0>}} \n\
\n\
sphere{ <0,0,0>, 1 \n\
            pigment {color <0.8,0.8,0> filter 0.9}\n\
            finish { \n\
                    ambient 0 \n\
                    diffuse 0.5 \n\
                    reflection .25 \n\
                    specular 1 \n\
                    roughness .01 \n\
    }}
";*/

}

double nuller(double x){
	if (abs(x)<1e-6)
		return 0;
	return x;
}

void PovRay::point(const Vector3& r, double rad) {
    dest << "sphere{ <" << nuller(r.x()) << "," << nuller(r.z()) << "," << nuller(r.y()) << ">, " << rad << " pigment {color <0.7,0.7,0>}}" << endl;
}

void PovRay::cylinder(const Vector3& rn, const Vector3& drn, double radc, const Vector3& rgb) {
    dest << "cylinder{ <" << rn.x() << "," << rn.z() << "," << rn.y() << ">, <" << drn.x() << "," << drn.z() << "," << drn.y() << "> ," << radc << " pigment {color <" << rgb << ">}}" << endl;
}

void PovRay::triangle(const Vector3& r1, const Vector3& r2, const Vector3& r3) {
dest << "triangle { <"
	<< nuller(r1.x()) << "," << nuller(r1.z()) << "," << nuller(r1.y()) << ">, <"
	<< nuller(r2.x()) << "," << nuller(r2.z()) << "," << nuller(r2.y()) << ">, <"
	<< nuller(r3.x()) << "," << nuller(r3.z()) << "," << nuller(r3.y()) <<
	">}" << endl;
}

void PovRay::begin_mesh() {
dest << "mesh {" << endl;
}

void PovRay::end_mesh() {
dest << "pigment {color <0,1,0> filter 0.5}\n\
		finish { \n\
				ambient 0 \n\
				diffuse 0.5 \n\
				reflection .25 \n\
				specular 1 \n\
				roughness .01 \n\
}}" << endl;
}

//#pragma package(smart_init)

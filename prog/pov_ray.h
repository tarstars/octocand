//---------------------------------------------------------------------------

#ifndef pov_rayH
#define pov_rayH
//---------------------------------------------------------------------------

#include <fstream>
#include <string>

#include "vector3.h"

class PovRay {
std::ofstream dest;

public:
	PovRay(std::string);

    void point(const Vector3&, double rad=0.005);
    void cylinder(const Vector3&, const Vector3&, double radc=0.005, const Vector3& rgb=Vector3(1,0,0));

	void triangle(const Vector3& r1, const Vector3& r2, const Vector3& r3);

	void begin_mesh();
	void end_mesh();

};

#endif

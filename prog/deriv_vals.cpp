//---------------------------------------------------------------------------


#pragma hdrstop

#include "deriv_vals.h"

#include <iostream>

using namespace std;

//---------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, const DerivVals& r){
	return r >> os;
}

#pragma package(smart_init)

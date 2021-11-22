#include<iostream>
#include<fstream>
#include<math.h>
#include"Helicities_ZZ.h"
using namespace std;

namespace sm_aazz {

    
//const double alpha_em = 1./137.036; // EM coupling at zero momentum (on shell scheme)


// compute the SM squared matrix element, including leptons, quarks and the W boson
double sqme(double s, double t, int exclude_loops){

  return dsigma_ZZ(s,t,exclude_loops);  

}

}

#include<math.h>
#include<iostream>
#include<fstream>

//#include "functions.h"



const double Pi=4*atan(1);
const double mW2=80.4*80.4+0.0001;// everything is up to power 2
const double mZ2=91.2*91.2+0.0001;
const double mh2=125*125+0.0001;
const double cW=sqrt(0.769);
const double sW=sqrt(0.231);
const double alpha=1/137.036;
const double Gammah=0.01;

double mup2=0.0023*0.0023;
double mdown2=0.0048*0.0048;
double mstrange2=0.095*0.095;
double mcharm2=1.275*1.275;
double mbottom2=4.18*4.18;
double mtop2=173.07*173.07;
double melectron2=0.000511*0.000511;
double mmuon2=0.1057*0.1057;
double mtau2=1.777*1.777;

double Qup=-1/3.;
double Qdown=2/3.;
double Qstrange=-1/3.;
double Qcharm=2/3.;
double Qbottom=-1/3.;
double Qtop=2/3.;
double Qelectron=-1;
double Qmuon=-1;
double Qtau=-1;

double Iup=1/2.;
double Idown=-1/2.;
double Istrange=-1/2.;
double Icharm=1/2.;
double Ibottom=-1/2.;
double Itop=1/2.;
double Ielectron=-1/2.;
double Imuon=-1/2.;
double Itau=-1/2.;

double reAs_pppm(double beta, double t, double u, double m)
{
     return -4;
}

double imAs_pppm(double beta, double t, double u, double m)
{
     return 0;
}

double reAs_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return 4-4*u*t/(s*s)*(pow(log(t/u),2)+Pi*Pi)+4*(t-u)/s*log(t/u);
}

double imAs_pppp(double beta, double t, double u, double m)
{
     return 0;
}

double reAs_pmpp(double beta, double t, double u, double m)
{
     
     return -4;
}
double imAs_pmpp(double beta, double t, double u, double m)
{
     
     return 0;
}

double reAs_pm00(double beta, double t, double u, double m)
{
     
     return 0;
}
double imAs_pm00(double beta, double t, double u, double m)
{
     
     return 0;
}



double reAs_pp00(double beta, double t, double u, double m)
{
     
     return 0;
}
double imAs_pp00(double beta, double t, double u, double m)
{
     
     return 0;
}


double reAs_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2*u*t),0.5);
     
     return delta*(-8*(t-u)/s+4*(t-u)*t*u/pow(s,3)*(pow(log(t/u),2)+Pi*Pi)+16*u*t/pow(s,2)*log(t/u));
}

double imAs_ppp0(double beta, double t, double u, double m)
{
     return 0;
}



double reAs_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2*u*t),0.5);
     
     return delta*(-8*u/s+4*t/u*pow(log(-s/t),2)+8*t/s*log(-s/t));
}

double imAs_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2*u*t),0.5);
     
     return delta*(4*t/u*(-2*Pi*log(-s/t))-8*t*Pi/s);
}


double reAs_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
 
     return 4-4*s*t/(u*u)*(pow(log(-s/t),2))+4*(s-t)/u*log(-s/t);
}

double imAs_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -4*s*t/(u*u)*(-2*Pi*log(-s/t))+4*(s-t)/u*(-Pi);
}




//W delat 
double re_dw_pppm(double beta, double t, double u, double m)
{
     return 0;
}
double im_dw_pppm(double beta, double t, double u, double m)
{
     return 0;
}

double re_dw_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double red=16*cW*cW*(pow(log(t/u),2)+Pi*Pi+s/u*log(-u/m)*log(-s/t)+s/t*log(-t/m)*log(-s/u));
     return red;
}
double im_dw_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double imd=16*cW*cW*(-Pi*(s/u*log(-u/m)+s/t*log(-t/m)));
     return imd;
}

double re_dw_pmpp(double beta, double t, double u, double m)
{
     
     return 0;
}

double im_dw_pmpp(double beta, double t, double u, double m)
{
     
     return 0;
}

double re_dw_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double red= -2*(pow(log(t/u),2)+s/u*pow(log(-t/s),2)+s/t*pow(log(-u/s),2))-4/s*log(s/m)*(t*log(-t/s)+u*log(-u/s));
     return red;
}

double im_dw_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double imd= 4*Pi*(log(s/m)+(1-t*t/(s*u))*log(-t/s)+(1-u*u/(s*t))*log(-u/s));
     return imd;
}

double re_dw_pp00(double beta, double t, double u, double m)
{    
     double s=4*mZ2/(1.0-beta*beta);
     
     return 2*pow(log(s/m),2)-2*Pi*Pi;
}

double im_dw_pp00(double beta, double t, double u, double m)
{    
     double s=4*mZ2/(1-beta*beta);
     
     return  -4*Pi*log(s/m);
}

double re_dw_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     
     double red=delta*(-2*(t-u)/s*(pow(log(s/m),2)+pow(log(t/u),2))-4*log(s/m)*log(t/u));
     //printf(" %lf\n",red);
     return red;
}

double im_dw_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     double imd= delta*(-2*(t-u)/s*(-2*Pi*log(s/m))+4*Pi*log(t/u));;
     return imd;
}

double re_dw_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     double red= delta*((-2*u/s+32*cW*cW*t/s)*(pow(log(s/m),2))+32*cW*cW*(u/s*log(-t/m)*log(-u/m)+t/s*pow(log(-t/m),2))-2*u/s*pow(log(t/u),2)+(-4*u*(t-u)/(s*s)+32*cW*cW*(t*t+s*s)/(s*s))*(log(s/m))*log(-t/m)-4*u/(s*s)*(u-t-8*cW*cW*t)*log(s/m)*log(-u/m));
     return red;
}

double im_dw_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     double imd= delta*((-2*u/s+32*cW*cW*t/s)*(-2*Pi*log(s/m))+(-4*u*(t-u)/(s*s)+32*cW*cW*(t*t+s*s)/(s*s))*(-Pi)*log(-t/m)-4*u/(s*s)*(u-t-8*cW*cW*t)*(-Pi)*log(-u/m));
     return imd;
}

double re_dw_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double red= 16*cW*cW*(pow(log(-t/s),2)+u/t*log(-t/m)*log(-u/s)+u/s*log(s/m)*log(u/t));
     return red;
}

double im_dw_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double imd= 16*cW*cW*Pi*(u/t*log(s/m)-u/s*log(-u/s)-(s*s+t*t)/(s*t)*log(-t/s));
     return imd;
}

//W amplitudes

double reAw_pppm(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pppm(beta,t,u,m)+re_dw_pppm(beta,t,u,m);
}

double imAw_pppm(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pppm(beta,t,u,m)+im_dw_pppm(beta,t,u,m);
}

double reAw_pppp(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pppp(beta,t,u,m)+re_dw_pppp(beta,t,u,m);
}

double imAw_pppp(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pppp(beta,t,u,m)+im_dw_pppp(beta,t,u,m);
}

double reAw_pmpp(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pmpp(beta,t,u,m)+re_dw_pmpp(beta,t,u,m);
}

double imAw_pmpp(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pmpp(beta,t,u,m)+im_dw_pmpp(beta,t,u,m);
}

double reAw_pm00(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pm00(beta,t,u,m)+re_dw_pm00(beta,t,u,m);
}

double imAw_pm00(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pm00(beta,t,u,m)+im_dw_pm00(beta,t,u,m);
}

double reAw_pp00(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pp00(beta,t,u,m)+re_dw_pp00(beta,t,u,m);
}

double imAw_pp00(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pp00(beta,t,u,m)+im_dw_pp00(beta,t,u,m);
}

double reAw_ppp0(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_ppp0(beta,t,u,m)+re_dw_ppp0(beta,t,u,m);
}

double imAw_ppp0(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_ppp0(beta,t,u,m)+im_dw_ppp0(beta,t,u,m);
}

double reAw_pmp0(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pmp0(beta,t,u,m)+re_dw_pmp0(beta,t,u,m);
}

double imAw_pmp0(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pmp0(beta,t,u,m)+im_dw_pmp0(beta,t,u,m);
}

double reAw_pmpm(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pmpm(beta,t,u,m)+re_dw_pmpm(beta,t,u,m);
}

double imAw_pmpm(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pmpm(beta,t,u,m)+im_dw_pmpm(beta,t,u,m);
}





//
//amplitude of helicity with prefactors for scalar
//

double reFs_pppm(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_pppm(beta,t,u,m);
}

double imFs_pppm(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_pppm(beta,t,u,m);
}

double reFs_pppp(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_pppp(beta,t,u,m);
}

double imFs_pppp(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_pppp(beta,t,u,m);
}

double reFs_pmpp(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_pmpp(beta,t,u,m);
}

double imFs_pmpp(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_pmpp(beta,t,u,m);
}

double reFs_pm00(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_pm00(beta,t,u,m);
}

double imFs_pm00(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_pm00(beta,t,u,m);
}

double reFs_pp00(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_pp00(beta,t,u,m);
}

double imFs_pp00(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_pp00(beta,t,u,m);
}

double reFs_ppp0(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_ppp0(beta,t,u,m);
}

double imFs_ppp0(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_ppp0(beta,t,u,m);
}

double reFs_pmp0(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_pmp0(beta,t,u,m);
}

double imFs_pmp0(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_pmp0(beta,t,u,m);
}

double reFs_pmpm(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*reAs_pmpm(beta,t,u,m);
}

double imFs_pmpm(double beta, double t, double isospin,double charge, double m)
{
     double g= (isospin-charge*sW*sW)/(sW*cW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*g*g*imAs_pmpm(beta,t,u,m);
}


//
//amplitude of helicity with prefactors for W boson
//

double reFw_pppm(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pppm(beta,t,u,m);
}

double imFw_pppm(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pppm(beta,t,u,m);
}

double reFw_pppp(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pppp(beta,t,u,m);
}

double imFw_pppp(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pppp(beta,t,u,m);
}

double reFw_pmpp(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pmpp(beta,t,u,m);
}

double imFw_pmpp(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pmpp(beta,t,u,m);
}

double reFw_pm00(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pm00(beta,t,u,m);
}

double imFw_pm00(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     //printf(" allons y f*gaimenrt %lf\n",sW*sW);
     //printf(" allons y f*gaimenrt %lf\n",im_dw_pm00(beta,t,u,m)/(sW*sW));
     return alpha*alpha/(sW*sW)*imAw_pm00(beta,t,u,m);
}

double reFw_pp00(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     //printf(" allons y f*gaimenrt %lf\n",reAs_pp00(beta,t,u,m));
     
     return alpha*alpha/(sW*sW)*reAw_pp00(beta,t,u,m);
}

double imFw_pp00(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pp00(beta,t,u,m);
}

double reFw_ppp0(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_ppp0(beta,t,u,m);
}

double imFw_ppp0(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_ppp0(beta,t,u,m);
}

double reFw_pmp0(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pmp0(beta,t,u,m);
}

double imFw_pmp0(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pmp0(beta,t,u,m);
}

double reFw_pmpm(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pmpm(beta,t,u,m);
}

double imFw_pmpm(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pmpm(beta,t,u,m);
}




//delta fermion 

double re_dvf_pppm(double beta, double t, double u, double m)
{
     return 0;
}

double im_dvf_pppm(double beta, double t, double u, double m)
{
     return 0;
}
double re_daf_pppm(double beta, double t, double u, double m)
{
     return 0;
}

double im_daf_pppm(double beta, double t, double u, double m)
{
     return 0;
}

double re_dvf_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -4*(pow(log(u/t),2)+Pi*Pi);
}

double im_dvf_pppp(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -4*(pow(log(u/t),2)+Pi*Pi);
}

double im_daf_pppp(double beta, double t, double u, double m)
{
     return 0;
}

double re_dvf_pmpp(double beta, double t, double u, double m)
{
     
     return 0;
}

double im_dvf_pmpp(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_pmpp(double beta, double t, double u, double m)
{ 
     return 0;
}

double im_daf_pmpp(double beta, double t, double u, double m)
{
     return 0;
}

double re_dvf_pp00(double beta, double t, double u, double m)
{
     
     return 0;
}

double im_dvf_pp00(double beta, double t, double u, double m)
{  
     
     return 0;
}

double re_daf_pp00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -8*m/mZ2*(pow(log(s/m),2)-Pi*Pi);
}

double im_daf_pp00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -8*m/mZ2*(-2*Pi*log(s/m));
}

double re_dvf_pm00(double beta, double t, double u, double m)
{
     return 0;
}

double im_dvf_pm00(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -8*m/mZ2*(t/u*pow(log(-t/s),2)+u/t*pow(log(-u/s),2));
}

double im_daf_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -8*m/mZ2*(t/u*2*Pi*log(-t/s)+u/t*2*Pi*log(-u/s));
}

double re_dvf_ppp0(double beta, double t, double u, double m)
{
     return 0;
}

double im_dvf_ppp0(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -pow(s*mZ2/(2*u*t),0.5)*8*m/mZ2*((u-t)/s*(pow(log(s/m),2)+pow(log(t/u),2))-2*log(s/m)*log(t/u));
}

double im_daf_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -pow(s*mZ2/(2*u*t),0.5)*8*m/mZ2*((u-t)/s*(-2*Pi*log(s/m))+2*Pi*log(t/u));
}

double re_dvf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(pow(log(-s/t),2));
}

double im_dvf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(-2*Pi*log(-s/t));
}
double re_daf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(pow(log(-s/t),2))+pow(s*mZ2/(2*u*t),0.5)*8*m*u/(mZ2*s)*(pow(log(s*m/(t*u)),2));
}

double im_daf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(-2*Pi*log(-s/t))+pow(s*mZ2/(2*u*t),0.5)*8*m*u/(mZ2*s)*(-2*Pi*log(s*m/(t*u)));
}

double re_dvf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -4*pow(log(-s/t),2);
}

double im_dvf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return 8*Pi*log(-s/t);
}
double re_daf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -4*pow(log(-s/t),2);;
}

double im_daf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return 8*Pi*log(-s/t);;
}



double reAvf_pppm(double beta, double t, double u, double m)
{
     return -2*reAs_pppm(beta,t,u,m)+re_dvf_pppm(beta,t,u,m);
}

double imAvf_pppm(double beta, double t, double u, double m)
{
     return -2*imAs_pppm(beta,t,u,m)+im_dvf_pppm(beta,t,u,m);
}

double reAvf_pppp(double beta, double t, double u, double m)
{
     //printf("ok ici 1 \n");
     return -2*reAs_pppp(beta,t,u,m)+re_dvf_pppp(beta,t,u,m);
}

double imAvf_pppp(double beta, double t, double u, double m)
{
     return -2*imAs_pppp(beta,t,u,m)+im_dvf_pppp(beta,t,u,m);
}

double reAvf_pmpp(double beta, double t, double u, double m)
{
     return -2*reAs_pmpp(beta,t,u,m)+re_dvf_pmpp(beta,t,u,m);
}

double imAvf_pmpp(double beta, double t, double u, double m)
{
     return  -2*imAs_pmpp(beta,t,u,m)+im_dvf_pmpp(beta,t,u,m);
}

double reAvf_pm00(double beta, double t, double u, double m)
{
     return -2*reAs_pm00(beta,t,u,m)+re_dvf_pm00(beta,t,u,m);
}

double imAvf_pm00(double beta, double t, double u, double m)
{
     return -2*imAs_pm00(beta,t,u,m)+im_dvf_pm00(beta,t,u,m);
}

double reAvf_pp00(double beta, double t, double u, double m)
{
     return -2*reAs_pp00(beta,t,u,m)+re_dvf_pp00(beta,t,u,m);
}

double imAvf_pp00(double beta, double t, double u, double m)
{
     return -2*imAs_pp00(beta,t,u,m)+im_dvf_pp00(beta,t,u,m);
}

double reAvf_ppp0(double beta, double t, double u, double m)
{
     return -2*reAs_ppp0(beta,t,u,m)+re_dvf_ppp0(beta,t,u,m);
}

double imAvf_ppp0(double beta, double t, double u, double m)
{
     return -2*imAs_ppp0(beta,t,u,m)+im_dvf_ppp0(beta,t,u,m);
}

double reAvf_pmp0(double beta, double t, double u, double m)
{
     return -2*reAs_pmp0(beta,t,u,m)+re_dvf_pmp0(beta,t,u,m);
}

double imAvf_pmp0(double beta, double t, double u, double m)
{
     return -2*imAs_pmp0(beta,t,u,m)+im_dvf_pmp0(beta,t,u,m);
}

double reAvf_pmpm(double beta, double t, double u, double m)
{
     return -2*reAs_pmpm(beta,t,u,m)+re_dvf_pmpm(beta,t,u,m);
}

double imAvf_pmpm(double beta, double t, double u, double m)
{
     return -2*imAs_pmpm(beta,t,u,m)+im_dvf_pmpm(beta,t,u,m);
}



double reAaf_pppm(double beta, double t, double u, double m)
{
     return -2*reAs_pppm(beta,t,u,m)+re_daf_pppm(beta,t,u,m);
}

double imAaf_pppm(double beta, double t, double u, double m)
{
     return -2*imAs_pppm(beta,t,u,m)+im_daf_pppm(beta,t,u,m);
}

double reAaf_pppp(double beta, double t, double u, double m)
{
     return -2*reAs_pppp(beta,t,u,m)+re_daf_pppp(beta,t,u,m);
}

double imAaf_pppp(double beta, double t, double u, double m)
{
     return -2*imAs_pppp(beta,t,u,m)+im_daf_pppp(beta,t,u,m);
}

double reAaf_pmpp(double beta, double t, double u, double m)
{
     return -2*reAs_pmpp(beta,t,u,m)+re_daf_pmpp(beta,t,u,m);
}

double imAaf_pmpp(double beta, double t, double u, double m)
{
     return  -2*imAs_pmpp(beta,t,u,m)+im_daf_pmpp(beta,t,u,m);
}

double reAaf_pm00(double beta, double t, double u, double m)
{
     return -2*reAs_pm00(beta,t,u,m)+re_daf_pm00(beta,t,u,m);
}

double imAaf_pm00(double beta, double t, double u, double m)
{
     return -2*imAs_pm00(beta,t,u,m)+im_daf_pm00(beta,t,u,m);
}

double reAaf_pp00(double beta, double t, double u, double m)
{
     return -2*reAs_pp00(beta,t,u,m)+re_daf_pp00(beta,t,u,m);
}

double imAaf_pp00(double beta, double t, double u, double m)
{
     return -2*imAs_pp00(beta,t,u,m)+im_daf_pp00(beta,t,u,m);
}

double reAaf_ppp0(double beta, double t, double u, double m)
{
     return -2*reAs_ppp0(beta,t,u,m)+re_daf_ppp0(beta,t,u,m);
}

double imAaf_ppp0(double beta, double t, double u, double m)
{
     return -2*imAs_ppp0(beta,t,u,m)+im_daf_ppp0(beta,t,u,m);
}

double reAaf_pmp0(double beta, double t, double u, double m)
{
     return -2*reAs_pmp0(beta,t,u,m)+re_daf_pmp0(beta,t,u,m);
}

double imAaf_pmp0(double beta, double t, double u, double m)
{
     return -2*imAs_pmp0(beta,t,u,m)+im_daf_pmp0(beta,t,u,m);
}

double reAaf_pmpm(double beta, double t, double u, double m)
{
     return -2*reAs_pmpm(beta,t,u,m)+re_daf_pmpm(beta,t,u,m);
}

double imAaf_pmpm(double beta, double t, double u, double m)
{
     return -2*imAs_pmpm(beta,t,u,m)+im_daf_pmpm(beta,t,u,m);
}




double reFf_pppm(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_pppm(beta,t,u,m)+gv*gv*reAvf_pppm(beta,t,u,m));
}

double imFf_pppm(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_pppm(beta,t,u,m)+gv*gv*imAvf_pppm(beta,t,u,m));
}

double reFf_pppp(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_pppp(beta,t,u,m)+gv*gv*reAvf_pppp(beta,t,u,m));
}

double imFf_pppp(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_pppp(beta,t,u,m)+gv*gv*imAvf_pppp(beta,t,u,m));
}

double reFf_pmpp(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_pmpp(beta,t,u,m)+gv*gv*reAvf_pmpp(beta,t,u,m));
}

double imFf_pmpp(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_pmpp(beta,t,u,m)+gv*gv*imAvf_pmpp(beta,t,u,m));
}

double reFf_pm00(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_pm00(beta,t,u,m)+gv*gv*reAvf_pm00(beta,t,u,m));
}

double imFf_pm00(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_pm00(beta,t,u,m)+gv*gv*imAvf_pm00(beta,t,u,m));
}

double reFf_pp00(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_pp00(beta,t,u,m)+gv*gv*reAvf_pp00(beta,t,u,m));
}

double imFf_pp00(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_pp00(beta,t,u,m)+gv*gv*imAvf_pp00(beta,t,u,m));
}

double reFf_ppp0(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_ppp0(beta,t,u,m)+gv*gv*reAvf_ppp0(beta,t,u,m));
}

double imFf_ppp0(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_ppp0(beta,t,u,m)+gv*gv*imAvf_ppp0(beta,t,u,m));
}

double reFf_pmp0(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_pmp0(beta,t,u,m)+gv*gv*reAvf_pmp0(beta,t,u,m));
}

double imFf_pmp0(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_pmp0(beta,t,u,m)+gv*gv*imAvf_pmp0(beta,t,u,m));
}

double reFf_pmpm(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*reAaf_pmpm(beta,t,u,m)+gv*gv*reAvf_pmpm(beta,t,u,m));
}

double imFf_pmpm(double beta, double t, double isospin,double charge, double m)
{
     double gv= (isospin-2*charge*sW*sW)/(2*sW*cW);
     double ga=isospin/(2*cW*sW);
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha*charge*charge*(ga*ga*imAaf_pmpm(beta,t,u,m)+gv*gv*imAvf_pmpm(beta,t,u,m));
}


//8 light fermions
const double fermions_m2 [8]={mup2,mdown2,mstrange2,mcharm2,mbottom2,melectron2,mmuon2,mtau2};
const double fermions_Q [8]={Qup,Qdown,Qstrange,Qcharm,Qbottom,Qelectron,Qmuon,Qtau};
const double fermions_I [8]={Iup,Idown,Istrange,Icharm,Ibottom,Ielectron,Imuon,Itau};

double reFf_pppm_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pppm_as(beta,t,u,m2)+gv*gv*reAvf_pppm_as(beta,t,u,m2));
     }
     return sum;
}

double imFf_pppm_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pppm_as(beta,t,u,m2)+gv*gv*imAvf_pppm_as(beta,t,u,m2));
     }
     return sum;
}

double reFf_pppp_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pppp_as(beta,t,u,m2)+gv*gv*reAvf_pppp_as(beta,t,u,m2));
     }
     return sum;
}

double imFf_pppp_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pppp_as(beta,t,u,m2)+gv*gv*imAvf_pppp_as(beta,t,u,m2));
     }
     return sum;
}

double reFf_pmpp_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pmpp_as(beta,t,u,m2)+gv*gv*reAvf_pmpp_as(beta,t,u,m2));
     }
     return sum;
}
double imFf_pmpp_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pmpp_as(beta,t,u,m2)+gv*gv*imAvf_pmpp_as(beta,t,u,m2));
     }
     return sum;
}

double reFf_pm00_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pm00_as(beta,t,u,m2)+gv*gv*reAvf_pm00_as(beta,t,u,m2));
     }
     return sum;
}

double imFf_pm00_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pm00_as(beta,t,u,m2)+gv*gv*imAvf_pm00_as(beta,t,u,m2));
     }
     return sum;
}
double reFf_pp00_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pp00_as(beta,t,u,m2)+gv*gv*reAvf_pp00_as(beta,t,u,m2));
     }
     return sum;
}

double imFf_pp00_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pp00_as(beta,t,u,m2)+gv*gv*imAvf_pp00_as(beta,t,u,m2));
     }
     return sum;
}

double reFf_ppp0_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_ppp0_as(beta,t,u,m2)+gv*gv*reAvf_ppp0_as(beta,t,u,m2));
     }
     return sum;
}

double imFf_ppp0_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_ppp0_as(beta,t,u,m2)+gv*gv*imAvf_ppp0_as(beta,t,u,m2));
     }
     return sum;
}

double reFf_pmp0_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pmp0_as(beta,t,u,m2)+gv*gv*reAvf_pmp0_as(beta,t,u,m2));
     }
     return sum;
}

double imFf_pmp0_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pmp0_as(beta,t,u,m2)+gv*gv*imAvf_pmp0_as(beta,t,u,m2));
     }
     return sum;
}

double reFf_pmpm_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pmpm_as(beta,t,u,m2)+gv*gv*reAvf_pmpm_as(beta,t,u,m2));
     }
     return sum;
}

double imFf_pmpm_as(double beta, double t)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     double isospin;
     double m2;
     double charge;

     double ga;
     double gv;
     double sum=0;
     for (int i=0;i<8;i++){
         isospin=fermions_I[i];
         charge=fermions_Q[i];
         m2=fermions_m2[i]; 

         ga=isospin/(2*cW*sW);
         gv=(isospin-2*charge*sW*sW)/(2*sW*cW);
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pmpm_as(beta,t,u,m2)+gv*gv*imAvf_pmpm_as(beta,t,u,m2));
     }
     return sum;
}



//Higgs
double ReHi(double s, double Qi, double Nci, double spin, double mi)
{
     double fi=0;
     double st=s/(4.*mi);
     if (spin==0)
     {
         //fi=1/st*(1+2*mi*ReC(st,mi));
     }
     else if (spin==0.5)
     {
         //fi=-2/st*(1-2*mi*(st-1)*ReC(st,mi));
     }
     else if (spin==1)
     {
         //fi=1/st*(3+mh2/(2.*mi))-(16*mi-(mh2+6*mi)/st)*ReC(st,mi);
     }
     else 
     {
         printf("error in the higgs channel spin not allowed");
     }
     
     return Nci*Qi*Qi*fi;
}
    
double ImHi(double s, double Qi, double Nci, double spin, double mi)
{
     double fi=0;
     double st=s/(4*mi);
     if (spin==0)
     {
         //fi=1/st*(2*mi*ImC(st,mi));
     }
     else if (spin==0.5)
     {
         //fi=-2/st*(-2*mi*(st-1)*ImC(st,mi));
     }
     else if (spin==1)
     {
         //fi=-(16*mi-(mh2+6*mi)/st)*ImC(st,mi);
     }
     else 
     {
         printf("error in the higgs channel spin not allowed");
     }
     
     return Nci*Qi*Qi*fi;
}

double sumreHi(double s)//here we  sum over the SM contribution
{
     double wbos=1;// if 0 do not take into account 1 yes
     double tquark=1;
     double uquark=0;
     double tau=1;
     double sum=0;
     
     if (wbos==1)
     {
         sum+=ReHi(s,+1,1,1,mW2);//+ReHi(s,-1,1,1,mW2);
     }
     if (tquark==1)
     {
         sum+=ReHi(s,-2/3.,3,1/2.,173.34*173.34);
     }
     if (tau==1)
     {
         sum+=ReHi(s,-1.,1,1/2.,1.77*1.77);
     }
     if (uquark==1)
     {
         sum+=ReHi(s,-2/3.,3,1/2.,0.0022*0.0022);
     }

     return sum;
}
double sumimHi(double s)//here we  sum over the SM contribution
{
     double wbos=1;// if 0 do not take into account 1 yes
     double tquark=1;
     double uquark=0;
     double tau=1;
     double sum=0;
     
     if (wbos==1)
     {
         sum+=ImHi(s,+1,1,1,mW2);//+ImHi(s,-1,1,1,mW2);
     }
     if (tquark==1)
     {
         sum+=ImHi(s,-2/3.,3,1/2.,173.34*173.34);
     }
     if (tau==1)
     {
         sum+=ImHi(s,-1.,1,1/2.,1.77*1.77);
     }
     if (uquark==1)
     {
         sum+=ImHi(s,-2/3.,3,1/2.,0.0022*0.0022);
     }
     return sum;
}


double reFh(double beta,double l1,double l2,double l3,double l4)//,double l1,double l2,double l3,double l4)//here we  sum over the SM contribution
{
     double s=4*mZ2/(1-beta*beta);
     double hel=(1+l1*l2)/2.*((1+l3*l4)*l3*l4/2.-(1+beta*beta)/(1-beta*beta)*(1-l3*l3)*(1-l4*l4));
     
     return -alpha*alpha/(2*sW*sW*cW*cW)*s/((s-mh2)*(s-mh2)+mh2*Gammah*Gammah)*((s-mh2)*sumreHi(s)+pow(mh2,0.5)*Gammah*sumimHi(s))*hel;
}
double imFh(double beta,double l1,double l2,double l3,double l4)//,double l1,double l2,double l3,double l4)//here we  sum over the SM contribution
{
     double s=4*mZ2/(1-beta*beta);
     double hel=(1+l1*l2)/2.*((1+l3*l4)*l3*l4/2.-(1+beta*beta)/(1-beta*beta)*(1-l3*l3)*(1-l4*l4));
     //printf(" hel %lf\n",hel/s);
     return -alpha*alpha/(2*sW*sW*cW*cW)*s/((s-mh2)*(s-mh2)+mh2*Gammah*Gammah)*((s-mh2)*sumimHi(s)-pow(mh2,0.5)*Gammah*sumreHi(s))*hel;
}

double reFh_pppp(double beta)//here we  sum over the SM contribution
{
     return reFh(beta,1,1,1,1);
}
double imFh_pppp(double beta)//here we  sum over the SM contribution
{
     return imFh(beta,1,1,1,1);
}

double reFh_pp00(double beta)//here we  sum over the SM contribution
{
     return reFh(beta,1,1,0,0);
}
double imFh_pp00(double beta)//here we  sum over the SM contribution
{
     return imFh(beta,1,1,0,0);
}



double dsigma_ZZ(double s, double t,int exclude_loop)
{
     if (s/(4*mZ2)<1)
     {
          return 0;
     }
     
     if (s<4*mZ2)
     {
          return 0;
     }
     double beta=pow(1-4*mZ2/(s+0.00001),0.5);
     if (t>mZ2-s/2.*(1-beta))          //mZ2-s/2.*(1-beta)
     {
          return 0;
     }
     if (t<mZ2-s/2.*(1+beta))       //mZ2-s/2.*(1+beta)
     {
          return 0;
     }
     if (2*mZ2-s-t>mZ2-s/2.*(1-beta))
     {
          return 0;
     }
     if (2*mZ2-s-t<mZ2-s/2.*(1+beta))
     {
          return 0;
     }
     if (t*(2*mZ2-s-t)<mZ2*mZ2)
     {
          return 0;
     }
     double coeff_W=1;
     double coeff_f=1;
     double coeff_exotic_f=0;
     double coeff_Higgs=1;

     if (exclude_loop==1)//no W boson loop
     {
         coeff_W=0;
     }
     if (exclude_loop==2)//no fermion loop
     {
         coeff_f=0;
     }
     if (exclude_loop==3)//add exotic fermion contrib
     {
         coeff_exotic_f=1;
     }
     if (exclude_loop==4)//no Higgs channel
     {
         coeff_Higgs=0;
     }

     //double beta=pow(1-4*mZ2/(s+0.00001),0.5);
     //printf(" masse %lf\n",mq);
     //double t=mZ2-s/2*(1-pow(1-4*mZ2/(s+0.0001),0.5)*cos(Pi*theta/180.0));
     double u=2*mZ2-s-t;

     double F2_pppp=pow(coeff_Higgs*reFh_pppp(beta)+coeff_W*reFw_pppp(beta,t,mW2)+coeff_f*reFf_pppp_as(beta,t)+coeff_f*reFf_pppp(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_Higgs*imFh_pppp(beta)+coeff_W*imFw_pppp(beta,t,mW2) +coeff_f*imFf_pppp_as(beta,t)+coeff_f*imFf_pppp(beta,t,Itop,Qtop,mtop2),2);

     double F2_pppp_mbeta=pow(coeff_Higgs*reFh_pppp(-beta)+coeff_W*reFw_pppp(-beta,t,mW2)+coeff_f*reFf_pppp_as(-beta,t)+coeff_f*reFf_pppp(-beta,t,Itop,Qtop,mtop2),2)+pow(coeff_Higgs*imFh_pppp(-beta)+coeff_W*imFw_pppp(-beta,t,mW2)+coeff_f*imFf_pppp_as(-beta,t)+coeff_f*imFf_pppp(-beta,t,Itop,Qtop,mtop2),2);
     
     double F2_pppm=pow(coeff_W*reFw_pppm(beta,t,mW2)+coeff_f*reFf_pppm_as(beta,t)+coeff_f*reFf_pppm(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pppm(beta,t,mW2)+coeff_f*imFf_pppm_as(beta,t)+coeff_f*imFf_pppm(beta,t,Itop,Qtop,mtop2),2);
     double F2_pppm_mbeta=pow(coeff_W*reFw_pppm(-beta,t,mW2)+coeff_f*reFf_pppm_as(-beta,t)+coeff_f*reFf_pppm(-beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pppm(-beta,t,mW2)+coeff_f*imFf_pppm_as(-beta,t)+coeff_f*imFf_pppm(-beta,t,Itop,Qtop,mtop2),2);

     double F2_pmpm=pow(coeff_W*reFw_pmpm(beta,t,mW2)+coeff_f*reFf_pmpm_as(beta,t)+coeff_f*reFf_pmpm(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmpm(beta,t,mW2)+coeff_f*imFf_pmpm_as(beta,t)+coeff_f*imFf_pmpm(beta,t,Itop,Qtop,mtop2),2);
     double F2_pmpm_mbeta=pow(coeff_W*reFw_pmpm(-beta,t,mW2)+coeff_f*reFf_pmpm_as(-beta,t)+coeff_f*reFf_pmpm(-beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmpm(-beta,t,mW2)+coeff_f*imFf_pmpm_as(-beta,t)+coeff_f*imFf_pmpm(-beta,t,Itop,Qtop,mtop2),2);

     double F2_pmpp=pow(coeff_W*reFw_pmpp(beta,t,mW2)+coeff_f*reFf_pmpp_as(beta,t)+coeff_f*reFf_pmpp(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmpp(beta,t,mW2)+coeff_f*imFf_pmpp_as(beta,t)+coeff_f*imFf_pmpp(beta,t,Itop,Qtop,mtop2),2);
     double F2_pmpp_mbeta=pow(coeff_W*reFw_pmpp(-beta,t,mW2)+coeff_f*reFf_pmpp_as(-beta,t)+coeff_f*reFf_pmpp(-beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmpp(-beta,t,mW2)+coeff_f*imFf_pmpp_as(-beta,t)+coeff_f*imFf_pmpp(-beta,t,Itop,Qtop,mtop2),2);

     double F2_pp00=pow(coeff_Higgs*reFh_pp00(beta)+coeff_W*reFw_pp00(beta,t,mW2)+coeff_f*reFf_pp00_as(beta,t)+coeff_f*reFf_pp00(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_Higgs*imFh_pp00(beta)+coeff_W*imFw_pp00(beta,t,mW2)+coeff_f*imFf_pp00_as(beta,t)+coeff_f*imFf_pp00(beta,t,Itop,Qtop,mtop2),2);
     double F2_pm00=pow(coeff_W*reFw_pm00(beta,t,mW2)+coeff_f*reFf_pm00_as(beta,t)+coeff_f*reFf_pm00(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pm00(beta,t,mW2)+coeff_f*imFf_pm00_as(beta,t)+coeff_f*imFf_pm00(beta,t,Itop,Qtop,mtop2),2);    

     double F2_ppp0=pow(coeff_W*reFw_ppp0(beta,t,mW2)+coeff_f*reFf_ppp0_as(beta,t)+coeff_f*reFf_ppp0(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_ppp0(beta,t,mW2)+coeff_f*imFf_ppp0_as(beta,t)+coeff_f*imFf_ppp0(beta,t,Itop,Qtop,mtop2),2);
     double F2_ppp0_mbeta=pow(coeff_W*reFw_ppp0(-beta,t,mW2)+coeff_f*reFf_ppp0_as(-beta,t)+coeff_f*reFf_ppp0(-beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_ppp0(-beta,t,mW2)+coeff_f*imFf_ppp0_as(-beta,t)+coeff_f*imFf_ppp0(-beta,t,Itop,Qtop,mtop2),2);
     double F2_ppp0_u=pow(coeff_W*reFw_ppp0(beta,u,mW2)+coeff_f*reFf_ppp0_as(beta,u)+coeff_f*reFf_ppp0(beta,u,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_ppp0(beta,u,mW2)+coeff_f*imFf_ppp0_as(beta,u)+coeff_f*imFf_ppp0(beta,u,Itop,Qtop,mtop2),2);
     double F2_ppp0_mbeta_u=pow(coeff_W*reFw_ppp0(-beta,u,mW2)+coeff_f*reFf_ppp0_as(-beta,u)+coeff_f*reFf_ppp0(-beta,u,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_ppp0(-beta,u,mW2)+coeff_f*imFf_ppp0_as(-beta,u)+coeff_f*imFf_ppp0(-beta,u,Itop,Qtop,mtop2),2);

     double F2_pmp0=pow(coeff_W*reFw_pmp0(beta,t,mW2)+coeff_f*reFf_pmp0_as(beta,t)+coeff_f*reFf_pmp0(beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmp0(beta,t,mW2)+coeff_f*imFf_pmp0_as(beta,t)+coeff_f*imFf_pmp0(beta,t,Itop,Qtop,mtop2),2);
     double F2_pmp0_mbeta=pow(coeff_W*reFw_pmp0(-beta,t,mW2)+coeff_f*reFf_pmp0_as(-beta,t)+coeff_f*reFf_pmp0(-beta,t,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmp0(-beta,t,mW2)+coeff_f*imFf_pmp0_as(-beta,t)+coeff_f*imFf_pmp0(-beta,t,Itop,Qtop,mtop2),2);
     double F2_pmp0_u=pow(coeff_W*reFw_pmp0(beta,u,mW2)+coeff_f*reFf_pmp0_as(beta,u)+coeff_f*reFf_pmp0(beta,u,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmp0(beta,u,mW2)+coeff_f*imFf_pmp0_as(beta,u)+coeff_f*imFf_pmp0(beta,u,Itop,Qtop,mtop2),2);
     double F2_pmp0_mbeta_u=pow(coeff_W*reFw_pmp0(-beta,u,mW2)+coeff_f*reFf_pmp0_as(-beta,u)+coeff_f*reFf_pmp0(-beta,u,Itop,Qtop,mtop2),2)+pow(coeff_W*imFw_pmp0(-beta,u,mW2)+coeff_f*imFf_pmp0_as(-beta,u)+coeff_f*imFf_pmp0(-beta,u,Itop,Qtop,mtop2),2);

     double sum_X=F2_pppp+F2_pppp_mbeta+F2_pppm+F2_pppm_mbeta+F2_pmpp+F2_pmpp_mbeta+F2_pmpm+F2_pmpm_mbeta;
     double sum_Z=F2_pp00+F2_pm00;
     double sum_Y=F2_ppp0+F2_ppp0_mbeta+F2_pmp0+F2_pmp0_mbeta +F2_ppp0_u+F2_ppp0_mbeta_u+F2_pmp0_u+F2_pmp0_mbeta_u;

     double sum=beta/(128.0*s*Pi*100000.1)*(2*sum_X+2*sum_Z+2*sum_Y);
     return sum;
}






















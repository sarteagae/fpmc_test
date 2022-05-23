#include<math.h>
#include<iostream>
#include<fstream>

#include "functions_ZZ.h"


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

constexpr size_t size = 9;
double fermions_I[size] {0};
fermions_I[0] = Iup;
fermions_I[1] = Idown;
fermions_I[2] = Istrange;
fermions_I[3] = Icharm;
fermions_I[4] = Ibottom;
fermions_I[5] = Itop;
fermions_I[6] = Ielectron;
fermions_I[7] = Imuon;
fermions_I[8] = Itau;

double fermions_Q[size] {0};
fermions_Q[0] = Qup;
fermions_Q[1] = Qdown;
fermions_Q[2] = Qstrange;
fermions_Q[3] = Qcharm;
fermions_Q[4] = Qbottom;
fermions_Q[5] = Qtop;
fermions_Q[6] = Qelectron;
fermions_Q[7] = Qmuon;
fermions_Q[8] = Qtau;

double fermions_m2[size] {0};
fermions_m2[0] = mup2;
fermions_m2[1] = mdown2;
fermions_m2[2] = mstrange2;
fermions_m2[3] = mcharm2;
fermions_m2[4] = mbottom2;
fermions_m2[5] = mtop2;
fermions_m2[6] = melectron2;
fermions_m2[7] = mmuon2;
fermions_m2[8] = mtau2;

double reAs_pppm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);//define normalized var for loop functions
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double xt=ut;
     double x1=u1;
     double bracket_u=m*x/Y*ReE1(st,xt,m)+2*m*(1+mZ2*s2/(s4*x1))*ReCz(xt,m)+mZ2*Y/(s4*x1*x1)*(2*x/s-1)*ReBz(xt,m)+2*mZ2*mZ2*m/s4*ReDst(st,xt,m);
     
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=m*x/Y*ReE1(st,xt,m)+2*m*(1+mZ2*s2/(s4*x1))*ReCz(xt,m)+mZ2*Y/(s4*x1*x1)*(2*x/s-1)*ReBz(xt,m)+2*mZ2*mZ2*m/s4*ReDst(st,xt,m);
     
     
     double rea=-4*s2*Y/(t1*u1*s4)+4*s2*m*(s*s4-2*Y)/(s4*Y)*ReC(st,m)+4*s*s4*m/Y*ReCzz(st,m)+8*m*m*ReF(st,tt,ut,m)+4*(mZ2*Y-m*s*s4)/(s*s*s4)*ReE2(tt,ut,m)-8*mZ2*m*Y/(s*s4)*ReDut(tt,ut,m)-4*(bracket_t+bracket_u);
     return rea;
}

double imAs_pppm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double xt=ut;
     double x1=u1;
     double bracket_u=m*x/Y*ImE1(st,xt,m)+2*m*(1+mZ2*s2/(s4*x1))*ImCz(xt,m)+mZ2*Y/(s4*x1*x1)*(2*x/s-1)*ImBz(xt,m)+2*mZ2*mZ2*m/s4*ImDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=m*x/Y*ImE1(st,xt,m)+2*m*(1+mZ2*s2/(s4*x1))*ImCz(xt,m)+mZ2*Y/(s4*x1*x1)*(2*x/s-1)*ImBz(xt,m)+2*mZ2*mZ2*m/s4*ImDst(st,xt,m);
     
     double ima=4*s2*m*(s*s4-2*Y)/(s4*Y)*ImC(st,m)+4*s*s4*m/Y*ImCzz(st,m)+8*m*m*ImF(st,tt,ut,m)+4*(mZ2*Y-m*s*s4)/(s*s*s4)*ImE2(tt,ut,m)-8*mZ2*m*Y/(s*s4)*ImDut(tt,ut,m)-4*(bracket_t+bracket_u);
     return ima;
}

double reAs_pppp(double beta, double t, double u, double m)
{
     //printf("ok ici 4 \n");
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
      
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     //printf("ok ici 5 \n");
     double x=u;
     double xt=ut;
     double x1=u1;
     double bracket_u=2*m*mZ2*mZ2/s4*ReDst(st,xt,m)-(s2+beta*s)*(2*mZ2*Y+x1*(2*x1+s)*(x+mZ2))/(2*s4*s*x1*x1)*ReBz(xt,m)-2*m*(x1*(x-t)+Y)*(s2+beta*s)/(s4*x1*s)*ReCz(xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=2*m*mZ2*mZ2/s4*ReDst(st,xt,m)-(s2+beta*s)*(2*mZ2*Y+x1*(2*x1+s)*(x+mZ2))/(2*s4*s*x1*x1)*ReBz(xt,m)-2*m*(x1*(x-u)+Y)*(s2+beta*s)/(s4*x1*s)*ReCz(xt,m);
     //printf("ok ici 6 %lf \n", bracket_u);
     double rea=4*(mZ2*(2*Y-s*s4)+beta*s*Y)/(s4*t1*u1)+16*mZ2*m/s4*ReC(st,m)+8*m*m*ReF(st,tt,ut,m)+8*Y*m/(s*s4)*(s2+beta*s)*ReDut(tt,ut,m)-2*((s2+beta*s)*Y-4*s*mZ2*m)/(s*s*s4)*ReE2(tt,ut,m)+4*(bracket_u+bracket_t);
     //printf("ok ici 7 %lf \n",rea);
     return rea;
}

double imAs_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double xt=ut;
     double x1=u1;
     double bracket_u=2*m*mZ2*mZ2/s4*ImDst(st,xt,m)-(s2+beta*s)*(2*mZ2*Y+x1*(2*x1+s)*(x+mZ2))/(2*s4*s*x1*x1)*ImBz(xt,m)-2*m*(x1*(x-t)+Y)*(s2+beta*s)/(s4*x1*s)*ImCz(xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=2*m*mZ2*mZ2/s4*ImDst(st,xt,m)-(s2+beta*s)*(2*mZ2*Y+x1*(2*x1+s)*(x+mZ2))/(2*s4*s*x1*x1)*ImBz(xt,m)-2*m*(x1*(x-u)+Y)*(s2+beta*s)/(s4*x1*s)*ImCz(xt,m);
     
     double ima=16*mZ2*m/s4*ImC(st,m)+8*m*m*ImF(st,tt,ut,m)+8*Y*m/(s*s4)*(s2+beta*s)*ImDut(tt,ut,m)-2*((s2+beta*s)*Y-4*s*mZ2*m)/(s*s*s4)*ImE2(tt,ut,m)+4*(bracket_u+bracket_t);
     return ima;
}

double reAs_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double bracket_u=mZ2*(Y+2*x*mZ2)/(s4*x1*x1)*ReBz(xt,m)-2*m*s2*x/(s4*x1)*ReCz(xt,m)+2*m*mZ2*mZ2/s4*ReDst(st,xt,m)-x*(s4*m+mZ2*mZ2)/(s4*Y)*ReE1(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=mZ2*(Y+2*x*mZ2)/(s4*x1*x1)*ReBz(xt,m)-2*m*s2*x/(s4*x1)*ReCz(xt,m)+2*m*mZ2*mZ2/s4*ReDst(st,xt,m)-x*(s4*m+mZ2*mZ2)/(s4*Y)*ReE1(st,xt,m);
     
     double rea=-4*s2*Y/(s4*t1*u1)+4*(s4*m+mZ2*mZ2)/(s4*Y)*(s*s2*ReC(st,m)+(s*s4-2*Y)*ReCzz(st,m))-4*m*s2/(s*s4)*ReE2(tt,ut,m)+8*m*m*ReF(st,tt,ut,m)+4*(bracket_u+bracket_t);
     return rea;
}
double imAs_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double bracket_u=mZ2*(Y+2*x*mZ2)/(s4*x1*x1)*ImBz(xt,m)-2*m*s2*x/(s4*x1)*ImCz(xt,m)+2*m*mZ2*mZ2/s4*ImDst(st,xt,m)-x*(s4*m+mZ2*mZ2)/(s4*Y)*ImE1(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=mZ2*(Y+2*x*mZ2)/(s4*x1*x1)*ImBz(xt,m)-2*m*s2*x/(s4*x1)*ImCz(xt,m)+2*m*mZ2*mZ2/s4*ImDst(st,xt,m)-x*(s4*m+mZ2*mZ2)/(s4*Y)*ImE1(st,xt,m);
     
     double rea=4*(s4*m+mZ2*mZ2)/(s4*Y)*(s*s2*ImC(st,m)+(s*s4-2*Y)*ImCzz(st,m))-4*m*s2/(s*s4)*ImE2(tt,ut,m)+8*m*m*ImF(st,tt,ut,m)+4*(bracket_u+bracket_t);
     return rea;
}

double reAs_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt= ut;
     double bracket_u=2*mZ2/(s4*x1*x1)*(x*x+mZ2*mZ2)*ReBz(xt,m)-8*m*mZ2*Y/(s4*s*x1)*ReCz(xt,m)-s*mZ2*m/s4*ReDst(st,xt,m)+s*x*mZ2/(2*s4*Y)*ReE1(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=2*mZ2/(s4*x1*x1)*(x*x+mZ2*mZ2)*ReBz(xt,m)-8*m*mZ2*Y/(s4*s*x1)*ReCz(xt,m)-s*mZ2*m/s4*ReDst(st,xt,m)+s*x*mZ2/(2*s4*Y)*ReE1(st,xt,m);
     
     double rea=-16*mZ2*Y/(s4*t1*u1)+2*s*s*s2*mZ2/(s4*Y)*ReC(st,m)+2*s*mZ2/(s4*Y)*(s*s4-2*Y)*ReCzz(st,m)-4*(t-u)*(t-u)*mZ2*m/(s4*s)*ReDut(tt,ut,m)-4*(bracket_u+bracket_t);
     return rea;
}
double imAs_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double bracket_u=2*mZ2/(s4*x1*x1)*(x*x+mZ2*mZ2)*ImBz(xt,m)-8*m*mZ2*Y/(s4*s*x1)*ImCz(xt,m)-s*mZ2*m/s4*ImDst(st,xt,m)+s*x*mZ2/(2*s4*Y)*ImE1(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=2*mZ2/(s4*x1*x1)*(x*x+mZ2*mZ2)*ImBz(xt,m)-8*m*mZ2*Y/(s4*s*x1)*ImCz(xt,m)-s*mZ2*m/s4*ImDst(st,xt,m)+s*x*mZ2/(2*s4*Y)*ImE1(st,xt,m);
     
     double ima=2*s*s*s2*mZ2/(s4*Y)*ImC(st,m)+2*s*mZ2/(s4*Y)*(s*s4-2*Y)*ImCzz(st,m)-4*(t-u)*(t-u)*mZ2*m/(s4*s)*ImDut(tt,ut,m)-4*(bracket_u+bracket_t);
     return ima;
}



double reAs_pp00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double bracket_u=2*mZ2/(s4*s*x1*x1)*(2*mZ2*Y+x1*(x-t)*(x+mZ2))*ReBz(xt,m)-8*mZ2*mZ2*m/(s4*x1)*ReCz(xt,m)-s*mZ2*m/s4*ReDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=2*mZ2/(s4*s*x1*x1)*(2*mZ2*Y+x1*(x-u)*(x+mZ2))*ReBz(xt,m)-8*mZ2*mZ2*m/(s4*x1)*ReCz(xt,m)-s*mZ2*m/s4*ReDst(st,xt,m);
     
     double rea=-4*mZ2*m*(t-u)*(t-u)/(s*s4)*ReDut(tt,ut,m)+16*mZ2*Y/(s4*t1*u1)+32*mZ2*m/s4*ReC(st,m)+2*mZ2*(t-u)*(t-u)/(s*s*s4)*ReE2(tt,ut,m)-4*(bracket_u+bracket_t);
     return rea;
}
double imAs_pp00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double xt=ut;
     double x1=u1;
     double bracket_u=2*mZ2/(s4*s*x1*x1)*(2*mZ2*Y+x1*(x-t)*(x+mZ2))*ImBz(xt,m)-8*mZ2*mZ2*m/(s4*x1)*ImCz(xt,m)-s*mZ2*m/s4*ImDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=2*mZ2/(s4*s*x1*x1)*(2*mZ2*Y+x1*(x-u)*(x+mZ2))*ImBz(xt,m)-8*mZ2*mZ2*m/(s4*x1)*ImCz(xt,m)-s*mZ2*m/s4*ImDst(st,xt,m);
     
     double ima=-4*mZ2*m*(t-u)*(t-u)/(s*s4)*ImDut(tt,ut,m)+32*mZ2*m/s4*ImC(st,m)+2*mZ2*(t-u)*(t-u)/(s*s*s4)*ImE2(tt,ut,m)-4*(bracket_u+bracket_t);
     return ima;
}


double reAs_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double bracket_u=(1+beta)*Y/(s*x1*x1)*(s*mZ2-2*x*x1)*ReBz(xt,m)+2*m*(1+beta)*(x*x-mZ2*mZ2+Y)/x1*ReCz(xt,m)+m*(Y+x*x-mZ2*mZ2)*ReDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=(1+beta)*Y/(s*x1*x1)*(s*mZ2-2*x*x1)*ReBz(xt,m)+2*m*(1+beta)*(x*x-mZ2*mZ2+Y)/x1*ReCz(xt,m)+m*(Y+x*x-mZ2*mZ2)*ReDst(st,xt,m);
     
     double rea=delta*(-4*(t-u)/s4*((1+beta)*Y/(t1*u1)+2*m*ReC(st,m)-1/s*(Y*(1+beta)/(2*s)+beta*m)*ReE2(tt,ut,m)+(1+beta)*m*Y/s*ReDut(tt,ut,m))+4/s4*(bracket_t-bracket_u));
     return rea;
}

double imAs_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double bracket_u=(1+beta)*Y/(s*x1*x1)*(s*mZ2-2*x*x1)*ImBz(xt,m)+2*m*(1+beta)*(x*x-mZ2*mZ2+Y)/x1*ImCz(xt,m)+m*(Y+x*x-mZ2*mZ2)*ImDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=(1+beta)*Y/(s*x1*x1)*(s*mZ2-2*x*x1)*ImBz(xt,m)+2*m*(1+beta)*(x*x-mZ2*mZ2+Y)/x1*ImCz(xt,m)+m*(Y+x*x-mZ2*mZ2)*ImDst(st,xt,m);
     
     double ima=delta*(-4*(t-u)/s4*(2*m*ImC(st,m)-1/s*(Y*(1+beta)/(2*s)+beta*m)*ImE2(tt,ut,m)+(1+beta)*m*Y/s*ImDut(tt,ut,m))+4/s4*(bracket_t-bracket_u));
     return ima;
}



double reAs_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=(mZ2*Y-x*x1*(x+mZ2)+xbeta*(mZ2*Y-x*x1*x1))/(s4*x1*x1)*ReBz(xt,m)-((2*mZ2*mZ2+x*s2)*(2*m*Y+s*x*x)+xbeta*s*x*(4*m*Y+s*x*x))/(2*s4*Y*s*x)*ReE1(st,xt,m)+2*m*((2*mZ2*x1+s*x)*Y-xbeta*s*s*x*x)/(s*s4*x*x1)*ReCz(xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=(mZ2*Y-x*x1*(x+mZ2)+xbeta*(mZ2*Y-x*x1*x1))/(s4*x1*x1)*ReBz(xt,m)-((2*mZ2*mZ2+x*s2)*(2*m*Y+s*x*x)+xbeta*s*x*(4*m*Y+s*x*x))/(2*s4*Y*s*x)*ReE1(st,xt,m)+2*m*((2*mZ2*x1+s*x)*Y-xbeta*s*s*x*x)/(s*s4*x*x1)*ReCz(xt,m);
     
     double rea=delta*(-4*(u-t-s*beta)*Y/(s4*t1*u1)+4*(u-t+s*beta)/s4*ReBz(st,m)+2*s/(s4*Y)*((t-u)*(2*mZ2*mZ2-s2*s2)+beta*(4*m*Y+s*(t*t+u*u)))*ReC(st,m) +2*s*s2/(s4*Y)*((u-t)*s4+beta*(s*s4-2*Y))*ReCzz(st,m)+4*m*(t-u)/(s*s4)*ReE2(tt,ut,m)-4*(bracket_t-bracket_u));
     return rea;
}

double imAs_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=(mZ2*Y-x*x1*(x+mZ2)+xbeta*(mZ2*Y-x*x1*x1))/(s4*x1*x1)*ImBz(xt,m)-((2*mZ2*mZ2+x*s2)*(2*m*Y+s*x*x)+xbeta*s*x*(4*m*Y+s*x*x))/(2*s4*Y*s*x)*ImE1(st,xt,m)+2*m*((2*mZ2*x1+s*x)*Y-xbeta*s*s*x*x)/(s*s4*x*x1)*ImCz(xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=(mZ2*Y-x*x1*(x+mZ2)+xbeta*(mZ2*Y-x*x1*x1))/(s4*x1*x1)*ImBz(xt,m)-((2*mZ2*mZ2+x*s2)*(2*m*Y+s*x*x)+xbeta*s*x*(4*m*Y+s*x*x))/(2*s4*Y*s*x)*ImE1(st,xt,m)+2*m*((2*mZ2*x1+s*x)*Y-xbeta*s*s*x*x)/(s*s4*x*x1)*ImCz(xt,m);
     
     double ima=delta*(4*(u-t+s*beta)/s4*ImBz(st,m)+2*s/(s4*Y)*((t-u)*(2*mZ2*mZ2-s2*s2)+beta*(4*m*Y+s*(t*t+u*u)))*ImC(st,m) +2*s*s2/(s4*Y)*((u-t)*s4+beta*(s*s4-2*Y))*ImCzz(st,m)+4*m*(t-u)/(s*s4)*ImE2(tt,ut,m)-4*(bracket_t-bracket_u));
      x=u;
      x1=u1;
      xt=ut;
     xbeta=-beta;
     return ima;
}


double reAs_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=-x/(2*s4*Y*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-xbeta*s*x*(x*x-mZ2*mZ2+Y))*ReE1(st,xt,m)+m/(s4*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-2*xbeta*s*x*(x*x-mZ2*mZ2+Y))*ReDst(st,xt,m)+2*m/(s4*x1)*(mZ2*(s4+xbeta*s)-2*mZ2*Y/s-xbeta*s*x*(x*x-mZ2*mZ2)/Y)*ReCz(xt,m)+((mZ2*mZ2*(x-t)/(s4*x1*x1)-0.5)*(1+xbeta)+mZ2/s4*(1-2*mZ2*mZ2/(x1*x1))+2*x*xbeta/s4-x*x/(Y*s4)*(s4-xbeta*(x-t)))*ReBz(xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=-x/(2*s4*Y*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-xbeta*s*x*(x*x-mZ2*mZ2+Y))*ReE1(st,xt,m)+m/(s4*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-2*xbeta*s*x*(x*x-mZ2*mZ2+Y))*ReDst(st,xt,m)+2*m/(s4*x1)*(mZ2*(s4+xbeta*s)-2*mZ2*Y/s-xbeta*s*x*(x*x-mZ2*mZ2)/Y)*ReCz(xt,m)+((mZ2*mZ2*(x-u)/(s4*x1*x1)-0.5)*(1+xbeta)+mZ2/s4*(1-2*mZ2*mZ2/(x1*x1))+2*x*xbeta/s4-x*x/(Y*s4)*(s4-xbeta*(x-u)))*ReBz(xt,m);
     
     double rea=4*(s2*Y+beta*mZ2*s*(u-t))/(s4*t1*u1)-4*s2*(Y-s*(s4+beta*(u-t)))/(s4*Y)*ReBz(st,m)+4*s*s2/(s4*Y)*((m+s*(s*s4-Y+mZ2*mZ2)/(2*Y))*(s4+beta*(u-t))-s*s4-mZ2*mZ2)*ReC(st,m)+4*s/Y*((m+(t*t*(t*t-mZ2*mZ2+Y)+u*u*(u*u-mZ2*mZ2+Y)+2*Y*(s2*s2-mZ2*mZ2))/(2*Y*s4))*(s4+beta*(u-t))+mZ2*mZ2-s2*s2)*ReCzz(st,m)+8*m*(m-mZ2*Y/(s*s4))*ReDut(tt,ut,m) +4*(bracket_t+bracket_u);
     return rea;
}

double imAs_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=-x/(2*s4*Y*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-xbeta*s*x*(x*x-mZ2*mZ2+Y))*ImE1(st,xt,m)+m/(s4*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-2*xbeta*s*x*(x*x-mZ2*mZ2+Y))*ImDst(st,xt,m)+2*m/(s4*x1)*(mZ2*(s4+xbeta*s)-2*mZ2*Y/s-xbeta*s*x*(x*x-mZ2*mZ2)/Y)*ImCz(xt,m)+((mZ2*mZ2*(x-t)/(s4*x1*x1)-0.5)*(1+xbeta)+mZ2/s4*(1-2*mZ2*mZ2/(x1*x1))+2*x*xbeta/s4-x*x/(Y*s4)*(s4-xbeta*(x-t)))*ImBz(xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=-x/(2*s4*Y*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-xbeta*s*x*(x*x-mZ2*mZ2+Y))*ImE1(st,xt,m)+m/(s4*Y)*(2*m*s4*Y+s*s4*x*x-2*mZ2*mZ2*Y-2*xbeta*s*x*(x*x-mZ2*mZ2+Y))*ImDst(st,xt,m)+2*m/(s4*x1)*(mZ2*(s4+xbeta*s)-2*mZ2*Y/s-xbeta*s*x*(x*x-mZ2*mZ2)/Y)*ImCz(xt,m)+((mZ2*mZ2*(x-u)/(s4*x1*x1)-0.5)*(1+xbeta)+mZ2/s4*(1-2*mZ2*mZ2/(x1*x1))+2*x*xbeta/s4-x*x/(Y*s4)*(s4-xbeta*(x-u)))*ImBz(xt,m);
     
     double ima=-4*s2*(Y-s*(s4+beta*(u-t)))/(s4*Y)*ImBz(st,m)+4*s*s2/(s4*Y)*((m+s*(s*s4-Y+mZ2*mZ2)/(2*Y))*(s4+beta*(u-t))-s*s4-mZ2*mZ2)*ImC(st,m)+4*s/Y*((m+(t*t*(t*t-mZ2*mZ2+Y)+u*u*(u*u-mZ2*mZ2+Y)+2*Y*(s2*s2-mZ2*mZ2))/(2*Y*s4))*(s4+beta*(u-t))+mZ2*mZ2-s2*s2)*ImCzz(st,m)+8*m*(m-mZ2*Y/(s*s4))*ImDut(tt,ut,m) +4*(bracket_t+bracket_u);
     return ima;
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
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=8*mZ2*s*beta/s4*ReC(st,m)+4*(2*cW*cW*s4*(s2+beta*s)+mZ2*(s4+beta*s))/(s*s4)*ReE2(tt,ut,m)-4*cW*cW*(2*mZ2*mZ2+(4*mW2-s)*(s2+beta*s))*ReF(st,tt,ut,m)+2*s*mZ2/s4*(s4+beta*s2)*(ReDst(st,tt,m)+ReDst(st,ut,m));
     return red;
}
double im_dw_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     //printf("%lf %lf %lf %lf",ut,st,tt,m);
     double imd=0*8*mZ2*s*beta/s4*ImC(st,m)+0*4*(2*cW*cW*s4*(s2+beta*s)+mZ2*(s4+beta*s))/(s*s4)*ImE2(tt,ut,m)-4*cW*cW*(2*mZ2*mZ2+(4*mW2-s)*(s2+beta*s))*ImF(st,tt,ut,m)+0*2*s*mZ2/s4*(s4+beta*s2)*(ImDst(st,tt,m)+ImDst(st,ut,m));
     //printf("\n imd d %lf \n",-4*cW*cW*(2*mZ2*mZ2+(4*mW2-s)*(s2+beta*s))* ImF(st,tt,ut,m));
     return imd;
}

double re_dw_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     
     double bracket_u=4*(4*mW2*(x+mZ2)*(x+mZ2)+x*s4*mZ2)/(s4*Y)*ReE1(st,xt,m)+4*mZ2*Y/s4*ReDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     
     double bracket_t=4*(4*mW2*(x+mZ2)*(x+mZ2)+x*s4*mZ2)/(s4*Y)*ReE1(st,xt,m)+4*mZ2*Y/s4*ReDst(st,xt,m);
     
     double red= 4*s/(s4*Y)*(s2*s4*(4*mW2-mZ2)-8*Y*mW2)*ReC(st,m)+4*(4*mW2-mZ2)*(s*s4-2*Y)/Y*ReCzz(st,m)+8*mW2*(4*mW2-mZ2+2*Y/s4)*ReF(st,tt,ut,m)+(bracket_t+bracket_u);
     return red;
}

double im_dw_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     
     double bracket_u=4*(4*mW2*(x+mZ2)*(x+mZ2)+x*s4*mZ2)/(s4*Y)*ImE1(st,xt,m)+4*mZ2*Y/s4*ImDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     
     double bracket_t=4*(4*mW2*(x+mZ2)*(x+mZ2)+x*s4*mZ2)/(s4*Y)*ImE1(st,xt,m)+4*mZ2*Y/s4*ImDst(st,xt,m);
     
     double imd= 4*s/(s4*Y)*(s2*s4*(4*mW2-mZ2)-8*Y*mW2)*ImC(st,m)+4*(4*mW2-mZ2)*(s*s4-2*Y)/Y*ImCzz(st,m)+8*mW2*(4*mW2-mZ2+2*Y/s4)*ImF(st,tt,ut,m)+(bracket_t+bracket_u);
     return imd;
}

double re_dw_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     
     double bracket_u=4*x1/(s4*Y)*(8*mW2*(Y+(x+mZ2)*(x+mZ2))-s2*s4*x)*ReCz(xt,m)+2/(s4*Y)*(8*mW2*(x*x*s*s4-2*Y*x1*x1)+s*((x*x-mZ2*mZ2)*(x*x-mZ2*mZ2)-2*mZ2*x*x*s4-s*x*Y))*ReDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     
     double bracket_t=4*x1/(s4*Y)*(8*mW2*(Y+(x+mZ2)*(x+mZ2))-s2*s4*x)*ReCz(xt,m)+2/(s4*Y)*(8*mW2*(x*x*s*s4-2*Y*x1*x1)+s*((x*x-mZ2*mZ2)*(x*x-mZ2*mZ2)-2*mZ2*x*x*s4-s*x*Y))*ReDst(st,xt,m);
     
     double red= 2*s/(s4*Y)*(8*mW2*(s2*s4-4*Y)+s4*(s2*s2-2*Y))*ReC(st,m)+2/Y*(s2+8*mW2)*(s*s4-2*Y)*ReCzz(st,m)+4*mW2*(s2+8*mW2)*ReF(st,tt,ut,m)-2/s4*(s4+16*mW2)*ReE2(tt,ut,m)+(bracket_t+bracket_u);
     return red;
}

double im_dw_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double bracket_u=4*x1/(s4*Y)*(8*mW2*(Y+(x+mZ2)*(x+mZ2))-s2*s4*x)*ImCz(xt,m)+2/(s4*Y)*(8*mW2*(x*x*s*s4-2*Y*x1*x1)+s*((x*x-mZ2*mZ2)*(x*x-mZ2*mZ2)-2*mZ2*x*x*s4-s*x*Y))*ImDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     double bracket_t=4*x1/(s4*Y)*(8*mW2*(Y+(x+mZ2)*(x+mZ2))-s2*s4*x)*ImCz(xt,m)+2/(s4*Y)*(8*mW2*(x*x*s*s4-2*Y*x1*x1)+s*((x*x-mZ2*mZ2)*(x*x-mZ2*mZ2)-2*mZ2*x*x*s4-s*x*Y))*ImDst(st,xt,m);
     
     double imd= 2*s/(s4*Y)*(8*mW2*(s2*s4-4*Y)+s4*(s2*s2-2*Y))*ImC(st,m)+2/Y*(s2+8*mW2)*(s*s4-2*Y)*ImCzz(st,m)+4*mW2*(s2+8*mW2)*ImF(st,tt,ut,m)-2/s4*(s4+16*mW2)*ImE2(tt,ut,m)+(bracket_t+bracket_u);
     //printf(" test %lf\n",pow(10,-13)*2/(t)*(((u*u*u))));//8*mW2*(u*u*s*s4-2*Y*u1*u1)+;;*ImDst(st,ut,m)
     //printf(" test %lf\n",imd);
     return imd;
}

double re_dw_pp00(double beta, double t, double u, double m)
{    
     double s=4*mZ2/(1.0-beta*beta);
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red= 4*s*ReC(st,m)-4*mW2*(s+2*mZ2-8*mW2)*ReF(st,tt,ut,m)-4*(4*mW2-mZ2)/s*ReE2(tt,ut,m);
     //printf(" %lf\n",red);
     //printf(" test %lf\n",4*s*ReC(st,m));
     return red;
}

double im_dw_pp00(double beta, double t, double u, double m)
{    
     double s=4*mZ2/(1-beta*beta);
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd= 4*s*ImC(st,m)-4*mW2*(s+2*mZ2-8*mW2)*ImF(st,tt,ut,m)-4*(4*mW2-mZ2)/s*ImE2(tt,ut,m);
     return imd;
}

double re_dw_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=delta*(-(t-u)*(s4+beta*s)/(s*s4)*(2*s*ReC(st,m)+ReE2(tt,ut,m))+(s4+beta*s)/s4*((t*t-mZ2*mZ2+Y)*ReDst(st,tt,m)-(u*u-mZ2*mZ2+Y)*ReDst(st,ut,m)));
     //printf(" %lf\n",red);
     return red;
}

double im_dw_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd= delta*(-(t-u)*(s4+beta*s)/(s*s4)*(2*s*ImC(st,m)+ImE2(tt,ut,m))+(s4+beta*s)/s4*((t*t-mZ2*mZ2+Y)*ImDst(st,tt,m)-(u*u-mZ2*mZ2+Y)*ImDst(st,ut,m)));;
     return imd;
}

double re_dw_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=(t-x-xbeta*s)/s4*(2*mZ2*mZ2+2*x*x+s*x-8*cW*cW*Y)*ReDst(st,xt,m)-16*cW*cW/s4*(x+mZ2)*(1+xbeta)*ReE1(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=(u-x-xbeta*s)/s4*(2*mZ2*mZ2+2*x*x+s*x-8*cW*cW*Y)*ReDst(st,xt,m)-16*cW*cW/s4*(x+mZ2)*(1+xbeta)*ReE1(st,xt,m);
     
     double red= delta*(-2/s4*((u-t-beta*s)*s4+8*cW*cW*s*(u-t+beta*s4))*ReC(st,m)+16*cW*cW*(t-u-beta*s)*ReCzz(st,m)-64*cW*cW*mW2*beta*Y/s4*ReF(st,tt,ut,m)+8*cW*cW*Y/s4*(u-t-beta*s)*ReDut(tt,ut,m)+(t-u+beta*s)/s*ReE2(tt,ut,m)-(bracket_t-bracket_u));
     return red;
}

double im_dw_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2.0*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=(t-x-xbeta*s)/s4*(2*mZ2*mZ2+2*x*x+s*x-8*cW*Y)*ImDst(st,xt,m)-16*cW*cW/s4*(x+mZ2)*(1+xbeta)*ImE1(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=(u-x-xbeta*s)/s4*(2*mZ2*mZ2+2*x*x+s*x-8*cW*cW*Y)*ImDst(st,xt,m)-16*cW*cW/s4*(x+mZ2)*(1+xbeta)*ImE1(st,xt,m);
     //printf(" %lf\n",delta);
     double imd= delta*(-2/s4*((u-t-beta*s)*s4+8*cW*cW*s*(u-t+beta*s4))*ImC(st,m)+16*cW*cW*(t-u-beta*s)*ImCzz(st,m)-64*cW*cW*mW2*beta*Y/s4*ImF(st,tt,ut,m)+8*cW*cW*Y/s4*(u-t-beta*s)*ImDut(tt,ut,m)+(t-u+beta*s)/s*ImE2(tt,ut,m)-(bracket_t-bracket_u));
     return imd;
}

double re_dw_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=8*cW*cW*(s2+xbeta*s)/s4*ReE1(st,xt,m)+2*mZ2/s4*(s*(s4+xbeta*(x-t))-2*Y)*ReDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=8*cW*cW*(s2+xbeta*s)/s4*ReE1(st,xt,m)+2*mZ2/s4*(s*(s4+xbeta*(x-u))-2*Y)*ReDst(st,xt,m);
     
     double red= 16*cW*cW*s*(s2/s4*ReC(st,m)+ReCzz(st,m))+4*cW*cW/s4*(s*(s2+4*mW2)*(s4+beta*(t-u))-2*Y*s2)*ReF(st,tt,ut,m)+(bracket_t+bracket_u);
     return red;
}

double im_dw_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double x=u;
     double x1=u1;
     double xt=ut;
     double xbeta=-beta;
     double bracket_u=8*cW*cW*(s2+xbeta*s)/s4*ImE1(st,xt,m)+2*mZ2/s4*(s*(s4+xbeta*(x-t))-2*Y)*ImDst(st,xt,m);
     x=t;
     x1=t1;
     xt=tt;
     xbeta=beta;
     double bracket_t=8*cW*cW*(s2+xbeta*s)/s4*ImE1(st,xt,m)+2*mZ2/s4*(s*(s4+xbeta*(x-u))-2*Y)*ImDst(st,xt,m);
     
     double imd= 16*cW*cW*s*(s2/s4*ImC(st,m)+ImCzz(st,m))+4*cW*cW/s4*(s*(s2+4*mW2)*(s4+beta*(t-u))-2*Y*s2)*ImF(st,tt,ut,m)+(bracket_t+bracket_u);
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
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     //printf("ok ici 2 %lf %lf %lf \n",beta,delta,st);
     
     double red= 4*(s2+s*beta)*(m*ReF(st,tt,ut,m)-1/(2*s)*ReE2(tt,ut,m));
     //printf("ok ici 3 %lf \n",red);
     return red;
}

double im_dvf_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd= 4*(s2+s*beta)*(m*ImF(st,tt,ut,m)-1/(2.*s)*ImE2(tt,ut,m));
     return imd;
}

double re_daf_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red= -8*s*m/s4*(4*beta*ReC(st,m)+(s4+beta*s2)*(ReDst(st,tt,m)+ReDst(st,ut,m)))+4*m*(s2+s*beta+8*m)*ReF(st,tt,ut,m)-2/s*(s2+s*beta+8*m*(s4+s*beta)/s4)*ReE2(tt,ut,m);
     return red;
}

double im_daf_pppp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd= -8*s*m/s4*(4*beta*ImC(st,m)+(s4+beta*s2)*(ImDst(st,tt,m)+ImDst(st,ut,m)))+4*m*(s2+s*beta+8*m)*ImF(st,tt,ut,m)-2/s*(s2+s*beta+8*m*(s4+s*beta)/s4)*ImE2(tt,ut,m);
     return imd;
}

double re_dvf_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red= -4*mZ2/Y*(s*(s2*s4-2*Y)/s4*ReC(st,m)+(s*s4-2*Y)*ReCzz(st,m))-4*mZ2/(s4*Y)*((t+mZ2)*(t+mZ2)*ReE1(st,tt,m)+(u+mZ2)*(u+mZ2)*ReE1(st,ut,m))-8*mZ2*m*ReF(st,tt,ut,m);
     return red;
}

double im_dvf_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd= -4*mZ2/Y*(s*(s2*s4-2*Y)/s4*ImC(st,m)+(s*s4-2*Y)*ImCzz(st,m))-4*mZ2/(s4*Y)*((t+mZ2)*(t+mZ2)*ImE1(st,tt,m)+(u+mZ2)*(u+mZ2)*ImE1(st,ut,m))-8*mZ2*m*ImF(st,tt,ut,m);
     return imd;
}

double re_daf_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=8*s*(s2*(4*m-mZ2)/(2*Y)+mZ2/s4)*ReC(st,m)+4*(s*s4-2*Y)/Y*(4*m-mZ2)*ReCzz(st,m)+8*m*(4*m-mZ2)*ReF(st,tt,ut,m)+16*m*Y/s4*ReDut(tt,ut,m)-4/Y*((mZ2*(t+mZ2)*(t+mZ2)/s4+4*m*t)*ReE1(st,tt,m)+(mZ2*(u+mZ2)*(u+mZ2)/s4+4*m*u)*ReE1(st,ut,m));
     return red;
}

double im_daf_pmpp(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd= 8*s*(s2*(4*m-mZ2)/(2*Y)+mZ2/s4)*ImC(st,m)+4*(s*s4-2*Y)/Y*(4*m-mZ2)*ImCzz(st,m)+8*m*(4*m-mZ2)*ImF(st,tt,ut,m)+16*m*Y/s4*ImDut(tt,ut,m)-4/Y*((mZ2*(t+mZ2)*(t+mZ2)/s4+4*m*t)*ImE1(st,tt,m)+(mZ2*(u+mZ2)*(u+mZ2)/s4+4*m*u)*ImE1(st,ut,m));
     return imd;
}

double re_dvf_pp00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red= -8*m*mZ2*ReF(st,tt,ut,m)+4*mZ2/s*ReE2(tt,ut,m);
     return red;
}

double im_dvf_pp00(double beta, double t, double u, double m)
{  
     double s=4*mZ2/(1-beta*beta);
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd= -8*m*mZ2*ImF(st,tt,ut,m)+4*mZ2/s*ImE2(tt,ut,m);
     return imd;
}

double re_daf_pp00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=-8*m*(mZ2+2*s2*m/mZ2)*ReF(st,tt,ut,m)-4/s*(4*m-mZ2)*ReE2(tt,ut,m)-16*m*s/mZ2*ReC(st,m);
     return red;
}

double im_daf_pp00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=-8*m*(mZ2+2*s2*m/mZ2)*ImF(st,tt,ut,m)-4/s*(4*m-mZ2)*ImE2(tt,ut,m)-16*m*s/mZ2*ImC(st,m) ;
     return imd;
}

double re_dvf_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=-4*s*mZ2/(s4*Y)*(s2*s4-4*Y)*ReC(st,m)-4*mZ2/Y*(s*s4-2*Y)*ReCzz(st,m)-8*mZ2*m*ReF(st,tt,ut,m)-4*mZ2/(s4*Y)*((2*t*t+s*t+2*mZ2*mZ2)*ReE1(st,tt,m)+(2*u*u+s*u+2*mZ2*mZ2)*ReE1(st,ut,m));
     return red;
}

double im_dvf_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=-4*s*mZ2/(s4*Y)*(s2*s4-4*Y)*ImC(st,m)-4*mZ2/Y*(s*s4-2*Y)*ImCzz(st,m)-8*mZ2*m*ImF(st,tt,ut,m)-4*mZ2/(s4*Y)*((2*t*t+s*t+2*mZ2*mZ2)*ImE1(st,tt,m)+(2*u*u+s*u+2*mZ2*mZ2)*ImE1(st,ut,m));
     return imd;
}

double re_daf_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=-8*m/mZ2*(s*(s2*s2/Y-2)*ReC(st,m)+s2*(s*s4/Y-2)*ReCzz(st,m)-4*mZ2*Y/s4*ReDut(tt,ut,m)+(t*t+mZ2*mZ2)/Y*ReE1(st,tt,m)+(u*u+mZ2*mZ2)/Y*ReE1(st,ut,m)+(2*s2*m+mZ2*mZ2)*ReF(st,tt,ut,m))-4*mZ2/Y*(s/s4*(s2*s4-4*Y)*ReC(st,m)+(s*s4-2*Y)*ReCzz(st,m)+(2*t*t+s*t+2*mZ2*mZ2)/s4*ReE1(st,tt,m)+(2*u*u+s*u+2*mZ2*mZ2)/s4*ReE1(st,ut,m));
     return red;
}

double im_daf_pm00(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=-8*m/mZ2*(s*(s2*s2/Y-2)*ImC(st,m)+s2*(s*s4/Y-2)*ImCzz(st,m)-4*mZ2*Y/s4*ImDut(tt,ut,m)+(t*t+mZ2*mZ2)/Y*ImE1(st,tt,m)+(u*u+mZ2*mZ2)/Y*ImE1(st,ut,m)+(2*s2*m+mZ2*mZ2)*ImF(st,tt,ut,m))-4*mZ2/Y*(s/s4*(s2*s4-4*Y)*ImC(st,m)+(s*s4-2*Y)*ImCzz(st,m)+(2*t*t+s*t+2*mZ2*mZ2)/s4*ImE1(st,tt,m)+(2*u*u+s*u+2*mZ2*mZ2)/s4*ImE1(st,ut,m));
     return imd;
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
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=delta*(-4*(s4+beta*s)*m/(s4*mZ2)*((u-t)*(2*ReC(st,m)+1/s*ReE2(tt,ut,m))-(s*t-2*mZ2*t1)*ReDst(st,tt,m)+(s*u-2*mZ2*u1)*ReDst(st,ut,m)));
     return red;
}

double im_daf_ppp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=delta*(-4*(s4+beta*s)*m/(s4*mZ2)*((u-t)*(2*ImC(st,m)+1/s*ImE2(tt,ut,m))-(s*t-2*mZ2*t1)*ImDst(st,tt,m)+(s*u-2*mZ2*u1)*ImDst(st,ut,m)));
     return imd;
}

double re_dvf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=delta*(-4*s/s4*(t-u-beta*s4)*ReC(st,m)-4*(t-u-beta*s)*ReCzz(st,m)+16*beta*m*Y/s4*ReF(st,tt,ut,m)-4*(1+beta)*(t+mZ2)/s4*ReE1(st,tt,m)+4*(1-beta)*(u+mZ2)/s4*ReE1(st,ut,m));
     return red;
}

double im_dvf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=delta*(-4*s/s4*(t-u-beta*s4)*ImC(st,m)-4*(t-u-beta*s)*ImCzz(st,m)+16*beta*m*Y/s4*ImF(st,tt,ut,m)-4*(1+beta)*(t+mZ2)/s4*ImE1(st,tt,m)+4*(1-beta)*(u+mZ2)/s4*ImE1(st,ut,m));
     return imd;
}
double re_daf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=delta*(-4*m/(mZ2*s)*(t-u+beta*s)*(2*s*ReC(st,m)+ReE1(st,tt,m)+ReE1(st,ut,m)+(s+4*mZ2)*Y/s4*ReDut(tt,ut,m))+16*beta*m*Y/s4*ReF(st,tt,ut,m)-4*s/s4*(t-u-beta*s4)*ReC(st,m)-4*(t-u-beta*s)*ReCzz(st,m)-4*(1+beta)*(t+mZ2)/s4*ReE1(st,tt,m)+4*(1-beta)*(u+mZ2)/s4*ReE1(st,ut,m));
     return red;
}

double im_daf_pmp0(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=delta*(-4*m/(mZ2*s)*(t-u+beta*s)*(2*s*ImC(st,m)+ImE1(st,tt,m)+ImE1(st,ut,m)+(s+4*mZ2)*Y/s4*ImDut(tt,ut,m))+16*beta*m*Y/s4*ImF(st,tt,ut,m)-4*s/s4*(t-u-beta*s4)*ImC(st,m)-4*(t-u-beta*s)*ImCzz(st,m)-4*(1+beta)*(t+mZ2)/s4*ImE1(st,tt,m)+4*(1-beta)*(u+mZ2)/s4*ImE1(st,ut,m));
     return imd;
}

double re_dvf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=-4*s*(s2/s4*ReC(st,m)+ReCzz(st,m))-2/s4*((s2+beta*s)*ReE1(st,tt,m)+(s2-beta*s)*ReE1(st,ut,m))-4*s*m/s4*(s4+beta*(t-u))*ReF(st,tt,ut,m);
     return red;
}

double im_dvf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=-4*s*(s2/s4*ImC(st,m)+ImCzz(st,m))-2/s4*((s2+beta*s)*ImE1(st,tt,m)+(s2-beta*s)*ImE1(st,ut,m))-4*s*m/s4*(s4+beta*(t-u))*ImF(st,tt,ut,m);
     return imd;
}
double re_daf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double red=-4*s*(s2/s4*ReC(st,m)+ReCzz(st,m))-2/s4*((s2+beta*s)*ReE1(st,tt,m)+(s2-beta*s)*ReE1(st,ut,m))-4*s*m/s4*(s4+beta*(t-u))*(ReF(st,tt,ut,m)-2*ReDut(tt,ut,m))-16*m*Y/s4*ReDut(tt,ut,m);
     return red;
}

double im_daf_pmpm(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double s2=s-2*mZ2;
     double s4=s-4*mZ2;
     double u1=u-mZ2;
     double t1=t-mZ2;
     double Y=u*t-mZ2*mZ2;
     double delta=pow(s*mZ2/(2*Y),0.5);
     
     double st=s/(4.0*m);
     double tt=t/(4.0*m);
     double ut=u/(4.0*m);
     
     double imd=-4*s*(s2/s4*ImC(st,m)+ImCzz(st,m))-2/s4*((s2+beta*s)*ImE1(st,tt,m)+(s2-beta*s)*ImE1(st,ut,m))-4*s*m/s4*(s4+beta*(t-u))*(ImF(st,tt,ut,m)-2*ImDut(tt,ut,m))-16*m*Y/s4*ImDut(tt,ut,m);;
     return imd;
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


//asymptot w


double reAs_pppm_as(double beta, double t, double u, double m)
{
     return -4;
}

double imAs_pppm_as(double beta, double t, double u, double m)
{
     return 0;
}

double reAs_pppp_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return 4-4*u*t/(s*s)*(pow(log(t/u),2)+Pi*Pi)+4*(t-u)/s*log(t/u);
}

double imAs_pppp_as(double beta, double t, double u, double m)
{
     return 0;
}

double reAs_pmpp_as(double beta, double t, double u, double m)
{
     return -4;
}
double imAs_pmpp_as(double beta, double t, double u, double m)
{
     return 0;
}

double reAs_pm00_as(double beta, double t, double u, double m)
{
     return 0;
}
double imAs_pm00_as(double beta, double t, double u, double m)
{
     return 0;
}



double reAs_pp00_as(double beta, double t, double u, double m)
{
     return 0;
}
double imAs_pp00_as(double beta, double t, double u, double m)
{
     return 0;
}


double reAs_ppp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2*u*t),0.5);
     
     return delta*(-8*(t-u)/s+4*(t-u)*t*u/pow(s,3)*(pow(log(t/u),2)+Pi*Pi)+16*u*t/pow(s,2)*log(t/u));
}

double imAs_ppp0_as(double beta, double t, double u, double m)
{
     return 0;
}



double reAs_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2*u*t),0.5);
     
     return delta*(-8*u/s+4*t/u*pow(log(-s/t),2)+8*t/s*log(-s/t));
}

double imAs_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2*u*t),0.5);
     
     return delta*(4*t/u*(-2*Pi*log(-s/t))-8*t*Pi/s);
}


double reAs_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
 
     return 4-4*s*t/(u*u)*(pow(log(-s/t),2))+4*(s-t)/u*log(-s/t);
}

double imAs_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -4*s*t/(u*u)*(-2*Pi*log(-s/t))+4*(s-t)/u*(-Pi);
}

//W delat 
double re_dw_pppm_as(double beta, double t, double u, double m)
{
     return 0;
}
double im_dw_pppm_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_dw_pppp_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double red=16*cW*cW*(pow(log(t/u),2)+Pi*Pi+s/u*log(-u/m)*log(-s/t)+s/t*log(-t/m)*log(-s/u));
     return red;
}
double im_dw_pppp_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double imd=16*cW*cW*(-Pi*(s/u*log(-u/m)+s/t*log(-t/m)));
     return imd;
}

double re_dw_pmpp_as(double beta, double t, double u, double m)
{
     
     return 0;
}

double im_dw_pmpp_as(double beta, double t, double u, double m)
{
     
     return 0;
}

double re_dw_pm00_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double red= -2*(pow(log(t/u),2)+s/u*pow(log(-t/s),2)+s/t*pow(log(-u/s),2))-4/s*log(s/m)*(t*log(-t/s)+u*log(-u/s));
     return red;
}

double im_dw_pm00_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double imd= 4*Pi*(log(s/m)+(1-t*t/(s*u))*log(-t/s)+(1-u*u/(s*t))*log(-u/s));
     return imd;
}

double re_dw_pp00_as(double beta, double t, double u, double m)
{    
     double s=4*mZ2/(1.0-beta*beta);
     
     return 2*pow(log(s/m),2)-2*Pi*Pi;
}

double im_dw_pp00_as(double beta, double t, double u, double m)
{    
     double s=4*mZ2/(1-beta*beta);
     
     return  -4*Pi*log(s/m);
}

double re_dw_ppp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     
     double red=delta*(-2*(t-u)/s*(pow(log(s/m),2)+pow(log(t/u),2))-4*log(s/m)*log(t/u));
     //printf(" %lf\n",red);
     return red;
}

double im_dw_ppp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     double imd= delta*(-2*(t-u)/s*(-2*Pi*log(s/m))+4*Pi*log(t/u));;
     return imd;
}

double re_dw_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     double red= delta*((-2*u/s+32*cW*cW*t/s)*(pow(log(s/m),2))+32*cW*cW*(u/s*log(-t/m)*log(-u/m)+t/s*pow(log(-t/m),2))-2*u/s*pow(log(t/u),2)+(-4*u*(t-u)/(s*s)+32*cW*cW*(t*t+s*s)/(s*s))*(log(s/m))*log(-t/m)-4*u/(s*s)*(u-t-8*cW*cW*t)*log(s/m)*log(-u/m));
     return red;
}

double im_dw_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1.0-beta*beta);
     double delta=pow(s*mZ2/(2.0*u*t),0.5);
     
     double imd= delta*((-2*u/s+32*cW*cW*t/s)*(-2*Pi*log(s/m))+(-4*u*(t-u)/(s*s)+32*cW*cW*(t*t+s*s)/(s*s))*(-Pi)*log(-t/m)-4*u/(s*s)*(u-t-8*cW*cW*t)*(-Pi)*log(-u/m));
     return imd;
}

double re_dw_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     double red= 16*cW*cW*(pow(log(-t/s),2)+u/t*log(-t/m)*log(-u/s)+u/s*log(s/m)*log(u/t));
     return red;
}

double im_dw_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double imd= 16*cW*cW*Pi*(u/t*log(s/m)-u/s*log(-u/s)-(s*s+t*t)/(s*t)*log(-t/s));
     return imd;
}

//W amplitudes

double reAw_pppm_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pppm_as(beta,t,u,m)+re_dw_pppm_as(beta,t,u,m);
}

double imAw_pppm_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pppm_as(beta,t,u,m)+im_dw_pppm_as(beta,t,u,m);
}

double reAw_pppp_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pppp_as(beta,t,u,m)+re_dw_pppp_as(beta,t,u,m);
}

double imAw_pppp_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pppp_as(beta,t,u,m)+im_dw_pppp_as(beta,t,u,m);
}

double reAw_pmpp_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pmpp_as(beta,t,u,m)+re_dw_pmpp_as(beta,t,u,m);
}

double imAw_pmpp_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pmpp_as(beta,t,u,m)+im_dw_pmpp_as(beta,t,u,m);
}

double reAw_pm00_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pm00_as(beta,t,u,m)+re_dw_pm00_as(beta,t,u,m);
}

double imAw_pm00_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pm00_as(beta,t,u,m)+im_dw_pm00_as(beta,t,u,m);
}

double reAw_pp00_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pp00_as(beta,t,u,m)+re_dw_pp00_as(beta,t,u,m);
}

double imAw_pp00_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pp00_as(beta,t,u,m)+im_dw_pp00_as(beta,t,u,m);
}

double reAw_ppp0_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_ppp0_as(beta,t,u,m)+re_dw_ppp0_as(beta,t,u,m);
}

double imAw_ppp0_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_ppp0_as(beta,t,u,m)+im_dw_ppp0_as(beta,t,u,m);
}

double reAw_pmp0_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pmp0_as(beta,t,u,m)+re_dw_pmp0_as(beta,t,u,m);
}

double imAw_pmp0_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pmp0_as(beta,t,u,m)+im_dw_pmp0_as(beta,t,u,m);
}

double reAw_pmpm_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*reAs_pmpm_as(beta,t,u,m)+re_dw_pmpm_as(beta,t,u,m);
}

double imAw_pmpm_as(double beta, double t, double u, double m)
{
     return (12*cW*cW*cW*cW-4*cW*cW+1)/(4*cW*cW)*imAs_pmpm_as(beta,t,u,m)+im_dw_pmpm_as(beta,t,u,m);
}


double reFw_pppm_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pppm_as(beta,t,u,m);
}

double imFw_pppm_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pppm_as(beta,t,u,m);
}

double reFw_pppp_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pppp_as(beta,t,u,m);
}

double imFw_pppp_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pppp_as(beta,t,u,m);
}

double reFw_pmpp_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pmpp_as(beta,t,u,m);
}

double imFw_pmpp_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pmpp_as(beta,t,u,m);
}

double reFw_pm00_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pm00_as(beta,t,u,m);
}

double imFw_pm00_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pm00_as(beta,t,u,m);
}

double reFw_pp00_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     
     return alpha*alpha/(sW*sW)*reAw_pp00_as(beta,t,u,m);
}

double imFw_pp00_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pp00_as(beta,t,u,m);
}

double reFw_ppp0_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_ppp0_as(beta,t,u,m);
}

double imFw_ppp0_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_ppp0_as(beta,t,u,m);
}

double reFw_pmp0_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pmp0_as(beta,t,u,m);
}

double imFw_pmp0_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pmp0_as(beta,t,u,m);
}

double reFw_pmpm_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*reAw_pmpm_as(beta,t,u,m);
}

double imFw_pmpm_as(double beta, double t, double m)
{
     double s=4*mZ2/(1-beta*beta);
     double u=2*mZ2-s-t;
     return alpha*alpha/(sW*sW)*imAw_pmpm_as(beta,t,u,m);
}






//
//fermion
//

double reFf_pppm(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pppm(beta,t,u,m)+gv*gv*reAvf_pppm(beta,t,u,m));
     }
     return sum;
}

double imFf_pppm(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pppm(beta,t,u,m)+gv*gv*imAvf_pppm(beta,t,u,m));
     }
     return sum;
}

double reFf_pppp(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pppp(beta,t,u,m)+gv*gv*reAvf_pppp(beta,t,u,m));
     }
     return sum;
}


double imFf_pppp(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pppp(beta,t,u,m)+gv*gv*imAvf_pppp(beta,t,u,m));
     }
     return sum;
}

double reFf_pmpp(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pmpp(beta,t,u,m)+gv*gv*reAvf_pmpp(beta,t,u,m));
     }
     return sum;
}

double imFf_pmpp(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pmpp(beta,t,u,m)+gv*gv*imAvf_pmpp(beta,t,u,m));
     }
     return sum;
}

double reFf_pm00(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pm00(beta,t,u,m)+gv*gv*reAvf_pm00(beta,t,u,m));
     }
     return sum;
}

double imFf_pm00(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pm00(beta,t,u,m)+gv*gv*imAvf_pm00(beta,t,u,m));
     }
     return sum;
}

double reFf_pp00(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pp00(beta,t,u,m)+gv*gv*reAvf_pp00(beta,t,u,m));
     }
     return sum;
}

double imFf_pp00(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pp00(beta,t,u,m)+gv*gv*imAvf_pp00(beta,t,u,m));
     }
     return sum;
}

double reFf_ppp0(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_ppp0(beta,t,u,m)+gv*gv*reAvf_ppp0(beta,t,u,m));
     }
     return sum;
}

double imFf_ppp0(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_ppp0(beta,t,u,m)+gv*gv*imAvf_ppp0(beta,t,u,m));
     }
     return sum;
}

double reFf_pmp0(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pmp0(beta,t,u,m)+gv*gv*reAvf_pmp0(beta,t,u,m));
     }
     return sum;
}

double imFf_pmp0(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pmp0(beta,t,u,m)+gv*gv*imAvf_pmp0(beta,t,u,m));
     }
     return sum;
}

double reFf_pmpm(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*reAaf_pmpm(beta,t,u,m)+gv*gv*reAvf_pmpm(beta,t,u,m));
     }
     return sum;
}

double imFf_pmpm(double beta, double t)
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
         sum=sum+alpha*alpha*charge*charge*(ga*ga*imAaf_pmpm(beta,t,u,m)+gv*gv*imAvf_pmpm(beta,t,u,m));
     }
     return sum;
}

//asymptot fermion


double re_dvf_pppm_as(double beta, double t, double u, double m)
{
     return 0;
}

double im_dvf_pppm_as(double beta, double t, double u, double m)
{
     return 0;
}
double re_daf_pppm_as(double beta, double t, double u, double m)
{
     return 0;
}

double im_daf_pppm_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_dvf_pppp_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -4*(pow(log(u/t),2)+Pi*Pi);
}

double im_dvf_pppp_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_pppp_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -4*(pow(log(u/t),2)+Pi*Pi);
}

double im_daf_pppp_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_dvf_pmpp_as(double beta, double t, double u, double m)
{
     return 0;
}

double im_dvf_pmpp_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_pmpp_as(double beta, double t, double u, double m)
{ 
     return 0;
}

double im_daf_pmpp_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_dvf_pp00_as(double beta, double t, double u, double m)
{
     
     return 0;
}

double im_dvf_pp00_as(double beta, double t, double u, double m)
{  
     
     return 0;
}

double re_daf_pp00_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -8*m/mZ2*(pow(log(s/m),2)-Pi*Pi);
}

double im_daf_pp00_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     //printf(" hel %lf\n",(m));
     return -8*m/mZ2*(-2*Pi*log(s/m));
}

double re_dvf_pm00_as(double beta, double t, double u, double m)
{
     return 0;
}

double im_dvf_pm00_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_pm00_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -8*m/mZ2*(t/u*pow(log(-t/s),2)+u/t*pow(log(-u/s),2));
}

double im_daf_pm00_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return -8*m/mZ2*(t/u*2*Pi*log(-t/s)+u/t*2*Pi*log(-u/s));
}

double re_dvf_ppp0_as(double beta, double t, double u, double m)
{
     return 0;
}

double im_dvf_ppp0_as(double beta, double t, double u, double m)
{
     return 0;
}

double re_daf_ppp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -pow(s*mZ2/(2*u*t),0.5)*8*m/mZ2*((u-t)/s*(pow(log(s/m),2)+pow(log(t/u),2))-2*log(s/m)*log(t/u));
}

double im_daf_ppp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -pow(s*mZ2/(2*u*t),0.5)*8*m/mZ2*((u-t)/s*(-2*Pi*log(s/m))+2*Pi*log(t/u));
}

double re_dvf_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(pow(log(-s/t),2));
}

double im_dvf_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(-2*Pi*log(-s/t));
}
double re_daf_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(pow(log(-s/t),2))+pow(s*mZ2/(2*u*t),0.5)*8*m*u/(mZ2*s)*(pow(log(s*m/(t*u)),2));
}

double im_daf_pmp0_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -pow(s*mZ2/(2*u*t),0.5)*8*t/s*(-2*Pi*log(-s/t))+pow(s*mZ2/(2*u*t),0.5)*8*m*u/(mZ2*s)*(-2*Pi*log(s*m/(t*u)));
}

double re_dvf_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -4*pow(log(-s/t),2);
}

double im_dvf_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     
     return 8*Pi*log(-s/t);
}
double re_daf_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return -4*pow(log(-s/t),2);;
}

double im_daf_pmpm_as(double beta, double t, double u, double m)
{
     double s=4*mZ2/(1-beta*beta);
     return 8*Pi*log(-s/t);;
}


double reAvf_pppm_as(double beta, double t, double u, double m)
{
     return -2*reAs_pppm_as(beta,t,u,m)+re_dvf_pppm_as(beta,t,u,m);
}

double imAvf_pppm_as(double beta, double t, double u, double m)
{
     return -2*imAs_pppm_as(beta,t,u,m)+im_dvf_pppm_as(beta,t,u,m);
}

double reAvf_pppp_as(double beta, double t, double u, double m)
{
     return -2*reAs_pppp_as(beta,t,u,m)+re_dvf_pppp_as(beta,t,u,m);
}

double imAvf_pppp_as(double beta, double t, double u, double m)
{
     return -2*imAs_pppp_as(beta,t,u,m)+im_dvf_pppp_as(beta,t,u,m);
}

double reAvf_pmpp_as(double beta, double t, double u, double m)
{
     return -2*reAs_pmpp_as(beta,t,u,m)+re_dvf_pmpp_as(beta,t,u,m);
}

double imAvf_pmpp_as(double beta, double t, double u, double m)
{
     return  -2*imAs_pmpp_as(beta,t,u,m)+im_dvf_pmpp_as(beta,t,u,m);
}

double reAvf_pm00_as(double beta, double t, double u, double m)
{
     return -2*reAs_pm00_as(beta,t,u,m)+re_dvf_pm00_as(beta,t,u,m);
}

double imAvf_pm00_as(double beta, double t, double u, double m)
{
     return -2*imAs_pm00_as(beta,t,u,m)+im_dvf_pm00_as(beta,t,u,m);
}

double reAvf_pp00_as(double beta, double t, double u, double m)
{
     return -2*reAs_pp00_as(beta,t,u,m)+re_dvf_pp00_as(beta,t,u,m);
}

double imAvf_pp00_as(double beta, double t, double u, double m)
{
     return -2*imAs_pp00_as(beta,t,u,m)+im_dvf_pp00_as(beta,t,u,m);
}

double reAvf_ppp0_as(double beta, double t, double u, double m)
{
     return -2*reAs_ppp0_as(beta,t,u,m)+re_dvf_ppp0_as(beta,t,u,m);
}

double imAvf_ppp0_as(double beta, double t, double u, double m)
{
     return -2*imAs_ppp0_as(beta,t,u,m)+im_dvf_ppp0_as(beta,t,u,m);
}

double reAvf_pmp0_as(double beta, double t, double u, double m)
{
     return -2*reAs_pmp0_as(beta,t,u,m)+re_dvf_pmp0_as(beta,t,u,m);
}

double imAvf_pmp0_as(double beta, double t, double u, double m)
{
     return -2*imAs_pmp0_as(beta,t,u,m)+im_dvf_pmp0_as(beta,t,u,m);
}

double reAvf_pmpm_as(double beta, double t, double u, double m)
{
     return -2*reAs_pmpm_as(beta,t,u,m)+re_dvf_pmpm_as(beta,t,u,m);
}

double imAvf_pmpm_as(double beta, double t, double u, double m)
{
     return -2*imAs_pmpm_as(beta,t,u,m)+im_dvf_pmpm_as(beta,t,u,m);
}



double reAaf_pppm_as(double beta, double t, double u, double m)
{
     return -2*reAs_pppm_as(beta,t,u,m)+re_daf_pppm_as(beta,t,u,m);
}

double imAaf_pppm_as(double beta, double t, double u, double m)
{
     return -2*imAs_pppm_as(beta,t,u,m)+im_daf_pppm_as(beta,t,u,m);
}

double reAaf_pppp_as(double beta, double t, double u, double m)
{
     return -2*reAs_pppp_as(beta,t,u,m)+re_daf_pppp_as(beta,t,u,m);
}

double imAaf_pppp_as(double beta, double t, double u, double m)
{
     return -2*imAs_pppp_as(beta,t,u,m)+im_daf_pppp_as(beta,t,u,m);
}

double reAaf_pmpp_as(double beta, double t, double u, double m)
{
     return -2*reAs_pmpp_as(beta,t,u,m)+re_daf_pmpp_as(beta,t,u,m);
}

double imAaf_pmpp_as(double beta, double t, double u, double m)
{
     return  -2*imAs_pmpp_as(beta,t,u,m)+im_daf_pmpp_as(beta,t,u,m);
}

double reAaf_pm00_as(double beta, double t, double u, double m)
{
     return -2*reAs_pm00_as(beta,t,u,m)+re_daf_pm00_as(beta,t,u,m);
}

double imAaf_pm00_as(double beta, double t, double u, double m)
{
     return -2*imAs_pm00_as(beta,t,u,m)+im_daf_pm00_as(beta,t,u,m);
}

double reAaf_pp00_as(double beta, double t, double u, double m)
{
     return -2*reAs_pp00_as(beta,t,u,m)+re_daf_pp00_as(beta,t,u,m);
}

double imAaf_pp00_as(double beta, double t, double u, double m)
{
     return -2*imAs_pp00_as(beta,t,u,m)+im_daf_pp00_as(beta,t,u,m);
}

double reAaf_ppp0_as(double beta, double t, double u, double m)
{
     return -2*reAs_ppp0_as(beta,t,u,m)+re_daf_ppp0_as(beta,t,u,m);
}

double imAaf_ppp0_as(double beta, double t, double u, double m)
{
     return -2*imAs_ppp0_as(beta,t,u,m)+im_daf_ppp0_as(beta,t,u,m);
}

double reAaf_pmp0_as(double beta, double t, double u, double m)
{
     return -2*reAs_pmp0_as(beta,t,u,m)+re_daf_pmp0_as(beta,t,u,m);
}

double imAaf_pmp0_as(double beta, double t, double u, double m)
{
     return -2*imAs_pmp0_as(beta,t,u,m)+im_daf_pmp0_as(beta,t,u,m);
}

double reAaf_pmpm_as(double beta, double t, double u, double m)
{
     return -2*reAs_pmpm_as(beta,t,u,m)+re_daf_pmpm_as(beta,t,u,m);
}

double imAaf_pmpm_as(double beta, double t, double u, double m)
{
     return -2*imAs_pmpm_as(beta,t,u,m)+im_daf_pmpm_as(beta,t,u,m);
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
     for (int i=0;i<size;i++){
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
         fi=1/st*(1+2*mi*ReC(st,mi));
     }
     else if (spin==0.5)
     {
         fi=-2/st*(1-2*mi*(st-1)*ReC(st,mi));
     }
     else if (spin==1)
     {
         fi=1/st*(3+mh2/(2.*mi))-(16*mi-(mh2+6*mi)/st)*ReC(st,mi);
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
         fi=1/st*(2*mi*ImC(st,mi));
     }
     else if (spin==0.5)
     {
         fi=-2/st*(-2*mi*(st-1)*ImC(st,mi));
     }
     else if (spin==1)
     {
         fi=-(16*mi-(mh2+6*mi)/st)*ImC(st,mi);
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


double dsigma_h(double s)
{
     double beta=pow(1-4*mZ2/(s+0.00001),0.5);
    
     double sum=imFh_pppp(beta);
     
     return sum;
}

double dsigma_W(double s, double theta)
{
     double beta=pow(1-4*mZ2/(s+0.00001),0.5);
     //printf(" %lf\n",beta);
     double t=mZ2-s/2*(1-pow(1-4*mZ2/(s+0.0001),0.5)*cos(Pi*theta/180.0));
     double u=2*mZ2-s-t;
     //first with +beta
     double sum=reFw_pppp(beta,t,mW2);
     
     return sum;//*2
}

double dsigma_f(double s, double theta)
{
     double mq=25.000001;
     double charge=-2/3.0;
     double isospin=1/2.0;
     double beta=pow(1-4*mZ2/(s+0.00001),0.5);
     //printf(" masse %lf\n",mq);
     double t=mZ2-s/2*(1-pow(1-4*mZ2/(s+0.0001),0.5)*cos(Pi*theta/180.0));
     double u=2*mZ2-s-t;
     //printf(" t %lf\n",t);
     //first with +beta
     double sum= reFf_pppp_as(beta,  t);
     return sum;//*2
}


double dsigma_asymptot(double s, double t,int exclude_loop)
     double beta=pow(1-4*mZ2/(s+0.00001),0.5);
     double u=2*mZ2-s-t;
     
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

     double F2_pppp=pow(coeff_W*reFw_pppp_as(beta,t,mW2)+coeff_f*reFf_pppp_as(beta,t)+coeff_f*reFf_pppp_as(beta,t),2)+pow(coeff_W*imFw_pppp_as(beta,t,mW2) +coeff_f*imFf_pppp_as(beta,t)+coeff_f*imFf_pppp_as(beta,t),2);

     double F2_pppp_mbeta=pow(coeff_W*reFw_pppp_as(-beta,t,mW2)+coeff_f*reFf_pppp_as(-beta,t)+coeff_f*reFf_pppp_as(-beta,t),2)+pow(coeff_W*imFw_pppp_as(-beta,t,mW2)+coeff_f*imFf_pppp_as(-beta,t)+coeff_f*imFf_pppp_as(-beta,t),2);
     
     double F2_pppm=pow(coeff_W*reFw_pppm_as(beta,t,mW2)+coeff_f*reFf_pppm_as(beta,t)+coeff_f*reFf_pppm_as(beta,t),2)+pow(coeff_W*imFw_pppm_as(beta,t,mW2)+coeff_f*imFf_pppm_as(beta,t)+coeff_f*imFf_pppm_as(beta,t),2);
     double F2_pppm_mbeta=pow(coeff_W*reFw_pppm_as(-beta,t,mW2)+coeff_f*reFf_pppm_as(-beta,t)+coeff_f*reFf_pppm_as(-beta,t),2)+pow(coeff_W*imFw_pppm_as(-beta,t,mW2)+coeff_f*imFf_pppm_as(-beta,t)+coeff_f*imFf_pppm_as(-beta,t),2);

     double F2_pmpm=pow(coeff_W*reFw_pmpm_as(beta,t,mW2)+coeff_f*reFf_pmpm_as(beta,t)+coeff_f*reFf_pmpm_as(beta,t),2)+pow(coeff_W*imFw_pmpm_as(beta,t,mW2)+coeff_f*imFf_pmpm_as(beta,t)+coeff_f*imFf_pmpm_as(beta,t),2);
     double F2_pmpm_mbeta=pow(coeff_W*reFw_pmpm_as(-beta,t,mW2)+coeff_f*reFf_pmpm_as(-beta,t)+coeff_f*reFf_pmpm_as(-beta,t),2)+pow(coeff_W*imFw_pmpm_as(-beta,t,mW2)+coeff_f*imFf_pmpm_as(-beta,t)+coeff_f*imFf_pmpm_as(-beta,t),2);

     double F2_pmpp=pow(coeff_W*reFw_pmpp_as(beta,t,mW2)+coeff_f*reFf_pmpp_as(beta,t)+coeff_f*reFf_pmpp_as(beta,t),2)+pow(coeff_W*imFw_pmpp_as(beta,t,mW2)+coeff_f*imFf_pmpp_as(beta,t)+coeff_f*imFf_pmpp_as(beta,t),2);
     double F2_pmpp_mbeta=pow(coeff_W*reFw_pmpp_as(-beta,t,mW2)+coeff_f*reFf_pmpp_as(-beta,t)+coeff_f*reFf_pmpp_as(-beta,t),2)+pow(coeff_W*imFw_pmpp_as(-beta,t,mW2)+coeff_f*imFf_pmpp_as(-beta,t)+coeff_f*imFf_pmpp_as(-beta,t),2);

     double F2_pp00=pow(coeff_W*reFw_pp00_as(beta,t,mW2)+coeff_f*reFf_pp00_as(beta,t)+coeff_f*reFf_pp00_as(beta,t),2)+pow(coeff_W*imFw_pp00_as(beta,t,mW2)+coeff_f*imFf_pp00_as(beta,t)+coeff_f*imFf_pp00_as(beta,t),2);
     double F2_pm00=pow(coeff_W*reFw_pm00_as(beta,t,mW2)+coeff_f*reFf_pm00_as(beta,t)+coeff_f*reFf_pm00_as(beta,t),2)+pow(coeff_W*imFw_pm00_as(beta,t,mW2)+coeff_f*imFf_pm00_as(beta,t)+coeff_f*imFf_pm00_as(beta,t),2);    

     double F2_ppp0=pow(coeff_W*reFw_ppp0_as(beta,t,mW2)+coeff_f*reFf_ppp0_as(beta,t)+coeff_f*reFf_ppp0_as(beta,t),2)+pow(coeff_W*imFw_ppp0_as(beta,t,mW2)+coeff_f*imFf_ppp0_as(beta,t)+coeff_f*imFf_ppp0_as(beta,t),2);
     double F2_ppp0_mbeta=pow(coeff_W*reFw_ppp0_as(-beta,t,mW2)+coeff_f*reFf_ppp0_as(-beta,t)+coeff_f*reFf_ppp0_as(-beta,t),2)+pow(coeff_W*imFw_ppp0_as(-beta,t,mW2)+coeff_f*imFf_ppp0_as(-beta,t)+coeff_f*imFf_ppp0_as(-beta,t),2);
     double F2_ppp0_u=pow(coeff_W*reFw_ppp0_as(beta,u,mW2)+coeff_f*reFf_ppp0_as(beta,u)+coeff_f*reFf_ppp0_as(beta,u),2)+pow(coeff_W*imFw_ppp0_as(beta,u,mW2)+coeff_f*imFf_ppp0_as(beta,u)+coeff_f*imFf_ppp0_as(beta,u),2);
     double F2_ppp0_mbeta_u=pow(coeff_W*reFw_ppp0_as(-beta,u,mW2)+coeff_f*reFf_ppp0_as(-beta,u)+coeff_f*reFf_ppp0_as(-beta,u),2)+pow(coeff_W*imFw_ppp0_as(-beta,u,mW2)+coeff_f*imFf_ppp0_as(-beta,u)+coeff_f*imFf_ppp0_as(-beta,u),2);

     double F2_pmp0=pow(coeff_W*reFw_pmp0_as(beta,t,mW2)+coeff_f*reFf_pmp0_as(beta,t)+coeff_f*reFf_pmp0_as(beta,t),2)+pow(coeff_W*imFw_pmp0_as(beta,t,mW2)+coeff_f*imFf_pmp0_as(beta,t)+coeff_f*imFf_pmp0_as(beta,t),2);
     double F2_pmp0_mbeta=pow(coeff_W*reFw_pmp0_as(-beta,t,mW2)+coeff_f*reFf_pmp0_as(-beta,t)+coeff_f*reFf_pmp0_as(-beta,t),2)+pow(coeff_W*imFw_pmp0_as(-beta,t,mW2)+coeff_f*imFf_pmp0_as(-beta,t)+coeff_f*imFf_pmp0_as(-beta,t),2);
     double F2_pmp0_u=pow(coeff_W*reFw_pmp0_as(beta,u,mW2)+coeff_f*reFf_pmp0_as(beta,u)+coeff_f*reFf_pmp0_as(beta,u),2)+pow(coeff_W*imFw_pmp0_as(beta,u,mW2)+coeff_f*imFf_pmp0_as(beta,u)+coeff_f*imFf_pmp0_as(beta,u),2);
     double F2_pmp0_mbeta_u=pow(coeff_W*reFw_pmp0_as(-beta,u,mW2)+coeff_f*reFf_pmp0_as(-beta,u)+coeff_f*reFf_pmp0_as(-beta,u),2)+pow(coeff_W*imFw_pmp0_as(-beta,u,mW2)+coeff_f*imFf_pmp0_as(-beta,u)+coeff_f*imFf_pmp0_as(-beta,u),2);

     double sum_X=F2_pppp+F2_pppp_mbeta+F2_pppm+F2_pppm_mbeta+F2_pmpp+F2_pmpp_mbeta+F2_pmpm+F2_pmpm_mbeta;
     double sum_Z=F2_pp00+F2_pm00;
     double sum_Y=F2_ppp0+F2_ppp0_mbeta+F2_pmp0+F2_pmp0_mbeta +F2_ppp0_u+F2_ppp0_mbeta_u+F2_pmp0_u+F2_pmp0_mbeta_u;

     double sum=beta/(128.0*s*Pi)*(2*sum_X+2*sum_Z+2*sum_Y);
     
double dsigma_ZZ(double s, double t,int exclude_loop)
{    
     if (s<4*mZ2)
     {
          return 0;
     }
     if (t>0)
     {
          return 0;
     }
     if (2*mZ2-s-t>0)
     {
          return 0;
     }
     if (t*(2*mZ2-s-t)<mZ2*mZ2)
     {
          return 0;
     }
     if (s>1000*4*mZ2)//asymptotic regime
     {
          return dsigma_asymptot(s,t,exclude_loop);
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

     double beta=pow(1-4*mZ2/(s+0.00001),0.5);
     //printf(" masse %lf\n",mq);
     //double t=mZ2-s/2*(1-pow(1-4*mZ2/(s+0.0001),0.5)*cos(Pi*theta/180.0));
     double u=2*mZ2-s-t;
     

     double F2_pppp=pow(coeff_Higgs*reFh_pppp(beta)+coeff_W*reFw_pppp(beta,t,mW2)+coeff_f*reFf_pppp(beta,t)+coeff_f*reFf_pppp(beta,t),2)+pow(coeff_Higgs*imFh_pppp(beta)+coeff_W*imFw_pppp(beta,t,mW2) +coeff_f*imFf_pppp(beta,t)+coeff_f*imFf_pppp(beta,t),2);

     double F2_pppp_mbeta=pow(coeff_Higgs*reFh_pppp(-beta)+coeff_W*reFw_pppp(-beta,t,mW2)+coeff_f*reFf_pppp(-beta,t)+coeff_f*reFf_pppp(-beta,t),2)+pow(coeff_Higgs*imFh_pppp(-beta)+coeff_W*imFw_pppp(-beta,t,mW2)+coeff_f*imFf_pppp(-beta,t)+coeff_f*imFf_pppp(-beta,t),2);
     
     double F2_pppm=pow(coeff_W*reFw_pppm(beta,t,mW2)+coeff_f*reFf_pppm(beta,t)+coeff_f*reFf_pppm(beta,t),2)+pow(coeff_W*imFw_pppm(beta,t,mW2)+coeff_f*imFf_pppm(beta,t)+coeff_f*imFf_pppm(beta,t),2);
     double F2_pppm_mbeta=pow(coeff_W*reFw_pppm(-beta,t,mW2)+coeff_f*reFf_pppm(-beta,t)+coeff_f*reFf_pppm(-beta,t),2)+pow(coeff_W*imFw_pppm(-beta,t,mW2)+coeff_f*imFf_pppm(-beta,t)+coeff_f*imFf_pppm(-beta,t),2);

     double F2_pmpm=pow(coeff_W*reFw_pmpm(beta,t,mW2)+coeff_f*reFf_pmpm(beta,t)+coeff_f*reFf_pmpm(beta,t),2)+pow(coeff_W*imFw_pmpm(beta,t,mW2)+coeff_f*imFf_pmpm(beta,t)+coeff_f*imFf_pmpm(beta,t),2);
     double F2_pmpm_mbeta=pow(coeff_W*reFw_pmpm(-beta,t,mW2)+coeff_f*reFf_pmpm(-beta,t)+coeff_f*reFf_pmpm(-beta,t),2)+pow(coeff_W*imFw_pmpm(-beta,t,mW2)+coeff_f*imFf_pmpm(-beta,t)+coeff_f*imFf_pmpm(-beta,t),2);

     double F2_pmpp=pow(coeff_W*reFw_pmpp(beta,t,mW2)+coeff_f*reFf_pmpp(beta,t)+coeff_f*reFf_pmpp(beta,t),2)+pow(coeff_W*imFw_pmpp(beta,t,mW2)+coeff_f*imFf_pmpp(beta,t)+coeff_f*imFf_pmpp(beta,t),2);
     double F2_pmpp_mbeta=pow(coeff_W*reFw_pmpp(-beta,t,mW2)+coeff_f*reFf_pmpp(-beta,t)+coeff_f*reFf_pmpp(-beta,t),2)+pow(coeff_W*imFw_pmpp(-beta,t,mW2)+coeff_f*imFf_pmpp(-beta,t)+coeff_f*imFf_pmpp(-beta,t),2);

     double F2_pp00=pow(coeff_Higgs*reFh_pp00(beta)+coeff_W*reFw_pp00(beta,t,mW2)+coeff_f*reFf_pp00(beta,t)+coeff_f*reFf_pp00(beta,t),2)+pow(coeff_Higgs*imFh_pp00(beta)+coeff_W*imFw_pp00(beta,t,mW2)+coeff_f*imFf_pp00(beta,t)+coeff_f*imFf_pp00(beta,t),2);
     double F2_pm00=pow(coeff_W*reFw_pm00(beta,t,mW2)+coeff_f*reFf_pm00(beta,t)+coeff_f*reFf_pm00(beta,t),2)+pow(coeff_W*imFw_pm00(beta,t,mW2)+coeff_f*imFf_pm00(beta,t)+coeff_f*imFf_pm00(beta,t),2);    

     double F2_ppp0=pow(coeff_W*reFw_ppp0(beta,t,mW2)+coeff_f*reFf_ppp0(beta,t)+coeff_f*reFf_ppp0(beta,t),2)+pow(coeff_W*imFw_ppp0(beta,t,mW2)+coeff_f*imFf_ppp0(beta,t)+coeff_f*imFf_ppp0(beta,t),2);
     double F2_ppp0_mbeta=pow(coeff_W*reFw_ppp0(-beta,t,mW2)+coeff_f*reFf_ppp0(-beta,t)+coeff_f*reFf_ppp0(-beta,t),2)+pow(coeff_W*imFw_ppp0(-beta,t,mW2)+coeff_f*imFf_ppp0(-beta,t)+coeff_f*imFf_ppp0(-beta,t),2);
     double F2_ppp0_u=pow(coeff_W*reFw_ppp0(beta,u,mW2)+coeff_f*reFf_ppp0(beta,u)+coeff_f*reFf_ppp0(beta,u),2)+pow(coeff_W*imFw_ppp0(beta,u,mW2)+coeff_f*imFf_ppp0(beta,u)+coeff_f*imFf_ppp0(beta,u),2);
     double F2_ppp0_mbeta_u=pow(coeff_W*reFw_ppp0(-beta,u,mW2)+coeff_f*reFf_ppp0(-beta,u)+coeff_f*reFf_ppp0(-beta,u),2)+pow(coeff_W*imFw_ppp0(-beta,u,mW2)+coeff_f*imFf_ppp0(-beta,u)+coeff_f*imFf_ppp0(-beta,u),2);

     double F2_pmp0=pow(coeff_W*reFw_pmp0(beta,t,mW2)+coeff_f*reFf_pmp0(beta,t)+coeff_f*reFf_pmp0(beta,t),2)+pow(coeff_W*imFw_pmp0(beta,t,mW2)+coeff_f*imFf_pmp0(beta,t)+coeff_f*imFf_pmp0(beta,t),2);
     double F2_pmp0_mbeta=pow(coeff_W*reFw_pmp0(-beta,t,mW2)+coeff_f*reFf_pmp0(-beta,t)+coeff_f*reFf_pmp0(-beta,t),2)+pow(coeff_W*imFw_pmp0(-beta,t,mW2)+coeff_f*imFf_pmp0(-beta,t)+coeff_f*imFf_pmp0(-beta,t),2);
     double F2_pmp0_u=pow(coeff_W*reFw_pmp0(beta,u,mW2)+coeff_f*reFf_pmp0(beta,u)+coeff_f*reFf_pmp0(beta,u),2)+pow(coeff_W*imFw_pmp0(beta,u,mW2)+coeff_f*imFf_pmp0(beta,u)+coeff_f*imFf_pmp0(beta,u),2);
     double F2_pmp0_mbeta_u=pow(coeff_W*reFw_pmp0(-beta,u,mW2)+coeff_f*reFf_pmp0(-beta,u)+coeff_f*reFf_pmp0(-beta,u),2)+pow(coeff_W*imFw_pmp0(-beta,u,mW2)+coeff_f*imFf_pmp0(-beta,u)+coeff_f*imFf_pmp0(-beta,u),2);

     double sum_X=F2_pppp+F2_pppp_mbeta+F2_pppm+F2_pppm_mbeta+F2_pmpp+F2_pmpp_mbeta+F2_pmpm+F2_pmpm_mbeta;
     double sum_Z=F2_pp00+F2_pm00;
     double sum_Y=F2_ppp0+F2_ppp0_mbeta+F2_pmp0+F2_pmp0_mbeta +F2_ppp0_u+F2_ppp0_mbeta_u+F2_pmp0_u+F2_pmp0_mbeta_u;

     double sum=beta/(128.0*s*Pi)*(2*sum_X+2*sum_Z+2*sum_Y);
     return sum;
}

























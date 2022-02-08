#include<math.h>
#include <stdio.h>
#include<gsl/gsl_sf.h>


const double Pi = 4 * atan(1);
double mZ2=91.2*91.2;//Mz/(4.0*mi);
const double epsilon=0.000001;

double reaz=4.76;
double imaz=0;

double ReBz(double s,double m)//use for mw az()=1.749
{
     if (s>0)
     {
     return -log(4*s)+1.749;
     }
     else
     {
     return -log(-4*s)+1.749;
     }
}
double ImBz(double s,double m)
{
     if (s>0)
     {
     return Pi;
     }
     else
     {
     return 0;
     }
}

double ReC(double s,double m)
{
    return 1/(4.0*m*2.0*s)*(pow(log(4*s),2)-Pi*Pi);
}
double ImC(double s,double m)
{
    return 1/(4.0*m*s)*(-Pi)*log(4*s);
}

double ReCz(double u,double m)
{
    return 1/(4.0*m*u)*(1/2.0*pow(log(-4*u),2)+reaz/2.0);
}
double ImCz(double u,double m)
{
    return 1/(4.0*m*u)*(imaz/2.0);
}

double ReCzz(double s,double m)
{
    return ReC(s,m)+1/(2.0*4.0*m*s)*reaz;
}
double ImCzz(double s,double m)
{
    return ImC(s,m)+1/(2.0*4.0*m*s)*imaz;
}

double ReDut(double u, double t,double m)
{
    return 2/(4.0*m*4.0*m*u*t)*(log(-4*u)*log(-4*t)-Pi*Pi/2.0+reaz);
}


double ImDut(double u, double t,double m)
{
    return 2/(4.0*m*4.0*m*u*t)*(imaz);
}

double ReDst(double s, double t,double m)
{
    return 2/(4.0*m*4.0*m*s*t)*(log(4*s)*log(-4*t)-Pi*Pi/2.0+reaz);
}
double ImDst(double s, double t,double m)
{
    return 2/(4.0*m*4.0*m*s*t)*(-Pi*log(-4*t)+imaz);
}

double ReF(double s, double t, double u,double m)
{
    return 2/(4.0*m*4.0*m*s*u)*log(4*s)*log(-4*u)+2/(4.0*m*4.0*m*s*t)*log(4*s)*log(-4*t)+2/(4.0*m*4.0*m*t*u)*log(-4*u)*log(-4*t);
}
double ImF(double s, double t, double u,double m)
{
    return -2*Pi/(4.0*m*s)*(1/(4.0*m*u)*log(-4*u)+1/(4.0*m*t)*log(-4*t));
}

double ReE1(double s, double t,double m)
{
    return Pi*Pi-reaz+pow(log(-4*t),2)-2*log(4*s)*log(-4*t);
}

double ImE1(double s, double t,double m)
{
    return -imaz+2*Pi*log(-4*t);
}

double ReE2(double t, double u,double m)
{
    return Pi*Pi+pow(log(-4*t)-log(-4*u),2);
}
double ImE2(double t, double u,double m)
{
    return 0;
}

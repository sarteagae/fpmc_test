
#include<math.h>
#include <stdio.h>
#include<gsl/gsl_sf.h>



const double PI = 4 * atan(1);
double mZ2=91.2*91.2;
const double epsilon=0.000001;


double mod(double theta,double Modulo) //return theta into [-Modulo/2;+Modulo/2[ We assume explicitely the epsilon prescritpion here
{
    double thetap=fmod(theta,Modulo); //so here we are between -Modulo and +Modulo
    if (thetap>=Modulo/2.0)
    {
        thetap=thetap-Modulo;
    }
    if (thetap<-Modulo/2.0)
    {
        thetap=thetap+Modulo;
    }
    return thetap;
}

double reb(double sm)//take into account the epsilon presc; =sqrt(r) cos theta/2 ; s has no dim sm=s/4m
{
    double r= sqrt(pow(1-1/sm,2)+16*epsilon*epsilon/(sm*sm));
    double theta= mod(atan2(4*epsilon/sm,(1-1/sm)),2*PI);
    return sqrt(r)*cos(theta/2.0);
}

double imb(double sm)//take into account the epsilon presc; =sqrt(r) cos theta/2 ; s has no dim sm=s/4m
{
    double r= sqrt(pow(1-1/sm,2)+16*epsilon*epsilon/(sm*sm));
    double theta= mod(atan2(4*epsilon/sm,(1-1/sm)),2*PI);
    return sqrt(r)*sin(theta/2.0);
}
double rebz(double a,double b)//take into account the epsilon presc; =sqrt(r) cos theta/2 ; s has no dim sm=s/4m
{
    double r= sqrt(pow(1-4*a/b,2)+16*epsilon*epsilon/(b*b));
    double theta= mod(atan2(4*epsilon/b,(1-4*a/b) ),2*PI);
    return sqrt(r)*cos(theta/2.0);
}

double imbz(double a,double b)//usefull for spence function; take into account the epsilon presc; =sqrt(r) cos theta/2 ; a,b has no dim :/4m
{
    double r= sqrt(pow(1-4*a/b,2)+16*epsilon*epsilon/(b*b));
    double theta= mod(atan2(4*epsilon/b,(1-4*a/b)),2*PI);
    return sqrt(r)*sin(theta/2.0);
}

double beta(double sm,double m2)
{
    return sqrt(1-mZ2/(m2*sm));
}

double re_logz(double sm)
{
    double rebs=reb(sm);
    double imbs=imb(sm);
    
    return 0.5*log((pow(1+rebs,2)+pow(imbs,2))/(pow(1.0-rebs,2)+pow(imbs,2)));// teh 0.5 comes from the sqrt inside the log
}

double im_logz(double sm)
{
    double rebs=reb(sm);
    double imbs=imb(sm);
   
    return mod(atan2(-imbs,(-1-rebs))-atan2(-imbs,(1-rebs)),2*PI);
}


double reI1(double a, double b, double rey, double imy)//sum of 2 dilogs
{
    if (a==0)
    {
       printf("error I1");
    }
    
    double r;
    double theta;
    double rez = -b/a;
    double imz=epsilon/a;

    double resum ; // the real part of the sum
    //double imsum ;// the imaginary part
    gsl_sf_result result_re;
    gsl_sf_result result_im;

    r = sqrt(rey * rey + imy * imy) / sqrt(pow(rey -rez, 2) + pow(imy -imz, 2));
    theta = mod(atan2(imy, rey) - atan2(imy-imz , rey -rez),2*PI);

    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z1) )
    resum = result_re.val;
    //imsum = result_im.val;
    
    r = sqrt((rey-1) * (rey-1) + imy * imy) / sqrt(pow(rey -rez, 2) + pow(imy -imz, 2));
    theta = mod(atan2(imy, rey-1) - atan2(imy-imz , rey - rez),2*PI);

    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z1) )
    resum -= result_re.val;
    //imsum += result_im.val;
    return resum;
}

double imI1(double a, double b, double rey, double imy)//sum of 2 dilogs
{
    if (a==0 )
    {
       printf("error I1");
    }
    
    double r;
    double theta;
    double rez = -b/a;
    double imz=epsilon/a;

    //double resum ; // the real part of the sum
    double imsum ;// the imaginary part
    gsl_sf_result result_re;
    gsl_sf_result result_im;


    r = sqrt(rey * rey + imy * imy) / sqrt(pow(rey -rez, 2) + pow(imy -imz, 2));
    theta = mod(atan2(imy, rey) - atan2(imy-imz , rey - rez),2*PI);
    //printf("r %lf \n", r);
    //printf("theta %lf \n", theta);
    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z1) )
    //resum = result_re.val;
    imsum = result_im.val;
    //printf("resulti1 %lf \n", result_im.val);
    
    r = sqrt((rey-1) * (rey-1) + imy * imy) / sqrt(pow(rey -rez, 2) + pow(imy -imz, 2));
    theta = mod(atan2(imy, rey-1) - atan2(imy-imz , rey - rez),2*PI);
    //printf("r %lf \n", r);
    //printf("theta %lf \n", theta);
    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z1) )
    //resum += result_re.val;
    imsum -= result_im.val;
    //printf("resulti1 %lf \n", result_im.val);
    //printf("result imsumi1   %lf \n", imsum);
    return imsum;
}

double test()
{
    double r=4.24;
    double theta=0.785;

    

    double resum ; // the real part of the sum
    //double imsum ;// the imaginary part
    gsl_sf_result result_re;
    gsl_sf_result result_im;


    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z1) )
    resum = result_re.val;
    gsl_sf_complex_dilog_e(2*r, theta, &result_re, &result_im);//+Li2( (y-1)/(y-z1) )
    resum -= result_re.val;
    return resum;
}
double reI2(double a,double b,double rey,double imy)//sum of 4 dilog with 4 arguments different
{

    if ( b==0)
    {
       printf("error I2");
    }
    
    double r;
    double theta;

    double rez1 = 0.5 * (1 + rebz(a,b));
    double imz1 = 0.5 * imbz(a,b);

    double rez2 = 1-rez1;
    double imz2 = -imz1;

    double resum ; // the real part of the sum
    //double imsum ;// the imaginary part
    gsl_sf_result result_re;
    gsl_sf_result result_im;

    r = sqrt(rey * rey + imy * imy) / sqrt(pow(rey - rez1, 2) + pow(imy - imz1, 2));
    theta = mod(atan2(imy, rey) - atan2(imy - imz1, rey - rez1),2*PI);

    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z1) )
    resum = result_re.val;
    //imsum = result_im.val;

    r = sqrt((rey-1) * (rey-1) + imy * imy) / sqrt(pow(rey - rez1, 2) + pow(imy - imz1, 2));
    theta = mod(atan2(imy, rey-1) - atan2(imy - imz1, rey - rez1),2*PI);

    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( (y-1)/(y-z1) )
    resum -= result_re.val;
    //imsum -= result_im.val;

    //now same thing but with z2
    r = sqrt(rey * rey + imy * imy) / sqrt(pow(rey - rez2, 2) + pow(imy - imz2, 2));
    theta = mod(atan2(imy, rey) - atan2(imy - imz2, rey - rez2),2*PI);

    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z2) )
    resum += result_re.val;
    //imsum += result_im.val;

    r = sqrt((rey - 1) * (rey - 1) + imy * imy) / sqrt(pow(rey - rez2, 2) + pow(imy - imz2, 2));
    theta = mod(atan2(imy, rey - 1) - atan2(imy - imz2, rey - rez2),2*PI);

    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( (y-1)/(y-z2) )
    resum -= result_re.val;
    //imsum -= result_im.val;

    return resum;
}


double imI2(double a,double b,double rey,double imy)//sum of 4 dilog with 4 arguments different
{
    if ( b==0)
    {
       printf("error I2");
    }

    double r;
    double theta;
    //printf("\n a %lf \n", a);
    //printf("b %lf \n\n", b);
    double rez1 = 0.5 * (1 + rebz(a,b));
    double imz1 = 0.5 * imbz(a,b);
    //printf("\n z1 %lf \n", rez1);
    //printf("z1 %lf \n\n", imz1);
    double rez2 = 1-rez1;
    double imz2 = -imz1;
    //printf("\n z2 %lf \n", rez2);
    //printf("z2 %lf \n\n", imz2);
    //double resum ; // the real part of the sum
    double imsum ;// the imaginary part
    gsl_sf_result result_re;
    gsl_sf_result result_im;


    r = sqrt(rey * rey + imy * imy) / sqrt(pow(rey - rez1, 2) + pow(imy - imz1, 2));
    theta = mod(atan2(imy, rey) - atan2(imy - imz1, rey - rez1),2*PI);
    //printf("r1 %lf \n", r);
    //printf("theta %lf \n", theta);
    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z1) )
    //resum = result_re.val;
    imsum = result_im.val;
    //printf("result %lf \n", result_im.val);

    r = sqrt((rey-1) * (rey-1) + imy * imy) / sqrt(pow(rey - rez1, 2) + pow(imy - imz1, 2));
    theta = mod( atan2(imy, rey-1) - atan2(imy - imz1, rey - rez1),2*PI);
    //printf("r2 %lf \n", r);
    //printf("theta %lf \n", theta);
    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( (y-1)/(y-z1) )
    //resum -= result_re.val;
    imsum -= result_im.val;
    //printf("result %lf \n", result_im.val);
    
    //now same thing but with z2
    r = sqrt(rey * rey + imy * imy) / sqrt(pow(rey - rez2, 2) + pow(imy - imz2, 2));
    theta = mod(atan2(imy, rey) - atan2(imy - imz2, rey - rez2),2*PI);
    //printf("r3 %lf \n", r);
    //printf("theta %lf \n", theta);
    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( y/(y-z2) f)
    //resum += result_re.val;
    imsum += result_im.val;
    //printf("result %lf \n", result_im.val);

    r = sqrt((rey - 1) * (rey - 1) + imy * imy) / sqrt(pow(rey - rez2, 2) + pow(imy - imz2, 2));
    theta = mod( atan2(imy, rey - 1) - atan2(imy - imz2, rey - rez2),2*PI);
    //printf("r4 %lf \n", r);
    //printf("theta %lf \n", theta);
    gsl_sf_complex_dilog_e(r, theta, &result_re, &result_im);//+Li2( (y-1)/(y-z2) )
    //resum -= result_re.val;
    imsum -= result_im.val;
    //printf("result %lf \n", result_im.val);
    //printf("result imsumi2   %lf \n", imsum);
    return imsum;
}

double re_log_dut(double eps,double x)//eps is a small param and x is a real one; =re (log (1+x+i eps))
{
    return 0.5*log(pow(1+x,2)+eps*eps);//0.5 for the sqrt 
}

double im_log_dut(double eps,double x)//eps is a small param and x is a real one; =im (log (1+x+i eps))
{
    return mod(atan2(eps,(1+x)),2*PI);
}
//functions Passarino Veltman


double ReB1(double sm)
{   
    return -reb(sm) * re_logz(sm) + imb(sm) * im_logz(sm);
}

double ImB1(double sm)
{   
    return -reb(sm) * im_logz(sm) - imb(sm) * re_logz(sm);
}

double ReBz(double sm, double m2)
{
    return ReB1(sm) - ReB1(mZ2/(4.0*m2));
}

double ImBz(double sm, double m2)
{
    return  ImB1(sm) - ImB1(mZ2/(4.0*m2));
}


double ReC(double sm, double m2)
{
    //printf("dans la fonction %lf\n", 1 / (8.0 * m+0.0001) * 1 / (s+0.0001) * (pow(logabs, 2) - pow(argb, 2))/2.0);
    double ReCs = 1 / (8.0 * m2*sm) * (pow(re_logz(sm), 2) - pow(im_logz(sm), 2));
    return ReCs;
}

double ImC(double sm, double m2)
{
    double ImCs= 1 / (8.0 * m2*sm)* 2 * re_logz(sm) * im_logz(sm);
    return ImCs;
}

double ReCz(double um, double m2)
{
    double ReCzu;
    double alph=1.0;

    ReCzu = 1/ (um - mZ2/(4.0*m2)) * (um * ReC(um,  m2) - alph*mZ2/(4.0*m2) * ReC(mZ2/(4.0*m2), m2));
    return ReCzu;
}
double ImCz(double um, double m2)
{
    double ImCzu;
    double alph=1.0;
    ImCzu = 1/ (um - mZ2/(4.0*m2)) * (um * ImC(um,m2) - alph*mZ2/(4.0*m2) * ImC(mZ2/(4.0*m2),m2));
    return ImCzu;
}

double ReCzz(double sm, double m2)
{
    double betas=beta(sm,m2);
    double alph=-(1+betas)/(1.0-betas);
    
    double y1=alph/(1.0+alph);
    double y2=alph*alph/(alph*alph-1.0);
    
    return (1/(4.0*m2*sm*betas)*(2*reI2(0.25,mZ2/(4.0*m2),y1,0)-reI2(0.25,sm,y2,0)));
}

double ImCzz(double sm, double m2)
{
    double betas=beta(sm,m2);
    double alph=-(1+betas)/(1.0-betas);
    
    double y1=alph/(1.0+alph);
    double y2=alph*alph/(alph*alph-1.0);
    
    return (1/(4.0*m2*sm*betas)*(2*imI2(0.25,mZ2/(4.0*m2),y1,0)-imI2(0.25,sm,y2,0)));
}


double ReDut(double um, double tm, double m2)//every variable is normalized
{
    double sm=2*mZ2/(4.0*m2)-um-tm;
    double Ym=um*tm-mZ2/(4.0*m2)*mZ2/(4.0*m2);
    double R1=pow(1+sm/Ym,0.5);

    return (1/(8.0*m2*m2*Ym*R1)*(2*reI2(0.25,mZ2/(4.0*m2),0.5*(1+R1),0)-reI2(0.25,um,0.5*(1+R1),0)-reI2(0.25,tm,0.5*(1+R1),0)+log((R1-1)/(R1+1.0))*(2*re_log_dut(-4*epsilon,mZ2*sm/(4.0*m2*Ym))-re_log_dut(-4*epsilon,um*sm/(Ym))-re_log_dut(-4*epsilon,tm*sm/(Ym)))));    
}

double ImDut(double um, double tm, double m2)//every variable is normalized
{
    double sm=2*mZ2/(4.0*m2)-um-tm;
    //printf("sm %lf",sm);
    double Ym=um*tm-mZ2/(4.0*m2)*mZ2/(4.0*m2);
    double R1=pow(1+sm/Ym,0.5);
    //printf("sum %lf \n",(8.0*m2*m2*Ym*R1));
    return (1/(8.0*m2*m2*Ym*R1))*(2*imI2(0.25,mZ2/(4.0*m2),0.5*(1+R1),0)-imI2(0.25,um,0.5*(1+R1),0)-imI2(0.25,tm,0.5*(1+R1),0)+log((R1-1)/(R1+1.0))*(2*im_log_dut(-4*epsilon,mZ2*sm/(4.0*m2*Ym))-im_log_dut(-4*epsilon,um*sm/(Ym))-im_log_dut(-4*epsilon,tm*sm/(Ym))));    
}

double ReDst(double sm, double tm, double m2)
{
    double um=2*mZ2/(4.0*m2)-sm-tm;
    double t1=tm-mZ2/(4.0*m2);
    double u1=um-mZ2/(4.0*m2);
    //printf("u %lf\n",u);
    
    
    double Ym=um*tm-mZ2/(4.0*m2)*mZ2/(4.0*m2);
    //printf("%lf\n",u*Y);
    double R2=pow(1+Ym/(sm*tm*tm),0.5);
    //printf("%lf\n",R2);
    
    double betas= beta(sm,m2);
    double beta1=tm/t1;
    double beta2=0.5*(1-betas);
    //printf("beta1  %lf\n",beta1);
    //printf("rebeta2  %lf\n",rebeta2);
    //printf("imbs  %lf\n",imbs);
    //printf("imbeta2  %lf\n\n",imbeta2);
    
    
    double y1p=-tm/(2.0*t1)*(1+R2);
    //printf("rey1p%lf\n",rey1p);
    double y1m=-tm/(2.0*t1)*(1-R2);
    //printf("%lf\n",rey1m);
    
    double y2p=(-mZ2/(4.0*m2)+t1*betas-tm*R2)/(2.0*t1+sm*(1+betas));
    //printf("rey2p%lf\n",rey2p);
    double y2m=(-mZ2/(4.0*m2)+t1*betas+tm*R2)/(2.0*t1+sm*(1+betas));
    
    
    //printf("tesr %lf\n",mZ2/(4.0*m2));
     
    
    double sumreI1=-reI1(sm,0,beta1+y1p,0)+reI1(sm,0,beta1+y1m,0)+reI1(-u1,-t1,y1p/(1-beta1),0)-reI1(-u1,-t1,y1m/(1-beta1),0)-reI1(t1,-t1,-y1p/beta1,0)+reI1(t1,-t1,-y1m/beta1,0)+reI1(sm,0,beta2+y2p,0)-reI1(sm,0,beta2+y2m,0)-reI1(-u1,-t1,y2p/(1-beta2),0)+reI1(-u1,-t1,y2m/(1-beta2),0)+reI1(t1,-t1,-y2p/beta2,0)-reI1(t1,-t1,-y2m/beta2,0);
        
    double sumreI2=-reI2(0.25,mZ2/(4.0*m2),y1p/(1-beta1),0)+reI2(0.25,mZ2/(4.0*m2),y1m/(1-beta1),0)+reI2(0.25,tm,-y1p/beta1,0)-reI2(0.25,tm,-y1m/beta1,0) - reI2(0.25,sm,y2p+beta2,0)+reI2(0.25,sm,y2m+beta2,0)+reI2(0.25,mZ2/(4.0*m2),y2p/(1-beta2),0)-reI2(0.25,mZ2/(4.0*m2),y2m/(1-beta2),0)-reI2(0.25,mZ2/(4.0*m2),-y2p/beta2,0)+reI2(0.25,mZ2/(4.0*m2),-y2m/beta2,0);
    
    double redst=-1/(16.0*m2*m2*sm*tm*R2)*(sumreI1+sumreI2);
    return redst;
    
    
}



double ImDst(double sm, double tm, double m2)
{
    double um=2*mZ2/(4.0*m2)-sm-tm;
    double t1=tm-mZ2/(4.0*m2);
    double u1=um-mZ2/(4.0*m2);
    //printf("u %lf\n",um);
    //printf("t %lf\n",tm);
    //printf("s %lf\n",sm);
    //printf("u %lf\n",u1);
    
    
    double Ym=um*tm-mZ2/(4.0*m2)*mZ2/(4.0*m2);
    //printf("ym%lf\n",Ym);
    double R2=pow(1+Ym/(sm*tm*tm),0.5);
    //printf("%lf\n",R2);
    
    double betas= beta(sm,m2);
    double beta1=tm/t1;
    double beta2=0.5*(1-betas);
    //printf("beta1  %lf\n",beta1);
    //printf("rebeta2  %lf\n",rebeta2);
    //printf("imbs  %lf\n",imbs);
    //printf("imbeta2  %lf\n\n",imbeta2);
    
    
    double y1p=-tm/(2.0*t1)*(1+R2);
    //printf("rey1p%lf\n",rey1p);
    double y1m=-tm/(2.0*t1)*(1-R2);
    //printf("%lf\n",rey1m);
    
    double y2p=(-mZ2/(4.0*m2)+t1*betas-tm*R2)/(2.0*t1+sm*(1+betas));
    //printf("rey2p%lf\n",rey2p);
    double y2m=(-mZ2/(4.0*m2)+t1*betas+tm*R2)/(2.0*t1+sm*(1+betas));
    //printf("y2p+beta2%lf\n",beta2);
     
    
    double sumimI1=-imI1(sm,0,beta1+y1p,0)+imI1(sm,0,beta1+y1m,0)+imI1(-u1,-t1,y1p/(1-beta1),0)-imI1(-u1,-t1,y1m/(1-beta1),0)-imI1(t1,-t1,-y1p/beta1,0)+imI1(t1,-t1,-y1m/beta1,0)+imI1(sm,0,beta2+y2p,0)-imI1(sm,0,beta2+y2m,0)-imI1(-u1,-t1,y2p/(1-beta2),0)+imI1(-u1,-t1,y2m/(1-beta2),0)+imI1(t1,-t1,-y2p/beta2,0)-imI1(t1,-t1,-y2m/beta2,0);
        
    double sumimI2=-imI2(0.25,mZ2/(4.0*m2),y1p/(1-beta1),0)+imI2(0.25,mZ2/(4.0*m2),y1m/(1-beta1),0)+imI2(0.25,tm,-y1p/beta1,0)-imI2(0.25,tm,-y1m/beta1,0) - imI2(0.25,sm,y2p+beta2,0)+imI2(0.25,sm,y2m+beta2,0)+imI2(0.25,mZ2/(4*m2),y2p/(1-beta2),0)-imI2(0.25,mZ2/(4.0*m2),y2m/(1-beta2),0)-imI2(0.25,mZ2/(4.0*m2),-y2p/beta2,0)+imI2(0.25,mZ2/(4.0*m2),-y2m/beta2,0);
    //printf("result imsumi2ds   %lf \n", sumimI2);
    double imdst=-1/(16.0*m2*m2*sm*tm*R2)*(sumimI1+sumimI2);
    return imdst;
}

double ReF(double sm, double tm, double um, double m2)
{
    return ReDut(tm,um,  m2)+ReDst(sm,tm,  m2)+ReDst(sm,um,  m2);
}
double ImF(double sm, double tm, double um, double m2)
{
    //printf("dut %lf \n",pow(10,13)*ImDut(tm,um,  m2));
    //printf("dst %lf \n",pow(10,13)*ImDst(sm,tm,  m2));
    //printf("dsu %lf \n",pow(10,13)*ImDst(sm,um,  m2));
    return ImDut(tm,um,  m2)+ImDst(sm,tm,  m2)+ImDst(sm,um,  m2);
}

double ReE1(double sm,double tm, double m2)
{
    return 2*(4*m2*tm-mZ2)*ReCz(tm, m2)-16*m2*m2*sm*tm*ReDst(sm,tm, m2);//return 2*t/(t-mZ/(4*m))*ReCz(t, m)-16*m*m*s*t*ReDst(s,t, m);
}
double ImE1(double sm,double tm, double m2)
{
    return 2*(4*m2*tm-mZ2)*ImCz(tm,  m2)-16*m2*m2*sm*tm*ImDst(sm,tm, m2);
}

double ReE2(double tm, double um, double m2)
{
    return 2*(4*m2*tm-mZ2)*ReCz(tm,  m2)+2*(4*m2*um-mZ2)*ReCz(um,m2)-(16*m2*m2*tm*um-mZ2*mZ2)*ReDut(tm,um,  m2);//2*t/(t-mZ/(4*m))*ReCz(t,  m)+2*u/(u-mZ/(4*m))*ReCz(u,  m)-16*m*m*(t*u-mZ/(4*m)*mZ/(4*m))*ReDut(t,u,  m);
}
double ImE2(double tm, double um, double m2)
{
    return 2*(4*m2*tm-mZ2)*ImCz(tm,  m2)+2*(4*m2*um-mZ2)*ImCz(um, m2)-(16*m2*m2*tm*um-mZ2*mZ2)*ImDut(tm,um, m2);
}




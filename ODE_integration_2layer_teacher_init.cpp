#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdlib>
using namespace std;

double I0(double qkk,double qll,double qkl)
{
    return 2/M_PI*asin(qkl/sqrt((1+qkk)*(1+qll)));
}

double I1(double qkk,double qmm,double qkm,double qkl_r,double qlm_r)
{
    return 2/M_PI/sqrt((1+qkk)*(1+qmm)-pow(qkm,2))*(qlm_r-qkm*qkl_r/(1+qkk));
}

double I2(double qkk,double qll,double qmm,double qnn,double qkl,double qkm,double qkn,double qlm,double qln,double qmn)
{
    double L4(0),L0(0),L1(0),L2(0);
    L4=(1+qkk)*(1+qll)-pow(qkl,2);
    L0=L4*qmn-qlm*qln*(1+qkk)-qkm*qkn*(1+qll)+qkl*qkm*qln+qkl*qkn*qlm;
    L1=L4*(1+qmm)-pow(qlm,2)*(1+qkk)-pow(qkm,2)*(1+qll)+2*qkl*qkm*qlm;
    L2=L4*(1+qnn)-pow(qln,2)*(1+qkk)-pow(qkn,2)*(1+qll)+2*qkl*qkn*qln;
    return 4/pow(M_PI,2)/sqrt(L4)*asin(L0/sqrt(L1*L2));
}

double error(double q[],double r[],double v[],double t[],double h[],int K,int A)
{
    double s1(0),s2(0),s3(0);
    for(int k=0; k<K; k++)
    {
        s1+=pow(v[k],2)*asin(q[(k+1)*(k+2)/2-1]/(1+q[(k+1)*(k+2)/2-1]));
        for (int l=0; l<k; l++)
        {
            s1+=2*v[k]*v[l]*asin(q[k*(k+1)/2+l]/sqrt((1+q[(k+1)*(k+2)/2-1])*(1+q[(l+1)*(l+2)/2-1])));
        }
        for (int a=0;a<A;a++)
        {
            s2+=h[a]*v[k]*asin(r[k*A+a]/sqrt((1+q[(k+1)*(k+2)/2-1])*(1+t[(a+1)*(a+2)/2-1])));
        }
    }
    for(int a=0; a<A; a++)
    {
        s3+=pow(h[a],2)*asin(t[(a+1)*(a+2)/2-1]/(1+t[(a+1)*(a+2)/2-1]));
        for (int b=0; b<a; b++)
        {
            s1+=2*h[a]*h[b]*asin(t[a*(a+1)/2+b]/sqrt((1+t[(a+1)*(a+2)/2-1])*(1+t[(b+1)*(b+2)/2-1])));
        }
    }
    return (s1-2*s2+s3)/M_PI;
}

void initQ(double q[],int K,double diagterm)
{
    for (int k=0; k<K; k++)
    {
        q[(k+1)*(k+2)/2-1]=diagterm;
        for (int l=0; l<k; l++)
        {
            q[k*(k+1)/2+l]=0;
        }
    }
}

void initR(double r[],int K,int A)
{
    for (int k=0; k<K*A; k++)
    {
        r[k]=0;
    }
}

void initV(double v[],int K,double sigma0)
{
    switch(K)
    {
    case 1:
        v[0]=0.79788456;
        break;
    case 2:
        v[0]=1.12837917;
        v[1]=0.46738995;
        break;
    case 3:
        v[0]=1.32638676;
        v[1]=0.73236399;
        v[2]=0.33490294;
        break;
    case 5:
        v[0]=1.56983372;
        v[1]=1.04430504;
        v[2]=0.71195013;
        v[3]=0.44764142;
        v[4]=0.2156925;
        break;
    case 10:
        v[0]=1.88071569;
        v[1]=1.42436122;
        v[2]=1.15086379;
        v[3]=0.94411641;
        v[4]=0.77198173;
        v[5]=0.62074666;
        v[6]=0.48316487;
        v[7]=0.35485397;
        v[8]=0.23288872;
        v[9]=0.11515255;
        break;
    default:
        cout << "Initialisation of v for this K is unknown" << endl;
        break;
    }
    if (sigma0!=1)
    {
        for (int k=0; k<K; k++)
        {
            v[k]=v[k]*sigma0;
        }
    }
}

void perm(double w[],double v[],int K)
{
    for (int k=0;k<K-1;k++)
    {
        w[k]=v[k+1];
    }
    w[K-1]=v[0];
}

void compute_f(double fqr[],double fqi[],double fr[],double fv[],double q[],double qr[],double qi[],double r[],double v[],int K,double t[],double h[],int A,double rho,double aw,double av)
{
    double c1(0),c2(0),c3(0);
    for (int k=0; k<K; k++)
    {
        fv[k]=0;
        for (int a=0;a<A;a++)
        {
            fv[k]+=av*h[a]*I0(q[(k+1)*(k+2)/2-1],t[(a+1)*(a+2)/2-1],r[k*A+a]);
        }
        for (int l=0; l<=k; l++)
        {
            fv[k]-=av*v[l]*I0(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[k*(k+1)/2+l]);
        }
        for (int l=k+1; l<K; l++)
        {
            fv[k]-=av*v[l]*I0(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[l*(l+1)/2+k]);
        }
        for (int a=0;a<A;a++)
        {
            fr[k*A+a]=0;
            for (int b=0;b<=a;b++)
            {
                fr[k*A+a]+=aw*v[k]*h[b]*I1(q[(k+1)*(k+2)/2-1],t[(b+1)*(b+2)/2-1],r[k*A+b],r[k*A+a],t[a*(a+1)/2+b]);
            }
            for (int b=a+1;b<A;b++)
            {
                fr[k*A+a]+=aw*v[k]*h[b]*I1(q[(k+1)*(k+2)/2-1],t[(b+1)*(b+2)/2-1],r[k*A+b],r[k*A+a],t[b*(b+1)/2+a]);
            }
            for (int l=0; l<=k; l++)
            {
                fr[k*A+a]-=aw*v[k]*v[l]*I1(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[k*(k+1)/2+l],r[k*A+a],r[l*A+a]);
            }
            for (int l=k+1; l<K; l++)
            {
                fr[k*A+a]-=aw*v[k]*v[l]*I1(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[l*(l+1)/2+k],r[k*A+a],r[l*A+a]);
            }
        }
        for (int l=0; l<=k; l++)
        {
            fqr[k*(k+1)/2+l]=0;
            fqi[k*(k+1)/2+l]=0;
            c1=0;
            for (int a=0;a<A;a++)
            {
                fqr[k*(k+1)/2+l]+=aw*v[k]*h[a]*I1(q[(k+1)*(k+2)/2-1],t[(a+1)*(a+2)/2-1],r[k*A+a],qr[k*(k+1)/2+l],r[l*A+a])+aw*v[l]*h[a]*I1(q[(l+1)*(l+2)/2-1],t[(a+1)*(a+2)/2-1],r[l*A+a],qr[k*(k+1)/2+l],r[k*A+a]);
                fqi[k*(k+1)/2+l]+=aw*v[k]*h[a]*I1(q[(k+1)*(k+2)/2-1],t[(a+1)*(a+2)/2-1],r[k*A+a],qi[k*(k+1)/2+l],0)+aw*v[l]*h[a]*I1(q[(l+1)*(l+2)/2-1],t[(a+1)*(a+2)/2-1],r[l*A+a],qi[k*(k+1)/2+l],0);
                for (int b=0;b<=a;b++)
                {
                    c1+=pow(aw,2)*v[k]*v[l]*h[a]*h[b]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],t[(a+1)*(a+2)/2-1],t[(b+1)*(b+2)/2-1],q[k*(k+1)/2+l],r[k*A+a],r[k*A+b],r[l*A+a],r[l*A+b],t[a*(a+1)/2+b]);
                }
                for (int b=a+1;b<A;b++)
                {
                    c1+=pow(aw,2)*v[k]*v[l]*h[a]*h[b]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],t[(a+1)*(a+2)/2-1],t[(b+1)*(b+2)/2-1],q[k*(k+1)/2+l],r[k*A+a],r[k*A+b],r[l*A+a],r[l*A+b],t[b*(b+1)/2+a]);
                }
            }
            fqr[k*(k+1)/2+l]+=rho*c1;
            fqi[k*(k+1)/2+l]+=(1-rho)*c1;
            for (int m=0; m<=l; m++)
            {
                c2=0;
                for (int a=0;a<A;a++)
                {
                    c2+=2*pow(aw,2)*v[k]*v[l]*v[m]*h[a]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],t[(a+1)*(a+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],r[k*A+a],q[l*(l+1)/2+m],r[l*A+a],r[m*A+a]);
                }
                fqr[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qr[k*(k+1)/2+l],qr[l*(l+1)/2+m])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[l*(l+1)/2+m],qr[k*(k+1)/2+l],qr[k*(k+1)/2+m])+rho*c2;
                fqi[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qi[k*(k+1)/2+l],qi[l*(l+1)/2+m])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[l*(l+1)/2+m],qi[k*(k+1)/2+l],qi[k*(k+1)/2+m])+(1-rho)*c2;
                c3=0;
                for (int n=0; n<=m; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[l*(l+1)/2+m],q[l*(l+1)/2+n],q[m*(m+1)/2+n]);
                }
                for (int n=m+1; n<=l; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[l*(l+1)/2+m],q[l*(l+1)/2+n],q[n*(n+1)/2+m]);
                }
                for (int n=l+1; n<=k; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[l*(l+1)/2+m],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                }
                for (int n=k+1; n<K; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[n*(n+1)/2+k],q[l*(l+1)/2+m],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                }
                fqr[k*(k+1)/2+l]+=rho*c3;
                fqi[k*(k+1)/2+l]+=(1-rho)*c3;
            }
            for (int m=l+1; m<=k; m++)
            {
                c2=0;
                for (int a=0;a<A;a++)
                {
                    c2+=2*pow(aw,2)*v[k]*v[l]*v[m]*h[a]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],t[(a+1)*(a+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],r[k*A+a],q[m*(m+1)/2+l],r[l*A+a],r[m*A+a]);
                }
                fqr[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qr[k*(k+1)/2+l],qr[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qr[k*(k+1)/2+l],qr[k*(k+1)/2+m])+rho*c2;
                fqi[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qi[k*(k+1)/2+l],qi[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qi[k*(k+1)/2+l],qi[k*(k+1)/2+m])+(1-rho)*c2;
                c3=0;
                for (int n=0; n<=l; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[l*(l+1)/2+n],q[m*(m+1)/2+n]);
                }
                for (int n=l+1; n<=m; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[m*(m+1)/2+n]);
                }
                for (int n=m+1; n<=k; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                }
                for (int n=k+1; n<K; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[n*(n+1)/2+k],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                }
                fqr[k*(k+1)/2+l]+=rho*c3;
                fqi[k*(k+1)/2+l]+=(1-rho)*c3;
            }
            for (int m=k+1; m<K; m++)
            {
                c2=0;
                for (int a=0;a<A;a++)
                {
                    c2+=2*pow(aw,2)*v[k]*v[l]*v[m]*h[a]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],t[(a+1)*(a+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],r[k*A+a],q[m*(m+1)/2+l],r[l*A+a],r[m*A+a]);
                }
                fqr[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+k],qr[k*(k+1)/2+l],qr[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qr[k*(k+1)/2+l],qr[m*(m+1)/2+k])+rho*c2;
                fqi[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+k],qi[k*(k+1)/2+l],qi[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qi[k*(k+1)/2+l],qi[m*(m+1)/2+k])+(1-rho)*c2;
                c3=0;
                for (int n=0; n<=l; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[l*(l+1)/2+n],q[m*(m+1)/2+n]);
                }
                for (int n=l+1; n<=k; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[m*(m+1)/2+n]);
                }
                for (int n=k+1; n<=m; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[n*(n+1)/2+k],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[m*(m+1)/2+n]);
                }
                for (int n=m+1; n<K; n++)
                {
                    c3+=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[n*(n+1)/2+k],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                }
                fqr[k*(k+1)/2+l]+=rho*c3;
                fqi[k*(k+1)/2+l]+=(1-rho)*c3;
            }
        }
    }
}

void AdamsBashforth(double qr[],double qi[],double r[],double v[],double t[],double h[],int K,int A,double rho,double aw,double av,int Tf,double dt,double dtpoints,string errFileName,string qrFileName,string qiFileName,string rFileName,string vFileName)
{
    int nsteps(Tf/dt),skipsteps(dtpoints/dt);
    double q[K*(K+1)/2],fqrprev[K*(K+1)/2],fqiprev[K*(K+1)/2],frprev[K*A],fvprev[K],fqrnew[K*(K+1)/2],fqinew[K*(K+1)/2],frnew[K*A],fvnew[K];
    ofstream errFile(errFileName.c_str());
    ofstream qrFile(qrFileName.c_str());
    ofstream qiFile(qiFileName.c_str());
    ofstream rFile(rFileName.c_str());
    ofstream vFile(vFileName.c_str());
    if(errFile && qrFile && qiFile && rFile && vFile)
    {
        for (int k=0; k<K*(K+1)/2; k++)
        {
            q[k]=qr[k]+qi[k];
        }
        errFile << error(q,r,v,t,h,K,A) << endl;
        for (int k=0; k<K*(K+1)/2; k++)
        {
            qrFile << qr[k] << " ";
            qiFile << qi[k] << " ";
        }
        qrFile << endl;
        qiFile << endl;
        for (int k=0; k<K*A; k++)
        {
            rFile << r[k] << " ";
        }
        rFile << endl;
        for (int k=0; k<K; k++)
        {
            vFile << v[k] << " ";
        }
        vFile << endl;
        compute_f(fqrprev,fqiprev,frprev,fvprev,q,qr,qi,r,v,K,t,h,A,rho,aw,av);
        for (int k=0; k<K*(K+1)/2; k++)
        {
            qr[k]+=dt*fqrprev[k];
            qi[k]+=dt*fqiprev[k];
            q[k]=qr[k]+qi[k];
        }
        for (int k=0; k<K*A; k++)
        {
            r[k]+=dt*frprev[k];
        }
        for (int k=0; k<K; k++)
        {
            v[k]+=dt*fvprev[k];
        }
        if (dtpoints==dt)
        {
            errFile << error(q,r,v,t,h,K,A) << endl;
            for (int k=0; k<K*(K+1)/2; k++)
            {
                qrFile << qr[k] << " ";
                qiFile << qi[k] << " ";
            }
            qrFile << endl;
            qiFile << endl;
            for (int k=0; k<K*A; k++)
            {
                rFile << r[k] << " ";
            }
            rFile << endl;
            for (int k=0; k<K; k++)
            {
                vFile << v[k] << " ";
            }
            vFile << endl;
        }
        for (int i=1; i<nsteps; i++)
        {
            compute_f(fqrnew,fqinew,frnew,fvnew,q,qr,qi,r,v,K,t,h,A,rho,aw,av);
            for (int k=0; k<K*(K+1)/2; k++)
            {
                qr[k]+=dt*(3*fqrnew[k]-fqrprev[k])/2;
                qi[k]+=dt*(3*fqinew[k]-fqiprev[k])/2;
                q[k]=qr[k]+qi[k];
                fqrprev[k]=fqrnew[k];
                fqiprev[k]=fqinew[k];
            }
            for (int k=0; k<K*A; k++)
            {
                r[k]+=dt*(3*frnew[k]-frprev[k])/2;
                frprev[k]=frnew[k];
            }
            for (int k=0; k<K; k++)
            {
                v[k]+=dt*(3*fvnew[k]-fvprev[k])/2;
                fvprev[k]=fvnew[k];
            }
            if ((i+1)%skipsteps==0)
            {
                cout << (i+1)/skipsteps << endl;
                errFile << error(q,r,v,t,h,K,A) << endl;
                for (int k=0; k<K*(K+1)/2; k++)
                {
                    qrFile << qr[k] << " ";
                    qiFile << qi[k] << " ";
                }
                qrFile << endl;
                qiFile << endl;
                for (int k=0; k<K*A; k++)
                {
                    rFile << r[k] << " ";
                }
                rFile << endl;
                for (int k=0; k<K; k++)
                {
                    vFile << v[k] << " ";
                }
                vFile << endl;
            }
        }
    }
    else
    {
        cout << "ERROR: cannot open file" << endl;
    }
}

void AdamsBashforth_t1(const int Tpre[3],double qrt[][3],double qit[][3],double vt[][3],double qr[],double qi[],double r[],double v[],double t[],double h[],int K,int A,double rho,double aw,double av,int Tf,double dt,double dtpoints,string errFileName,string qrFileName,string qiFileName,string rFileName,string vFileName)
{
    int nsteps(Tf/dt),skipsteps(dtpoints/dt),c(0);
    double q[K*(K+1)/2],fqrprev[K*(K+1)/2],fqiprev[K*(K+1)/2],frprev[K*A],fvprev[K],fqrnew[K*(K+1)/2],fqinew[K*(K+1)/2],frnew[K*A],fvnew[K];
    ofstream errFile(errFileName.c_str());
    ofstream qrFile(qrFileName.c_str());
    ofstream qiFile(qiFileName.c_str());
    ofstream rFile(rFileName.c_str());
    ofstream vFile(vFileName.c_str());
    if(errFile && qrFile && qiFile && rFile && vFile)
    {
        for (int k=0; k<K*(K+1)/2; k++)
        {
            q[k]=qr[k]+qi[k];
        }
        errFile << error(q,r,v,t,h,K,A) << endl;
        for (int k=0; k<K*(K+1)/2; k++)
        {
            qrFile << qr[k] << " ";
            qiFile << qi[k] << " ";
        }
        qrFile << endl;
        qiFile << endl;
        for (int k=0; k<K*A; k++)
        {
            rFile << r[k] << " ";
        }
        rFile << endl;
        for (int k=0; k<K; k++)
        {
            vFile << v[k] << " ";
        }
        vFile << endl;
        compute_f(fqrprev,fqiprev,frprev,fvprev,q,qr,qi,r,v,K,t,h,A,rho,aw,av);
        for (int k=0; k<K*(K+1)/2; k++)
        {
            qr[k]+=dt*fqrprev[k];
            qi[k]+=dt*fqiprev[k];
            q[k]=qr[k]+qi[k];
        }
        for (int k=0; k<K*A; k++)
        {
            r[k]+=dt*frprev[k];
        }
        for (int k=0; k<K; k++)
        {
            v[k]+=dt*fvprev[k];
        }
        if (dtpoints==dt)
        {
            errFile << error(q,r,v,t,h,K,A) << endl;
            for (int k=0; k<K*(K+1)/2; k++)
            {
                qrFile << qr[k] << " ";
                qiFile << qi[k] << " ";
            }
            qrFile << endl;
            qiFile << endl;
            for (int k=0; k<K*A; k++)
            {
                rFile << r[k] << " ";
            }
            rFile << endl;
            for (int k=0; k<K; k++)
            {
                vFile << v[k] << " ";
            }
            vFile << endl;
        }
        for (int i=1; i<nsteps; i++)
        {
            compute_f(fqrnew,fqinew,frnew,fvnew,q,qr,qi,r,v,K,t,h,A,rho,aw,av);
            for (int k=0; k<K*(K+1)/2; k++)
            {
                qr[k]+=dt*(3*fqrnew[k]-fqrprev[k])/2;
                qi[k]+=dt*(3*fqinew[k]-fqiprev[k])/2;
                q[k]=qr[k]+qi[k];
                fqrprev[k]=fqrnew[k];
                fqiprev[k]=fqinew[k];
            }
            for (int k=0; k<K*A; k++)
            {
                r[k]+=dt*(3*frnew[k]-frprev[k])/2;
                frprev[k]=frnew[k];
            }
            for (int k=0; k<K; k++)
            {
                v[k]+=dt*(3*fvnew[k]-fvprev[k])/2;
                fvprev[k]=fvnew[k];
            }
            if ((i+1)%skipsteps==0)
            {
                cout << (i+1)/skipsteps << endl;
                errFile << error(q,r,v,t,h,K,A) << endl;
                for (int k=0; k<K*(K+1)/2; k++)
                {
                    qrFile << qr[k] << " ";
                    qiFile << qi[k] << " ";
                }
                qrFile << endl;
                qiFile << endl;
                for (int k=0; k<K*A; k++)
                {
                    rFile << r[k] << " ";
                }
                rFile << endl;
                for (int k=0; k<K; k++)
                {
                    vFile << v[k] << " ";
                }
                vFile << endl;
                if ((i+1)==Tpre[c]/dt)
                {
                    for (int k=0; k<K*(K+1)/2; k++)
                    {
                        qrt[k][c]=qr[k];
                        qit[k][c]=qi[k];
                    }
                    for (int k=0; k<K; k++)
                    {
                        vt[k][c]=v[k];
                    }
                    c+=1;
                }
            }
        }
    }
    else
    {
        cout << "ERROR: cannot open file" << endl;
    }
}

void main2(int K,int A,double rho)
{
    int Tf(10000);
    //int K(3),A(3);
    double sigma0(1),aw(1),av(1),dt(0.1),dtpoints(10);
    double qr[K*(K+1)/2];//={0.50039316,  0.01088364,  0.49960063, -0.03300625,  0.02103842, 0.45028161};//={0.46263299};//={0.10167589, 0.00325954, 0.11671978, 0.00081399, 0.00292362, 0.11548411};
    double qi[K*(K+1)/2];//={0.51315949, -0.00210669,  0.48445879, -0.013012  ,  0.00280359, 0.47636638};//={0.52022594};//={8.86334676e-01, -4.33296290e-04,  8.61724268e-01, -7.38692673e-03, -1.50315355e-02,  9.24392065e-01};
    double r[K*A];//={0.01555814,  0.00771486, 0.02811012,  0.03416652,-0.02251075, -0.05969664};//={-0.0005068 , -0.07583292,  0.04702699};//={-0.000060384,  0.0001706711, -0.00070461, -0.000189088, -0.000519145,  0.000919399, 0.0000207166, -0.0001229755, -0.0000724259};
    double v[K];//={-0.72379809};//={0.19253684, 0.77432526, 0.71316453};//={-2.08101922, -0.3100416 ,  0.79853657};;//={-0.20267418,  1.10420871,  0.77422588};//{0.44808509, 0.56382914, 0.42514543};
    double t[A*(A+1)/2];//={0.52173389, 0.00647602, 0.45166125};//={0.4453226 ,  0.05280237,  0.49768827, -0.08251821, -0.013367  ,0.54396388};//={0.0947562 , -0.00681219,  0.09496094, -0.00891996, -0.01067338, 0.0899326};
    double h1[A],h2[A];//={-1.17965869,  0.11397837,  0.42285759}//={0.88276946, 1.01772186, 1.30290494};//h0[A],h1[A]={0.67829355, 1.315615, 0.300041},h2[A]={0.88725544, 1.19158723, 0.51130675};
    double qrt[K*(K+1)/2],qit[K*(K+1)/2];
    string errFileName("results/err_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+".txt");
    string qrFileName("results/qr_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+".txt");
    string qiFileName("results/qi_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+".txt");
    string rFileName("results/r_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+".txt");
    string vFileName("results/v_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+".txt");
    initQ(qr,K,rho*sigma0);
    initQ(qi,K,(1-rho)*sigma0);
    initQ(t,A,rho);
    initR(r,K,A);
    initV(v,K,sigma0);
    initV(h1,A,1);
    initV(h2,A,1);
    AdamsBashforth(qr,qi,r,v,t,h1,K,A,rho,aw,av,Tf,dt,dtpoints,errFileName,qrFileName,qiFileName,rFileName,vFileName);
    for (int k=0;k<K*(K+1)/2;k++)
    {
        qrt[k]=qr[k];
        qit[k]=qi[k];
    }
    errFileName="results/err_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_twv.txt";
    qrFileName="results/qr_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_twv.txt";
    qiFileName="results/qi_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_twv.txt";
    rFileName="results/r_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_twv.txt";
    vFileName="results/v_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_twv.txt";
    initR(r,K,A);
    AdamsBashforth(qr,qi,r,v,t,h2,K,A,rho,aw,av,Tf,dt,dtpoints,errFileName,qrFileName,qiFileName,rFileName,vFileName);
    errFileName="results/err_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_tw.txt";
    qrFileName="results/qr_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_tw.txt";
    qiFileName="results/qi_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_tw.txt";
    rFileName="results/r_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_tw.txt";
    vFileName="results/v_A"+to_string(A)+"_K"+to_string(K)+"_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_sigma0"+to_string(sigma0)+"_Tf"+to_string(Tf)+"_tw.txt";
    for (int k=0;k<K*(K+1)/2;k++)
    {
        qr[k]=qrt[k];
        qi[k]=qit[k];
    }
    initR(r,K,A);
    initV(v,K,sigma0);
    AdamsBashforth(qr,qi,r,v,t,h2,K,A,rho,aw,av,Tf,dt,dtpoints,errFileName,qrFileName,qiFileName,rFileName,vFileName);
}

int main(int argc,char *argv[])
{
    int A(3);
    if(argc==1)
    {
        A=3;
    }
    else
    {
        A = atoi(argv[1]);
    }
    //double aw[7]={0.1,0.2,0.5,1,2,5};
    //double rho[5]={0.1,0.25,0.5,0.75,1};
    //double av[7]={0.02,0.05,0.1,0.2,0.5,1,2};
    int K[5]={1,2,3,5,10};
    for (int i=0;i<1;i++)
    {
        main2(K[i],A,0.1);
        main2(K[i],A,0.5);
    }
    //double v[10][3]={{-1.53221889,  1.13410588,  0.3884811},{-0.63593858, -0.57983309,  0.8632137},
                    //{-1.07025762,  1.8069262 ,  0.65021007},{0.28890048, 1.5267316 , 1.52741228},
                    //{-0.58388372,  0.47888286,  0.09093863},{0.2358605 , -0.20459371,  0.76424595},
                    //{0.12548513, 1.36244078, 2.06754005},{0.62415929,  0.33575737, -1.15723042},
                    //{0.88276946, 1.01772186, 1.30290494},{-0.69126245, -1.63396207,  0.37421204}};
    //double h2[10][3]={{1.58866146, 0.43675005, 0.58254912},{-0.07455767,  0.27028633, -0.22681131},
                    //{0.55474851,  1.88143563, -0.24043112},{0.58255734, -0.90889123,  1.25104802},
                    //{-0.2241376 ,  0.64409892,  0.23844966},{0.75296138, 0.70715091, 0.44594371},
                    //{-0.73214019,  0.21978161, -0.89188138},{2.65527279,  0.71451529, -1.77304739},
                    //{0.19253684, 0.77432526, 0.71316453},{-0.33309497,  1.32210356,  1.13402017}};
    //for (int k=0;k<10;k++)
    //{
        //main2(v[k],0.1,k);
        //main2(v[k],0.5,k);
    //}
    return 0;
}

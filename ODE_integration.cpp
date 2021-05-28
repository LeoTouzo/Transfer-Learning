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

double error(double q[],double r[],double v[],double t,int K)
{
    double s1(0),s2(0);
    for(int k=0; k<K; k++)
    {
        s1+=pow(v[k],2)*asin(q[(k+1)*(k+2)/2-1]/(1+q[(k+1)*(k+2)/2-1]));
        s2+=v[k]*asin(r[k]/sqrt((1+q[(k+1)*(k+2)/2-1])*(1+t)));
        for (int l=0; l<k; l++)
        {
            s1+=2*v[k]*v[l]*asin(q[k*(k+1)/2+l]/sqrt((1+q[(k+1)*(k+2)/2-1])*(1+q[(l+1)*(l+2)/2-1])));
        }
    }
    return (asin(t/(1+t))+s1-2*s2)/M_PI;
}

void initQ(double qr[],double qi[],int K,double rho,double sigma0)
{
    for (int k=0; k<K; k++)
    {
        qr[(k+1)*(k+2)/2-1]=rho*sigma0;
        qi[(k+1)*(k+2)/2-1]=(1-rho)*sigma0;
        for (int l=0; l<k; l++)
        {
            qr[k*(k+1)/2+l]=0;
            qi[k*(k+1)/2+l]=0;
        }
    }
}

void initR(double r[],int K,double rho,double sigma0)
{
    for (int k=0; k<K; k++)
    {
        r[k]=0;
    }
}

void initV(double v[],int K,double sigma0)
{
    switch(K)
    {
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

void compute_f(double fqr[],double fqi[],double fr[],double fv[],double q[],double qr[],double qi[],double r[],double v[],int K,double t,double rho,double aw,double av)
{
    double a(0),b(0),c(0);
    for (int k=0; k<K; k++)
    {
        fv[k]=av*I0(q[(k+1)*(k+2)/2-1],t,r[k]);
        for (int l=0; l<=k; l++)
        {
            fv[k]-=av*v[l]*I0(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[k*(k+1)/2+l]);
        }
        for (int l=k+1; l<K; l++)
        {
            fv[k]-=av*v[l]*I0(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[l*(l+1)/2+k]);
        }
        fr[k]=aw*v[k]*I1(q[(k+1)*(k+2)/2-1],t,r[k],r[k],t);
        for (int l=0; l<=k; l++)
        {
            fr[k]-=aw*v[k]*v[l]*I1(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[k*(k+1)/2+l],r[k],r[l]);
        }
        for (int l=k+1; l<K; l++)
        {
            fr[k]-=aw*v[k]*v[l]*I1(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[l*(l+1)/2+k],r[k],r[l]);
        }
        for (int l=0; l<=k; l++)
        {
            a=pow(aw,2)*v[k]*v[l]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],t,t,q[k*(k+1)/2+l],r[k],r[k],r[l],r[l],t);
            fqr[k*(k+1)/2+l]=aw*v[k]*I1(q[(k+1)*(k+2)/2-1],t,r[k],qr[k*(k+1)/2+l],r[l])+aw*v[l]*I1(q[(l+1)*(l+2)/2-1],t,r[l],qr[k*(k+1)/2+l],r[k])+rho*a;
            fqi[k*(k+1)/2+l]=aw*v[k]*I1(q[(k+1)*(k+2)/2-1],t,r[k],qi[k*(k+1)/2+l],0)+aw*v[l]*I1(q[(l+1)*(l+2)/2-1],t,r[l],qi[k*(k+1)/2+l],0)+(1-rho)*a;
            for (int m=0; m<=l; m++)
            {
                b=2*pow(aw,2)*v[k]*v[l]*v[m]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],t,q[k*(k+1)/2+l],q[k*(k+1)/2+m],r[k],q[l*(l+1)/2+m],r[l],r[m]);
                fqr[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qr[k*(k+1)/2+l],qr[l*(l+1)/2+m])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[l*(l+1)/2+m],qr[k*(k+1)/2+l],qr[k*(k+1)/2+m])+rho*b;
                fqi[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qi[k*(k+1)/2+l],qi[l*(l+1)/2+m])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[l*(l+1)/2+m],qi[k*(k+1)/2+l],qi[k*(k+1)/2+m])+(1-rho)*b;
                for (int n=0; n<=m; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[l*(l+1)/2+m],q[l*(l+1)/2+n],q[m*(m+1)/2+n]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=m+1; n<=l; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[l*(l+1)/2+m],q[l*(l+1)/2+n],q[n*(n+1)/2+m]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=l+1; n<=k; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[l*(l+1)/2+m],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=k+1; n<K; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[n*(n+1)/2+k],q[l*(l+1)/2+m],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
            }
            for (int m=l+1; m<=k; m++)
            {
                b=2*pow(aw,2)*v[k]*v[l]*v[m]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],t,q[k*(k+1)/2+l],q[k*(k+1)/2+m],r[k],q[m*(m+1)/2+l],r[l],r[m]);
                fqr[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qr[k*(k+1)/2+l],qr[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qr[k*(k+1)/2+l],qr[k*(k+1)/2+m])+rho*b;
                fqi[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[k*(k+1)/2+m],qi[k*(k+1)/2+l],qi[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qi[k*(k+1)/2+l],qi[k*(k+1)/2+m])+(1-rho)*b;
                for (int n=0; n<=l; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[l*(l+1)/2+n],q[m*(m+1)/2+n]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=l+1; n<=m; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[m*(m+1)/2+n]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=m+1; n<=k; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=k+1; n<K; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[k*(k+1)/2+m],q[n*(n+1)/2+k],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
            }
            for (int m=k+1; m<K; m++)
            {
                b=2*pow(aw,2)*v[k]*v[l]*v[m]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],t,q[k*(k+1)/2+l],q[m*(m+1)/2+k],r[k],q[m*(m+1)/2+l],r[l],r[m]);
                fqr[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+k],qr[k*(k+1)/2+l],qr[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qr[k*(k+1)/2+l],qr[m*(m+1)/2+k])+rho*b;
                fqi[k*(k+1)/2+l]-=aw*v[k]*v[m]*I1(q[(k+1)*(k+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+k],qi[k*(k+1)/2+l],qi[m*(m+1)/2+l])+aw*v[l]*v[m]*I1(q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[m*(m+1)/2+l],qi[k*(k+1)/2+l],qi[m*(m+1)/2+k])+(1-rho)*b;
                for (int n=0; n<=l; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[l*(l+1)/2+n],q[m*(m+1)/2+n]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=l+1; n<=k; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[k*(k+1)/2+n],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[m*(m+1)/2+n]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=k+1; n<=m; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[n*(n+1)/2+k],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[m*(m+1)/2+n]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
                for (int n=m+1; n<K; n++)
                {
                    c=pow(aw,2)*v[k]*v[l]*v[m]*v[n]*I2(q[(k+1)*(k+2)/2-1],q[(l+1)*(l+2)/2-1],q[(m+1)*(m+2)/2-1],q[(n+1)*(n+2)/2-1],q[k*(k+1)/2+l],q[m*(m+1)/2+k],q[n*(n+1)/2+k],q[m*(m+1)/2+l],q[n*(n+1)/2+l],q[n*(n+1)/2+m]);
                    fqr[k*(k+1)/2+l]+=rho*c;
                    fqi[k*(k+1)/2+l]+=(1-rho)*c;
                }
            }
        }
    }
}

void AdamsBashforth(double qr[],double qi[],double r[],double v[],double t,int K,double rho,double aw,double av,int Tf,double dt,double dtpoints,string errFileName,string qrFileName,string qiFileName,string rFileName,string vFileName)
{
    int nsteps(Tf/dt),skipsteps(dtpoints/dt);
    double q[K*(K+1)/2],fqrprev[K*(K+1)/2],fqiprev[K*(K+1)/2],frprev[K],fvprev[K],fqrnew[K*(K+1)/2],fqinew[K*(K+1)/2],frnew[K],fvnew[K];
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
        errFile << error(q,r,v,t,K) << endl;
        for (int k=0; k<K*(K+1)/2; k++)
        {
            qrFile << qr[k] << " ";
            qiFile << qi[k] << " ";
        }
        qrFile << endl;
        qiFile << endl;
        for (int k=0; k<K; k++)
        {
            rFile << r[k] << " ";
            vFile << v[k] << " ";
        }
        rFile << endl;
        vFile << endl;
        compute_f(fqrprev,fqiprev,frprev,fvprev,q,qr,qi,r,v,K,t,rho,aw,av);
        for (int k=0; k<K*(K+1)/2; k++)
        {
            qr[k]+=dt*fqrprev[k];
            qi[k]+=dt*fqiprev[k];
            q[k]=qr[k]+qi[k];
        }
        for (int k=0; k<K; k++)
        {
            r[k]+=dt*frprev[k];
            v[k]+=dt*fvprev[k];
        }
        if (dtpoints==dt)
        {
            errFile << error(q,r,v,t,K) << endl;
            for (int k=0; k<K*(K+1)/2; k++)
            {
                qrFile << qr[k] << " ";
                qiFile << qi[k] << " ";
            }
            qrFile << endl;
            qiFile << endl;
            for (int k=0; k<K; k++)
            {
                rFile << r[k] << " ";
                vFile << v[k] << " ";
            }
            rFile << endl;
            vFile << endl;
        }
        for (int i=1; i<nsteps; i++)
        {
            compute_f(fqrnew,fqinew,frnew,fvnew,q,qr,qi,r,v,K,t,rho,aw,av);
            for (int k=0; k<K*(K+1)/2; k++)
            {
                qr[k]+=dt*(3*fqrnew[k]-fqrprev[k])/2;
                qi[k]+=dt*(3*fqinew[k]-fqiprev[k])/2;
                q[k]=qr[k]+qi[k];
                fqrprev[k]=fqrnew[k];
                fqiprev[k]=fqinew[k];
            }
            for (int k=0; k<K; k++)
            {
                r[k]+=dt*(3*frnew[k]-frprev[k])/2;
                v[k]+=dt*(3*fvnew[k]-fvprev[k])/2;
                frprev[k]=frnew[k];
                fvprev[k]=fvnew[k];
            }
            if ((i+1)%skipsteps==0)
            {
                cout << (i+1)/skipsteps << endl;
                errFile << error(q,r,v,t,K) << endl;
                for (int k=0; k<K*(K+1)/2; k++)
                {
                    qrFile << qr[k] << " ";
                    qiFile << qi[k] << " ";
                }
                qrFile << endl;
                qiFile << endl;
                for (int k=0; k<K; k++)
                {
                    rFile << r[k] << " ";
                    vFile << v[k] << " ";
                }
                rFile << endl;
                vFile << endl;
            }
        }
    }
    else
    {
        cout << "ERROR: cannot open file" << endl;
    }
}

void main2(double aw,double av)
{
    int const K(10),Tf(10000);
    int npoints(0);
    double rho(0.5),sigma0(1),dt(0.1),dtpoints(10);
    double qr[K*(K+1)/2],qi[K*(K+1)/2],q[K*(K+1)/2],r[K],v[K],t(0);//,qrt[K*(K+1)/2],qit[K*(K+1)/2];
    string errFileName("results/err_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_Tf"+to_string(Tf)+".txt"),qrFileName("results/qr_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_Tf"+to_string(Tf)+".txt"),qiFileName("results/qi_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_Tf"+to_string(Tf)+".txt"),rFileName("results/r_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_Tf"+to_string(Tf)+".txt"),vFileName("results/v_rho"+to_string(rho)+"_aw"+to_string(aw)+"_av"+to_string(av)+"_Tf"+to_string(Tf)+".txt");
    npoints=Tf/dtpoints+1;
    t=rho;
    initQ(qr,qi,K,rho,sigma0);
    initR(r,K,rho,sigma0);
    initV(v,K,sigma0);
    AdamsBashforth(qr,qi,r,v,t,K,rho,aw,av,Tf,dt,dtpoints,errFileName,qrFileName,qiFileName,rFileName,vFileName);
    //for (int k=0; k<K*(K+1)/2; k++)
    //{
        //qrt[k]=qr[k];
        //qit[k]=qi[k];
    //}
    //initR(r,K,rho,sigma0);
    //errFileName="results/freezeW_err_rho01_av005_Tf"+to_string(Tf)+"_twv.txt";
    //qrFileName="results/freezeW_qr_rho01_av005_Tf"+to_string(Tf)+"_twv.txt";
    //qiFileName="results/freezeW_qi_rho01_av005_Tf"+to_string(Tf)+"_twv.txt";
    //rFileName="results/freezeW_r_rho01_av005_Tf"+to_string(Tf)+"_twv.txt";
    //vFileName="results/freezeW_v_rho01_av005_Tf"+to_string(Tf)+"_twv.txt";
    //AdamsBashforth(qr,qi,r,v,t,K,rho,aw,av,Tf,dt,dtpoints,errFileName,qrFileName,qiFileName,rFileName,vFileName);
    //for (int k=0; k<K*(K+1)/2; k++)
    //{
        //qr[k]=qrt[k];
        //qi[k]=qit[k];
    //}
    //initR(r,K,rho,sigma0);
    //initV(v,K,sigma0);
    //errFileName="results/freezeW_err_rho01_av005_Tf"+to_string(Tf)+"_tw.txt";
    //qrFileName="results/freezeW_qr_rho01_av005_Tf"+to_string(Tf)+"_tw.txt";
    //qiFileName="results/freezeW_qi_rho01_av005_Tf"+to_string(Tf)+"_tw.txt";
    //rFileName="results/freezeW_r_rho01_av005_Tf"+to_string(Tf)+"_tw.txt";
    //vFileName="results/freezeW_v_rho01_av005_Tf"+to_string(Tf)+"_tw.txt";
    //AdamsBashforth(qr,qi,r,v,t,K,rho,aw,av,Tf,dt,dtpoints,errFileName,qrFileName,qiFileName,rFileName,vFileName);
    //return 0;
}

int main(int argc,char *argv[])
{
    double aw(0.5),av[7]={0.02,0.05,0.1,0.2,0.5,1,2};
    if(argc==1)
    {
        aw = 0.5;
    }
    else
    {
        aw = atof(argv[1]);
    }
    for (int k=0;k<7;k++)
    {
       main2(aw,av[k]);
    }
    return 0;
}

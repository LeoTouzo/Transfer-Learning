# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 10:43:39 2021

@author: leoto
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.special import erf
from scipy.special import binom
from scipy.integrate import quad

def I0(qkk,qll,qkl):
    return 2/np.pi*np.arcsin(qkl/np.sqrt((1+qkk)*(1+qll)))

def I1(qkk,qmm,qkm,qkl_r,qlm_r):
    return 2/np.pi/np.sqrt((1+qkk)*(1+qmm)-qkm**2)*(qlm_r-qkm*qkl_r/(1+qkk))

def I2(qkk,qll,qmm,qnn,qkl,qkm,qkn,qlm,qln,qmn):
    L4=(1+qkk)*(1+qll)-qkl**2
    L0=L4*qmn-qlm*qln*(1+qkk)-qkm*qkn*(1+qll)+qkl*qkm*qln+qkl*qkn*qlm
    L1=L4*(1+qmm)-qlm**2*(1+qkk)-qkm**2*(1+qll)+2*qkl*qkm*qlm
    L2=L4*(1+qnn)-qln**2*(1+qkk)-qkn**2*(1+qll)+2*qkl*qkn*qln
    return 4/np.pi**2/np.sqrt(L4)*np.arcsin(L0/np.sqrt(L1*L2))

def error(q,r,v,t,K):
    s1,s2=0,0
    for k in range(K):
        s1+=v[k]**2*np.arcsin(q[k,k]/(1+q[k,k]))
        s2+=v[k]*np.arcsin(r[k]/np.sqrt((1+q[k,k])*(1+t)))
        for l in range(k+1,K):
            s1+=2*v[k]*v[l]*np.arcsin(q[k,l]/np.sqrt((1+q[k,k])*(1+q[l,l])))
    return (np.arcsin(t/(1+t))+s1-2*s2)/np.pi

def error2(q,r,v,t,K):
    s1,s2=0,0
    for k in range(K):
        s1+=v[k]**2*np.arcsin(q[k,k]/(1+q[k,k]))
        s2+=v[k]*np.arcsin(np.sign(r[k])*np.minimum(np.abs(r[k]),np.sqrt(q[k,k]))/np.sqrt((1+q[k,k])*(1+t)))
        for l in range(k+1,K):
            s1+=2*v[k]*v[l]*np.arcsin(q[k,l]/np.sqrt((1+q[k,k])*(1+q[l,l])))
    return (np.arcsin(t/(1+t))+s1-2*s2)/np.pi

def to_y(qr,qi,r,v,K):
    y=np.zeros(K*(K+3))
    for k in range(K):
        for l in range (k,K):
            y[int(k*(K-(k-1)/2)+l-k)]=qr[k,l]
            y[int(K*(K+1)/2+k*(K-(k-1)/2)+l-k)]=qi[k,l]
        y[K*(K+1)+k]=r[k]
        y[K*(K+2)+k]=v[k]
    return y

def from_y(y,K):
    qr,qi,r,v=np.zeros((K,K)),np.zeros((K,K)),np.zeros(K),np.zeros(K)
    for k in range(K):
        for l in range(k,K):
            qr[k,l]=y[int(k*(K-(k-1)/2)+l-k)]
            qr[l,k]=y[int(k*(K-(k-1)/2)+l-k)]
            qi[k,l]=y[int(K*(K+1)/2+k*(K-(k-1)/2)+l-k)]
            qi[l,k]=y[int(K*(K+1)/2+k*(K-(k-1)/2)+l-k)]
        r[k]=y[K*(K+1)+k]
        v[k]=y[K*(K+2)+k]
    return qr,qi,r,v

def f(time,y,t,K,rho,aw,av):
    qr,qi,r,v=from_y(y,K)
    q=qr+qi
    fqr,fqi,fr,fv=np.zeros((K,K)),np.zeros((K,K)),np.zeros(K),np.zeros(K)
    for k in range(K):
        fv[k]=av*(I0(q[k,k],t,r[k])-np.sum(v*I0(q[k,k],np.diag(q),q[k])))
        fr[k]=aw*v[k]*(I1(q[k,k],t,r[k],r[k],t)-np.sum(v*I1(q[k,k],np.diag(q),q[k],r[k],r)))
        for l in range(k,K):
            fqr[k,l]=aw*v[k]*(I1(q[k,k],t,r[k],qr[k,l],r[l])-np.sum(v*I1(q[k,k],np.diag(q),q[k],qr[k,l],qr[l])))+aw*v[l]*(I1(q[l,l],t,r[l],qr[l,k],r[k])-np.sum(v*I1(q[l,l],np.diag(q),q[l],qr[l,k],qr[k])))+rho*aw**2*v[k]*v[l]*(I2(q[k,k],q[l,l],t,t,q[k,l],r[k],r[k],r[l],r[l],t)-2*np.sum(v*I2(q[k,k],q[l,l],np.diag(q),t,q[k,l],q[k],r[k],q[l],r[l],r))+np.sum(np.dot(np.array([v]).T,np.array([v]))*I2(q[k,k],q[l,l],np.array([np.diag(q) for i in range(K)]).T,np.array([np.diag(q) for i in range(K)]),q[k,l],np.array([q[k] for i in range(K)]).T,np.array([q[k] for i in range(K)]),np.array([q[l] for i in range(K)]).T,np.array([q[l] for i in range(K)]),q)))
            fqi[k,l]=aw*v[k]*(I1(q[k,k],t,r[k],qi[k,l],0)-np.sum(v*I1(q[k,k],np.diag(q),q[k],qi[k,l],qi[l])))+aw*v[l]*(I1(q[l,l],t,r[l],qi[l,k],0)-np.sum(v*I1(q[l,l],np.diag(q),q[l],qi[l,k],qi[k])))+(1-rho)*aw**2*v[k]*v[l]*(I2(q[k,k],q[l,l],t,t,q[k,l],r[k],r[k],r[l],r[l],t)-2*np.sum(v*I2(q[k,k],q[l,l],np.diag(q),t,q[k,l],q[k],r[k],q[l],r[l],r))+np.sum(np.dot(np.array([v]).T,np.array([v]))*I2(q[k,k],q[l,l],np.array([np.diag(q) for i in range(K)]).T,np.array([np.diag(q) for i in range(K)]),q[k,l],np.array([q[k] for i in range(K)]).T,np.array([q[k] for i in range(K)]),np.array([q[l] for i in range(K)]).T,np.array([q[l] for i in range(K)]),q)))
    return to_y(fqr,fqi,fr,fv,K)

def f2(y,t,K,rho,aw,av):
    qr,qi,r,v=from_y(y,K)
    q=qr+qi
    fqr,fqi,fr,fv=np.zeros((K,K)),np.zeros((K,K)),np.zeros(K),np.zeros(K)
    for k in range(K):
        fv[k]=av*(I0(q[k,k],t,r[k])-np.sum(v*I0(q[k,k],np.diag(q),q[k])))
        fr[k]=aw*v[k]*(I1(q[k,k],t,r[k],r[k],t)-np.sum(v*I1(q[k,k],np.diag(q),q[k],r[k],r)))
        for l in range(k,K):
            fqr[k,l]=aw*v[k]*(I1(q[k,k],t,r[k],qr[k,l],r[l])-np.sum(v*I1(q[k,k],np.diag(q),q[k],qr[k,l],qr[l])))+aw*v[l]*(I1(q[l,l],t,r[l],qr[l,k],r[k])-np.sum(v*I1(q[l,l],np.diag(q),q[l],qr[l,k],qr[k])))+rho*aw**2*v[k]*v[l]*(I2(q[k,k],q[l,l],t,t,q[k,l],r[k],r[k],r[l],r[l],t)-2*np.sum(v*I2(q[k,k],q[l,l],np.diag(q),t,q[k,l],q[k],r[k],q[l],r[l],r))+np.sum(np.dot(np.array([v]).T,np.array([v]))*I2(q[k,k],q[l,l],np.array([np.diag(q) for i in range(K)]).T,np.array([np.diag(q) for i in range(K)]),q[k,l],np.array([q[k] for i in range(K)]).T,np.array([q[k] for i in range(K)]),np.array([q[l] for i in range(K)]).T,np.array([q[l] for i in range(K)]),q)))
            fqi[k,l]=aw*v[k]*(I1(q[k,k],t,r[k],qi[k,l],0)-np.sum(v*I1(q[k,k],np.diag(q),q[k],qi[k,l],qi[l])))+aw*v[l]*(I1(q[l,l],t,r[l],qi[l,k],0)-np.sum(v*I1(q[l,l],np.diag(q),q[l],qi[l,k],qi[k])))+(1-rho)*aw**2*v[k]*v[l]*(I2(q[k,k],q[l,l],t,t,q[k,l],r[k],r[k],r[l],r[l],t)-2*np.sum(v*I2(q[k,k],q[l,l],np.diag(q),t,q[k,l],q[k],r[k],q[l],r[l],r))+np.sum(np.dot(np.array([v]).T,np.array([v]))*I2(q[k,k],q[l,l],np.array([np.diag(q) for i in range(K)]).T,np.array([np.diag(q) for i in range(K)]),q[k,l],np.array([q[k] for i in range(K)]).T,np.array([q[k] for i in range(K)]),np.array([q[l] for i in range(K)]).T,np.array([q[l] for i in range(K)]),q)))
    return to_y(fqr,fqi,fr,fv,K)

# def Euler(y0,t,K,rho,aw,av,Tf,dt,dtpoints):
#     nsteps=Tf//dt
#     npoints=Tf//dtpoints+1
#     skipsteps=dtpoints//dt
#     y=y0.copy()
#     ylist=np.zeros((len(y0),npoints))
#     ylist[:,0]=y0
#     for i in range(nsteps):
#         y+=dt*f(y,t,K,rho,aw,av)
#         if i%skipsteps==0:
#             ylist[:,(i+1)//skipsteps]=y
#     return ylist
    
def AdamsBashforth(y0,t,K,rho,aw,av,Tf,dt,dtpoints):
    nsteps=int(Tf/dt)
    npoints=Tf//dtpoints+1
    skipsteps=int(dtpoints/dt)
    y=y0.copy()
    fprev=f2(y0,t,K,rho,aw,av)
    y+=dt*fprev
    ylist=np.zeros((len(y0),npoints))
    ylist[:,0],ylist[:,1]=y0,y
    for i in range(1,nsteps):
        fnew=f2(y,t,K,rho,aw,av)
        y+=dt*(3/2*fnew-fprev/2)
        fprev=fnew
        if (i+1)%skipsteps==0:
            print((i+1)//skipsteps)
            ylist[:,(i+1)//skipsteps]=y
    return ylist

def F(x,sigma0):
    return erf(x/(np.sqrt(2)*sigma0))

def p(x,sigma0):
    return 2*np.exp(-x**2/(2*sigma0**2))/(np.sqrt(2*np.pi)*sigma0)

def pk(x,k,K,sigma0):
    return (K-k+1)*binom(K,k-1)*(1-F(x,sigma0))**(k-1)*F(x,sigma0)**(K-k)*p(x,sigma0)

def vaverage(K,sigma0):
    v=np.zeros(K)
    for k in range(1,K+1):
        v[k-1],err=quad(lambda x:x*pk(x,k,K,sigma0),0,np.inf)
    return v

K=10
#rho=0.1
aw=0.5
av=0.5
sigma0=1
Tf=10000
nsteps=Tf+1
T=np.linspace(0,Tf,nsteps)
#t=rho
#t1,t2=rho,rho
r0=np.zeros(K)
#qr0=rho*sigma0*np.eye(K)
#qi0=(1-rho)*sigma0*np.eye(K)
v0=vaverage(K,sigma0)
#y0=to_y(qr0,qi0,r0,v0,K)

rholist=[0.001,0.01,0.1,0.2,0.5,0.9,1]
errs,qrsols,qisols,rsols,vsols=[],[],[],[],[]
for rho in rholist:
    print(rho)
    t1,t2=rho,rho
    qr0=rho*sigma0*np.eye(K)
    qi0=(1-rho)*sigma0*np.eye(K)
    y0=to_y(qr0,qi0,r0,v0,K)
    ysol1=solve_ivp(f,(0,Tf),y0,t_eval=T,args=(t1,K,rho,aw,av))
    qrsol1,qisol1,rsol1,vsol1=np.zeros((K,K,nsteps)),np.zeros((K,K,nsteps)),np.zeros((K,nsteps)),np.zeros((K,nsteps))
    for i in range(nsteps):
        qrsol1[:,:,i],qisol1[:,:,i],rsol1[:,i],vsol1[:,i]=from_y(ysol1.y[:,i],K)
    qsol1=qrsol1+qisol1
    err1=error(qsol1,rsol1,vsol1,t1,K)
    #y0_2=to_y(qr0,qi0,r0,v0,K)
    y0twv=to_y(qrsol1[:,:,-1],qisol1[:,:,-1],r0,vsol1[:,-1],K)
    y0tw=to_y(qrsol1[:,:,-1],qisol1[:,:,-1],r0,v0,K)
    # ysol2=Euler(y0_2,t2,K,rho,aw,av,Tf,dt,dtpoints)
    # qrsol2,qisol2,rsol2,vsol2=np.zeros((K,K,nsteps)),np.zeros((K,K,nsteps)),np.zeros((K,nsteps)),np.zeros((K,nsteps))
    # for i in range(nsteps):
    #     qrsol2[:,:,i],qisol2[:,:,i],rsol2[:,i],vsol2[:,i]=from_y(ysol2[:,i],K)
    # qsol2=qrsol2+qisol2
    # err2=error(qsol2,rsol2,vsol2,t2,K)
    ysoltwv=solve_ivp(f,(0,Tf),y0twv,t_eval=T,args=(t2,K,rho,aw,av))
    qrsoltwv,qisoltwv,rsoltwv,vsoltwv=np.zeros((K,K,nsteps)),np.zeros((K,K,nsteps)),np.zeros((K,nsteps)),np.zeros((K,nsteps))
    for i in range(nsteps):
        qrsoltwv[:,:,i],qisoltwv[:,:,i],rsoltwv[:,i],vsoltwv[:,i]=from_y(ysoltwv.y[:,i],K)
    qsoltwv=qrsoltwv+qisoltwv
    errtwv=error(qsoltwv,rsoltwv,vsoltwv,t2,K)
    ysoltw=solve_ivp(f,(0,Tf),y0tw,t_eval=T,args=(t2,K,rho,aw,av))
    qrsoltw,qisoltw,rsoltw,vsoltw=np.zeros((K,K,nsteps)),np.zeros((K,K,nsteps)),np.zeros((K,nsteps)),np.zeros((K,nsteps))
    for i in range(nsteps):
        qrsoltw[:,:,i],qisoltw[:,:,i],rsoltw[:,i],vsoltw[:,i]=from_y(ysoltw.y[:,i],K)
    qsoltw=qrsoltw+qisoltw
    errtw=error(qsoltw,rsoltw,vsoltw,t2,K)
    errs.append([err1,errtwv,errtw])
    qrsols.append([qrsol1,qrsoltwv,qrsoltw])
    qisols.append([qisol1,qisoltwv,qisoltw])
    rsols.append([rsol1,rsoltwv,rsoltw])
    vsols.append([vsol1,vsoltwv,vsoltw])
errs=np.array(errs)
qrsols=np.array(qrsols)
qisols=np.array(qisols)
rsols=np.array(rsols)
vsols=np.array(vsols)

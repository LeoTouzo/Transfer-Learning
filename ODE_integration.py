# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 10:43:39 2021

@author: leoto
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.special import erf

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

def f(y,time,t,K,rho,aw,av):
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

K=5
rho=0.1
aw=0.5
av=0.5
sigma0=1
Tf=500
nsteps=101
T=np.linspace(0,Tf,nsteps)
# t=rho
# r0=np.zeros(K)
# qr0=rho*sigma0*np.ones((K,K))
# qi0=(1-rho)*sigma0*np.ones((K,K))
# v0=np.random.normal(0,sigma0,K)
y0=to_y(qr0,qi0,r0,v0,K)

ysol=odeint(f,y0,T,(t,K,rho,aw,av))
qrsol,qisol,rsol,vsol=np.zeros((K,K,nsteps)),np.zeros((K,K,nsteps)),np.zeros((K,nsteps)),np.zeros((K,nsteps))
for i in range(nsteps):
    qrsol[:,:,i],qisol[:,:,i],rsol[:,i],vsol[:,i]=from_y(ysol[i],K)
qsol=qrsol+qisol
err=error(qsol,rsol,vsol,t,K)
plt.plot(T,err,linestyle='dashed',label='predicted error')
plt.xlabel('time')
plt.ylabel('error')
plt.yscale('log')
plt.legend()
plt.figure()
plt.plot(T,qrsol[0,0],label='qr_1,1')
plt.plot(T,qisol[0,0],label='qi_1,1')
plt.legend()
plt.xlabel('time')
plt.figure()
plt.plot(T,rsol[0])
plt.xlabel('time')
plt.ylabel('r1')
plt.figure()
plt.plot(T,vsol[0])
plt.xlabel('time')
plt.ylabel('v1')

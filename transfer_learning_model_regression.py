# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 16:41:01 2021

@author: leoto
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

#error function activation
# def g(x):
#     return erf(x/np.sqrt(2))

# def dg(x):
#     return np.sqrt(2/np.pi)*np.exp(-x*x/2)

#ReLU activation
def g(x):
    return (x+np.abs(x))/2

def dg(x):
    return x>0

N=100
alpha1,alpha2=300,300
M1,M2=alpha1*N,alpha2*N
rho=0.1
K=10
sigma0=0.002
eta=0.2
# sigmaw=0.0002
# sigmav=0.0002
# eta_w=0.2
# eta_v=0.2
errorstep=10

#b=np.random.binomial(n=1,p=rho,size=N)
b=np.zeros(N)
for i in range(10):
    b[i]=1
J1=np.random.normal(0,1,N)*b
J2=np.random.normal(0,1,N)*b
W0=np.random.normal(0,sigma0,(K,N))
v0=np.random.normal(0,sigma0,K)

def eval_teacher(J,x):
    return np.dot(J,x)/np.sqrt(N)

def eval_student(W,v,x):
    s=0
    for k in range(len(v)):
        s+=v[k]*g(np.dot(W[k],x)/np.sqrt(N))
    return s

testsize=1000
xtest=np.random.normal(0,1,(N,testsize))
y1test=eval_teacher(J1,xtest)
y2test=eval_teacher(J2,xtest)

def test_error(xtest,ytest,W,v):
    return np.sum((ytest-eval_student(W,v,xtest))**2)/len(ytest)

def train_student(J,W0,v0,eta_w,eta_v,M,xtest,ytest,errorstep):
    N=len(J)
    W,v=W0.copy(),v0.copy()
    testerror=[test_error(xtest,ytest,W,v)]
    Wlist,vlist=[W.copy()],[v.copy()]
    for t in range(M):
        x=np.random.normal(0,1,N)
        lam=W.dot(x)/np.sqrt(N)
        delta=eval_student(W,v,x)-eval_teacher(J,x)
        W-=eta_w/np.sqrt(N)*delta*np.dot(np.array([v*dg(lam)]).T,np.array([x]))
        v-=eta_v/np.sqrt(N)*delta*g(lam)
        if t%errorstep==0:
            testerror.append(test_error(xtest,ytest,W,v))
            Wlist.append(W.copy())
            vlist.append(v.copy())
    return W,v,testerror,np.array(Wlist),np.array(vlist) #returns the final weights and the evolution of error and weights during training

def train_student_soft_transfer(J,W0,v0,W_old,v_old,beta_w,beta_v,eta_w,eta_v,M,xtest,ytest,errorstep):
    N=len(J)
    W,v=W0.copy(),v0.copy()
    testerror=[test_error(xtest,ytest,W,v)]
    for t in range(M):
        x=np.random.normal(0,1,N)
        lam=W.dot(x)/np.sqrt(N)
        delta=eval_student(W,v,x)-eval_teacher(J,x)
        W-=eta_w/np.sqrt(N)*delta*np.dot(np.array([v*dg(lam)]).T,np.array([x]))+beta_w*(W-W_old)
        v-=eta_v/np.sqrt(N)*delta*g(lam)+beta_v*(v-v_old)
        if t%errorstep==0:
            testerror.append(test_error(xtest,ytest,W,v))
    return W,v,testerror

step=np.linspace(0,M1,M1//errorstep+1)
W1,v1,testerror1,W1list,v1list=train_student(J1,W0,v0,eta,eta,M1,xtest,y1test,errorstep)
plt.plot(step,testerror1)
plt.xlabel('training step')
plt.ylabel('test error')
plt.yscale('log')
plt.figure()

#To optimize initialisation and learning rate
# step=np.linspace(0,M1,M1//errorstep+1)
# sigma0=[0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01]
# eta=[0.05,0.1,0.2,0.5]
# W1s=[]
# err=[]
# for et in eta:
#     print(et)
#     for s0 in sigma0:
#         print(s0)
#         W0=np.random.normal(0,s0,(K,N))
#         v0=np.random.normal(0,s0,K)
#         W1,v1,testerror1,W1list,v1list=train_student(J1,W0,v0,et,et,M1,xtest,y1test,errorstep)
#         plt.plot(step,testerror1,label='sigma0='+str(s0)+', eta='+str(et))
#         err.append(testerror1)
#         W1s.append(W1)
#     plt.xlabel('training step')
#     plt.ylabel('test error')
#     plt.yscale('log')
#     plt.legend()
#     plt.figure()

#Hard transfer
# W2_0,v2_0,testerror2_0,W2list_0,v2list_0=train_student(J2,W0,v0,eta,eta,M2,xtest,y2test,errorstep)
# W2_ht,v2_ht,testerror2_ht,W2list_ht,v2list_ht=train_student(J2,W1,v1,eta,eta,M2,xtest,y2test,errorstep)
# W2_htw,v2_htw,testerror2_htw,W2list_htw,v2list_htw=train_student(J2,W1,v0,eta,eta,M2,xtest,y2test,errorstep)
# W2_htv,v2_htv,testerror2_htv,W2list_htv,v2list_htv=train_student(J2,W0,v1,eta,eta,M2,xtest,y2test,errorstep)
# W02=np.random.normal(0,np.std(W1),(K,N))
# v02=np.random.normal(0,np.std(v1),K)
# W2_scaleW1,v2_scaleW1,testerror2_scaleW1,W2list_scaleW1,v2list_scaleW1=train_student(J2,W02,v0,eta,eta,M2,xtest,y2test,errorstep)
# W2_scalev1,v2_scalev1,testerror2_scalev1,W2list_scalev1,v2list_scalev1=train_student(J2,W0,v02,eta,eta,M2,xtest,y2test,errorstep)
# step=np.linspace(0,M2,M2//errorstep+1)
# plt.plot(step,testerror2_0,label='no transfer')
# plt.plot(step,testerror2_ht,label='hard transfer (W and v)')
# plt.plot(step,testerror2_htw,label='hard transfer (W)')
# plt.plot(step,testerror2_htv,label='hard transfer (v)')
# plt.plot(step,testerror2_scaleW1,label='transfer of std(W)')
# plt.plot(step,testerror2_scalev1,label='transfer of std(v)')
# plt.xlabel('training step')
# plt.ylabel('test error')
# plt.yscale('log')
# plt.legend()

#To plot the weights
# index=np.linspace(1,N,N)
# base_vec=np.eye(N)
# for k in range(K):
#     plt.scatter(index,np.dot(base_vec,W1[k]))

# plt.figure()
# plt.plot(step,W2list_0[:,0,maxJ[-1]],label='no transfer')
# plt.plot(step,W2list_ht[:,0,maxJ[-1]],label='hard transfer (W and v)')
# plt.plot(step,W2list_htw[:,0,maxJ[-1]],label='hard transfer (W)')
# plt.plot(step,W2list_htv[:,0,maxJ[-1]],label='hard transfer (v)')
# plt.xlabel('training step')
# plt.ylabel('W')
# plt.legend()

# plt.figure()
# plt.plot(step,W2list_0[:,0,maxJ[0]],label='no transfer')
# plt.plot(step,W2list_ht[:,0,maxJ[0]],label='hard transfer (W and v)')
# plt.plot(step,W2list_htw[:,0,maxJ[0]],label='hard transfer (W)')
# plt.plot(step,W2list_htv[:,0,maxJ[0]],label='hard transfer (v)')
# plt.xlabel('training step')
# plt.ylabel('W')
# plt.legend()

# plt.figure()
# plt.plot(step,v2list_0,label='no transfer')
# plt.plot(step,v2list_ht,label='hard transfer (W and v)')
# plt.plot(step,v2list_htw,label='hard transfer (W)')
# plt.plot(step,v2list_htv,label='hard transfer (v)')
# plt.xlabel('training step')
# plt.ylabel('v')
# plt.legend()

# plt.figure()

#Soft transfer
# W2_0,v2_0,testerror2_0,W2list,v2list=train_student(J2,W0,v0,eta,eta,M2,xtest,y2test,errorstep)
# W2_st,v2_st,testerror2_st=train_student_soft_transfer(J2,W0,v0,W1,v1,1e-4,0,eta,eta,M2,xtest,y2test,errorstep)
# W2_st,v2_st,testerror2_st2=train_student_soft_transfer(J2,W0,v0,W1,v1,0,1e-4,eta,eta,M2,xtest,y2test,errorstep)
# W2_st,v2_st,testerror2_st3=train_student_soft_transfer(J2,W0,v0,W1,v1,1e-4,1e-4,eta,eta,M2,xtest,y2test,errorstep)
# step=np.linspace(0,M2,M2//errorstep+1)
# plt.plot(step,testerror2_0,label='no transfer')
# plt.plot(step,testerror2_st,label='beta_w=1e-4')
# plt.plot(step,testerror2_st2,label='beta_v=1e-4')
# plt.plot(step,testerror2_st3,label='beta_w=beta_v=1e-4')
# plt.xlabel('training step')
# plt.ylabel('test error')
# plt.yscale('log')
# plt.legend()

#To average over several initialisations
# n_avg=10
# step=np.linspace(0,M2,M2//errorstep+1)
# err2_0,err2_ht,err2_htw,err2_htv=0,0,0,0
# best=np.zeros(4)
# for i in range(n_avg):
#     print(i)
#     # b=np.random.binomial(n=1,p=rho,size=N)
#     # J1=np.random.normal(0,1,N)*b
#     # J2=np.random.normal(0,1,N)*b
#     testsize=1000
#     xtest=np.random.normal(0,1,(N,testsize))
#     y1test=eval_teacher(J1,xtest)
#     y2test=eval_teacher(J2,xtest)
#     W0=np.random.normal(0,sigma0,(K,N))
#     v0=np.random.normal(0,sigma0,K)
#     W1,v1,testerror1,W1list,v1list=train_student(J1,W0,v0,eta,eta,M1,xtest,y1test,errorstep)
#     W2_0,v2_0,testerror2_0,W2list_0,v2list_0=train_student(J2,W0,v0,eta,eta,M2,xtest,y2test,errorstep)
#     W2_ht,v2_ht,testerror2_ht,W2list_ht,v2list_ht=train_student(J2,W1,v1,eta,eta,M2,xtest,y2test,errorstep)
#     W2_htw,v2_htw,testerror2_htw,W2list_htw,v2list_htw=train_student(J2,W1,v0,eta,eta,M2,xtest,y2test,errorstep)
#     W2_htv,v2_htv,testerror2_htv,W2list_htv,v2list_htv=train_student(J2,W0,v1,eta,eta,M2,xtest,y2test,errorstep)
#     err2_0+=np.array(testerror2_0)/n_avg
#     err2_ht+=np.array(testerror2_ht)/n_avg
#     err2_htw+=np.array(testerror2_htw)/n_avg
#     err2_htv+=np.array(testerror2_htv)/n_avg
#     finalerror=[testerror2_0[-1],testerror2_ht[-1],testerror2_htw[-1],testerror2_htv[-1]]
#     best[np.argmin(finalerror)]+=1
#     plt.plot(step,testerror2_0,label='no transfer')
#     plt.plot(step,testerror2_ht,label='hard transfer (W and v)')
#     plt.plot(step,testerror2_htw,label='hard transfer (W)')
#     plt.plot(step,testerror2_htv,label='hard transfer (v)')
#     plt.xlabel('training step')
#     plt.ylabel('test error')
#     plt.yscale('log')
#     plt.legend()
#     plt.figure()
# # sigma0=np.std(W1)
# # W0=np.random.normal(0,sigma0,(K,N))
# # v0=np.random.normal(0,sigma0,K)
# # W2_02,v2_02,testerror2_02,W2list_02,v2list_02=train_student(J2,W0,v0,eta_w,eta_v,M2,xtest,y2test,errorstep)
# plt.plot(step,err2_0,label='no transfer')
# plt.plot(step,err2_ht,label='hard transfer (W and v)')
# plt.plot(step,err2_htw,label='hard transfer (W)')
# plt.plot(step,err2_htv,label='hard transfer (v)')
# # plt.plot(step,testerror2_02,label='transfer of std(W)')
# plt.xlabel('training step')
# plt.ylabel('test error')
# plt.yscale('log')
# plt.legend()

# for i in range(7):
#     for j in range(4):
#         plt.plot(step,err[i,j],label='eta_v='+str(paramv[i,j])+', eta_w='+str(paramw[i,j]))
#     plt.xlabel('training step')
#     plt.ylabel('test error')
#     plt.yscale('log')
#     plt.legend()
#     plt.figure()

# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 09:33:33 2021

@author: leoto
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

#error function activation
def g(x):
    return erf(x/np.sqrt(2))

def dg(x):
    return np.sqrt(2/np.pi)*np.exp(-x*x/2)

#ReLU activation
# def g(x):
#     return (x+np.abs(x))/2

# def dg(x):
#     return x>0
    
def logloss(x):
    return np.log(1+np.exp(-x))

def dlogloss(x):
    return -1/(1+np.exp(x))

N=100
alpha1,alpha2=10000,500
M1,M2=alpha1*N,alpha2*N
rho=0.1
K=5
sigma0=0.01
eta=0.01
sigmaw=sigma0
sigmav=sigma0
eta_w=eta
eta_v=eta
errorstep=10

b=np.random.binomial(n=1,p=rho,size=N)
J1=np.random.normal(0,1,N)*b
J2=np.random.normal(0,1,N)*b
W0=np.random.normal(0,sigmaw,(K,N))
v0=np.random.normal(0,sigmav,K)

def eval_teacher(J,x):
    return np.sign(np.dot(J,x))

def eval_student(W,v,x):
    s=0
    for k in range(len(v)):
        s+=v[k]*g(np.dot(W[k],x)/np.sqrt(N))
    return s

testsize=1000
xtest=np.random.normal(0,1,(N,testsize))
y1test=eval_teacher(J1,xtest)
y2test=eval_teacher(J2,xtest)

def test_loss(xtest,ytest,W,v):
    return np.sum((ytest-eval_student(W,v,xtest))**2)/len(ytest)

def test_error(xtest,ytest,W,v):
    return -np.sum(ytest*np.sign(eval_student(W,v,xtest)))/2/len(ytest)+1/2

def train_student(J,W0,v0,eta_w,eta_v,M,xtest,ytest,errorstep):
    N=len(J)
    W,v=W0.copy(),v0.copy()
    testloss=[test_loss(xtest,ytest,W,v)]
    #testerror=[test_error(xtest,ytest,W,v)]
    Wlist,vlist=[W.copy()],[v.copy()]
    for t in range(M):
        x=np.random.normal(0,1,N)
        lam=W.dot(x)/np.sqrt(N)
        delta=eval_student(W,v,x)-eval_teacher(J,x)
        W-=eta_w/np.sqrt(N)*delta*np.dot(np.array([v*dg(lam)]).T,np.array([x]))
        v-=eta_v/np.sqrt(N)*delta*g(lam)
        if t%errorstep==0:
            testloss.append(test_loss(xtest,ytest,W,v))
            #testerror.append(test_error(xtest,ytest,W,v))
            Wlist.append(W.copy())
            vlist.append(v.copy())
    return W,v,testloss,np.array(Wlist),np.array(vlist) #,testerror

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

def train_student_logloss(J,W0,v0,eta_w,eta_v,M,xtest,ytest,errorstep):
    N=len(J)
    W,v=W0.copy(),v0.copy()
    testloss=[test_loss(xtest,ytest,W,v)]
    #testerror=[test_error(xtest,ytest,W,v)]
    Wlist,vlist=[W.copy()],[v.copy()]
    for t in range(M):
        x=np.random.normal(0,1,N)
        y=eval_teacher(J,x)
        lam=W.dot(x)/np.sqrt(N)
        delta=y*dlogloss(eval_student(W,v,x)*y)
        W-=eta_w/np.sqrt(N)*delta*np.dot(np.array([v*dg(lam)]).T,np.array([x]))
        v-=eta_v/np.sqrt(N)*delta*g(lam)
        if t%errorstep==0:
            testloss.append(test_loss(xtest,ytest,W,v))
            #testerror.append(test_error(xtest,ytest,W,v))
            Wlist.append(W.copy())
            vlist.append(v.copy())
    return W,v,testloss,np.array(Wlist),np.array(vlist) #,testerror
    
W1,v1,testerror1,W1list,v1list=train_student_logloss(J1,W0,v0,eta_w,eta_v,M1,xtest,y1test,errorstep)
step=np.linspace(0,M1,M1//errorstep+1)
plt.plot(step,testerror1)
plt.xlabel('training step')
plt.ylabel('test error')
plt.yscale('log')
#plt.figure()

#To optimize initialisation and learning rate
# step=np.linspace(0,M1,M1//errorstep+1)
# sigma0=[0.01,0.02,0.05,0.1,0.2,0.5]
# eta=[0.001,0.002,0.005,0.01,0.02,0.05]
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
# W2_0,v2_0,testerror2_0,W2list_0,v2list_0=train_student(J2,W0,v0,eta_w,eta_v,M2,xtest,y2test,errorstep)
# W2_ht,v2_ht,testerror2_ht,W2list_ht,v2list_ht=train_student(J2,W1,v1,eta_w,eta_v,M2,xtest,y2test,errorstep)
# W2_htw,v2_htw,testerror2_htw,W2list_htw,v2list_htw=train_student(J2,W1,v0,eta_w,eta_v,M2,xtest,y2test,errorstep)
# W2_htv,v2_htv,testerror2_htv,W2list_htv,v2list_htv=train_student(J2,W0,v1,eta_w,eta_v,M2,xtest,y2test,errorstep)
# step=np.linspace(0,M2,M2//errorstep+1)
# plt.plot(step,testerror2_0,label='no transfer')
# plt.plot(step,testerror2_ht,label='hard transfer (W and v)')
# plt.plot(step,testerror2_htw,label='hard transfer (W)')
# plt.plot(step,testerror2_htv,label='hard transfer (v)')
# plt.xlabel('training step')
# plt.ylabel('test error')
# plt.yscale('log')
# plt.legend()
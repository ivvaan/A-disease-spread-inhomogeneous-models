import numpy as np
from scipy.integrate import odeint

modeling_period=9.0
t_step=0.003
m_t=np.arange(0.0, modeling_period+t_step, t_step,dtype=np.float64)

def HCS_SIR(r0_0,I0,τ,vs0,vb0):
    '''
    Heterogeneous in Catching and Spreading
    SIR model
    '''
    def M(mu,v):return np.exp(mu+0.5*v)
    def a_s(mus,C,vs):
        Cpvs=C+vs
        return M(mus,vs)*(1.0+Cpvs*(1.0+0.5*Cpvs))
    M0=M(0.0,vs0+vb0)
    a_s0=r0_0*a_s(0.0,0.0,vs0)/(r0_0-1.0+1.0/τ)
    α=r0_0/M0/a_s0
    
    def ODE(y,t):
        S,I,A,mus,mub,C,vs,vb=y
        Cpvs=C+vs
        Cpvb=C+vb 
        minus_αMA=-α*M(mus+mub,Cpvs+Cpvb)*A

        dS = minus_αMA*S
        dI = -dS - I
        dA= -a_s(mus,C,vs)*dS-A/τ
        dmus = minus_αMA*Cpvs
        dmub = minus_αMA*Cpvb
        dC = dmus*Cpvb
        dvs = dmus*Cpvs
        dvb = dmub*Cpvb
        return [dS, dI,dA, dmus,dmub,dC,dvs,dvb]
    
    A0=a_s0*I0
    y0=[1.0-I0,I0,A0,0.0,0.0,0.0,vs0,vb0]
    
    sol = odeint(ODE,y0, m_t).T
    S,I,A,mus,mub,C,vs,vb=sol
    #dS=-ODE(sol,m_t)[0]
    return S,I

def HCS_SIR2(r0_0,I0,τ,vs0,vb0):
    '''
    Heterogeneous in Catching and Spreading
    SIR model splitted Infectious
    '''
    def M(mu,v):return np.exp(mu+0.5*v)
    def a_s(mus,C,vs):
        Cpvs=C+vs
        return M(mus,vs)*(1.0+Cpvs*(1.0+0.5*Cpvs))
    M0=M(0.0,vs0+vb0)
    r0=(r0_0-1.0+1.0/τ)
    a_s0=a_s(0.0,0.0,vs0)
    α=r0/M0/a_s0
    γ_a=1/τ
    γ_p=1/(1-τ)
    def ODE(y,t):
        S,I_a,I_p,A,mus,mub,C,vs,vb=y
        Cpvs=C+vs
        Cpvb=C+vb 
        minus_αMA=-α*M(mus+mub,Cpvs+Cpvb)*A

        dS = minus_αMA*S
        dI_a = -dS - γ_a*I_a
        dI_p = γ_a*I_a -  γ_p*I_p
        dA= -a_s(mus,C,vs)*dS-γ_a*A
        dmus = minus_αMA*Cpvs
        dmub = minus_αMA*Cpvb
        dC = dmus*Cpvb
        dvs = dmus*Cpvs
        dvb = dmub*Cpvb
        return [dS, dI_a,dI_p,dA, dmus,dmub,dC,dvs,dvb]
    
    r_ratio=r0_0/r0
    A0=a_s0*I0*r_ratio
    y0=[1.0-I0,I0*r_ratio,I0*(1.0-r_ratio),A0,0.0,0.0,0.0,vs0,vb0]
    
    sol = odeint(ODE,y0, m_t).T
    S,I_a,I_p,A,mus,mub,C,vs,vb=sol
    #dS=-ODE(sol,m_t)[0]
    return S,I_a+I_p
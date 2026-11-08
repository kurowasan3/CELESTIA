import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
# symbolic
k,n1,n2,V0,T,R,t,x=sp.symbols('k n1 n2 V0 T R t x',positive=True)
Phi_T=V0*T/n1
I1T=Phi_T/(k*n1)
tau=k*n2**2/R
I2p=Phi_T/(k*n2)              # flux continuity, I1=0
print('I1(T)=',sp.simplify(I1T),' I2(T+)=',sp.simplify(I2p), ' ratio',sp.simplify(I2p/I1T))
V1=-(n1/n2)*R*I2p            # V1=n1 dPhi/dt, n2 dPhi/dt=-R I2
print('switch voltage', sp.simplify(V0-V1), ' -> V0(1+T/tau)?', sp.simplify(V0-V1-V0*(1+T/tau)))
heat=sp.integrate(R*(I2p*sp.exp(-t/tau))**2,(t,0,sp.oo))
work=sp.integrate(V0*V0*t/(k*n1**2),(t,0,T))
print('heat',sp.simplify(heat),' work',sp.simplify(work), sp.simplify(heat-work))
# periodic steady state
dP=V0*T/n1; Ps=sp.symbols('Ps')
Ps_sol=sp.solve(sp.Eq(Ps,(Ps+dP)*sp.exp(-x)),Ps)[0]
E_cyc=sp.simplify(((Ps_sol+dP)**2-Ps_sol**2)/(2*k))
target=V0**2*T**2/(2*k*n1**2)*sp.cosh(x/2)/sp.sinh(x/2)
print('E per cycle check',sp.simplify((E_cyc-target).rewrite(sp.exp)))
# work by source per cycle directly
Wsrc=sp.integrate(V0*(Ps_sol+V0*t/n1)/(k*n1),(t,0,T))
print('source work check',sp.simplify((Wsrc-target).rewrite(sp.exp)))
# numeric simulation with slightly leaky coupling (independent of flux argument)
kk,N1,N2,v0,TT,RR=1e-6,100,50,10.0,1e-3,5.0
L1=kk*N1**2; L2=kk*N2**2*(1+1e-4); M=kk*N1*N2
Ai=np.linalg.inv(np.array([[L1,M],[M,L2]]))
def rhs(tt,y,on):
    i1,i2=y
    Rs=0.0 if on else 1e4
    Rd=0.0 if i2>0 else 1e4
    return Ai@np.array([(v0 if on else 0)-Rs*i1,-(RR+Rd)*i2])
tau_n=kk*N2**2/RR
s=solve_ivp(rhs,[0,TT],[0,0],args=(True,),method='Radau',rtol=1e-9,atol=1e-12)
print('num I1(T)',s.y[0,-1],' exact',v0*TT/(kk*N1**2),' I2',s.y[1,-1])
s2=solve_ivp(rhs,[TT,TT+2e-5],s.y[:,-1],args=(False,),method='Radau',rtol=1e-9,atol=1e-12,max_step=1e-7,dense_output=True)
print('num I2(T+2e-5)',s2.y[1,-1],' exact',(N1/N2)*v0*TT/(kk*N1**2)*np.exp(-2e-5/tau_n))
# switch voltage = V0 - V1, V1 = L1 dI1/dt + M dI2/dt ; evaluate after spike (t=T+5us)
yy=s2.sol(TT+5e-6); d=rhs(TT+5e-6,yy,False); V1=L1*d[0]+M*d[1]
print('switch V num',v0-V1,' exact (at T+5us)',v0+RR*TT*v0/(kk*N2**2)*np.exp(-5e-6/tau_n))
Toff=0.6e-3; per=TT+Toff; y=np.array([0.,0.]); E=[]
for c in range(30):
    a=solve_ivp(rhs,[0,TT],y,args=(True,),method='Radau',rtol=1e-9,atol=1e-12)
    b=solve_ivp(rhs,[0,Toff],a.y[:,-1],args=(False,),method='Radau',rtol=1e-9,atol=1e-12,max_step=1e-6,dense_output=True)
    tt=np.linspace(0,Toff,60001); i2=np.clip(b.sol(tt)[1],0,None)
    E.append(np.trapezoid(RR*i2**2,tt)); y=b.y[:,-1]
xx=Toff/tau_n
print('Pavg num',E[-1]/per,' formula',v0**2*TT**2/(2*kk*N1**2)/np.tanh(xx/2)/per)

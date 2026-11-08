# O20243-1 neutron moderation: elastic scattering off nucleus mass A m, isotropic in CM
import sympy as sp, numpy as np
A,v,Phi=sp.symbols('A v Phi',positive=True)
VG=v/(A+1); u=A*v/(A+1)           # CM velocity, neutron speed in CM
vx=VG+u*sp.cos(Phi); vy=u*sp.sin(Phi)
ratio=sp.simplify((vx**2+vy**2)/v**2); print('E\'/E =',sp.factor(ratio))
print('check form:',sp.simplify(ratio-(A**2+2*A*sp.cos(Phi)+1)/(A+1)**2))
# independent check: lab-frame momentum+energy conservation for head-on (Phi=pi)
v1,V2=sp.symbols('v1 V2')
s=sp.solve([sp.Eq(v,v1+A*V2),sp.Eq(v**2,v1**2+A*V2**2)],[v1,V2])
print('head-on solutions',s, ' -> alpha=',sp.simplify(ratio.subs(Phi,sp.pi)))
al=((A-1)/(A+1))**2
x=sp.symbols('x',positive=True)
mean=sp.integrate(x,(x,al,1))/(1-al); print('mean ratio',sp.simplify(mean))
xi=-sp.integrate(sp.log(x),(x,al,1))/(1-al)
xi_ans=1+al*sp.log(al)/(1-al)
print('xi check (A=12):',float(xi.subs(A,12)),float(xi_ans.subs(A,12)))
print('xi A->1 limit:',sp.limit(xi_ans,A,1))
# Monte Carlo
rng=np.random.default_rng(1)
for Av in [1,12]:
    c=rng.uniform(-1,1,2_000_000); rr=(Av**2+2*Av*c+1)/(Av+1)**2
    print('A',Av,'MC mean',rr.mean(),'MC <ln>',np.log(rr).mean(),' alpha',((Av-1)/(Av+1))**2)
L=np.log(2.0e6/0.025); print('ln ratio',L,'N_H',L/1,'N_C',L/float(xi_ans.subs(A,12)),'N_H naive',L/np.log(2))
# with given table values ln2=0.693, ln10=2.303, ln11=2.398, ln13=2.565
lnal=2*(2.398-2.565); alv=(11/13)**2; print('alpha',alv,'xi from table',1+alv*lnal/(1-alv))
print('ln(8e7) from table',3*0.693+7*2.303)
# MC simulation of slowing down histories (A=1): median collisions to go 2MeV->0.025eV
for Av in [1,12]:
    N=20000; E=np.full(N,2.0e6); cnt=np.zeros(N,int); alive=np.ones(N,bool)
    while alive.any():
        c=rng.uniform(-1,1,alive.sum()); E[alive]*= (Av**2+2*Av*c+1)/(Av+1)**2; cnt[alive]+=1
        alive=E>0.025
    print('A',Av,'mean collisions',cnt.mean(),'median',np.median(cnt))

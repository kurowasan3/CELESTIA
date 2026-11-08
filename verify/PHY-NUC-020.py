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
# with given table values ln2=0.693, ln10=2.303, ln11=2.398, ln13=2.565
lnal=2*(2.398-2.565); alv=(11/13)**2; print('alpha',alv,'xi from table',1+alv*lnal/(1-alv))

# ---- new (5): N'/N, Jensen inequality, heavy-nucleus limit ----
c_=(1+al)/2
diff_=sp.simplify(xi_ans+sp.log(c_))   # xi - log(2/(1+alpha)) > 0 ?
for Av in [1.0001,1.5,2,4,12,50,200,2000]:
    val=float(diff_.subs(A,Av)); print(' A=%g  xi - log(2/(1+alpha)) = %.3e   N\'/N = %.5f'%(Av,val,float((xi_ans/(-sp.log(c_))).subs(A,Av))))
dl=sp.symbols('delta',positive=True)
xi_d=1+(1-dl)*sp.log(1-dl)/dl
print('xi series in delta:',sp.series(xi_d,dl,0,3))
print('log(2/(1+alpha)) series:',sp.series(-sp.log(1-dl/2),dl,0,3))
print('N\'/N series:',sp.series(xi_d/(-sp.log(1-dl/2)),dl,0,2))
eA=sp.symbols('e_A',positive=True)  # 1/A
print('xi in 1/A:',sp.series(xi_ans.subs(A,1/eA),eA,0,2))
print('hydrogen N\'/N =',1/np.log(2))
# Monte Carlo: product of many ratios -> typical log decrement per collision = xi, not -log(mean)
for Av in [1,12]:
    Nc=400; trials=4000
    c=rng.uniform(-1,1,(trials,Nc)); rr=(Av**2+2*Av*c+1)/(Av+1)**2
    print(' A',Av,' median of -log(prod)/N =',np.median(-np.log(rr).sum(1)/Nc),' xi=',float(xi_ans.subs(A,Av)) if Av>1 else 1.0,' -log(mean)=',-np.log(rr.mean()))

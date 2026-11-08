# M54-1: straight-line process, unknown Cv inferred from heat-switch point
import sympy as sp, numpy as np
P0,V0,c,V,nR=sp.symbols('P0 V0 c V nR',positive=True)  # c = Cv/R
P=P0*(5-V/V0)
T=P*V/nR
Vt=sp.solve(sp.diff(T,V),V)[0]; print('Tmax at',Vt,' Tmax',sp.simplify(T.subs(V,Vt)))
dQdV=c*sp.diff(P*V,V)+P   # dQ = c d(PV) + P dV
Vs=sp.solve(dQdV,V)[0]; print('switch V',sp.simplify(Vs))
cs=sp.solve(sp.Eq(Vs,sp.Rational(35,12)*V0),c); print('c from Vs=35/12 V0:',cs)
cc=cs[0]
print('Vs - Vt =',sp.simplify(Vs-Vt))
Qf=lambda a,b: sp.integrate(dQdV.subs(c,cc),(V,a,b))
Qabs=sp.simplify(Qf(V0,Vs.subs(c,cc))); Qrel=sp.simplify(Qf(Vs.subs(c,cc),4*V0))
print('Q absorbed',Qabs,float(Qabs/(P0*V0)),' Q released(signed)',Qrel,float(Qrel/(P0*V0)))
# check via first law ΔU + W
dU=cc*(P0*4*V0-4*P0*V0); Wt=sp.integrate(P,(V,V0,4*V0)); print('net Q',sp.simplify(Qabs+Qrel),' dU+W',sp.simplify(dU+Wt))
# heat between Vt and Vs while T falls
Qmid=sp.simplify(Qf(Vt,Vs.subs(c,cc))); print('Q in Tmax->Vs',Qmid, ' dU there',sp.simplify(cc*(P*V).subs(V,Vs.subs(c,cc))-cc*(P*V).subs(V,Vt)))
# numeric check: discrete heat increments for monatomic vs diatomic
for cn in [1.5,2.5]:
    Vg=np.linspace(1,4,300001); Pg=5-Vg; U=cn*Pg*Vg
    dQ=np.diff(U)+0.5*(Pg[1:]+Pg[:-1])*np.diff(Vg)
    idx=np.where(dQ<0)[0][0]; print('c',cn,'switch numeric',Vg[idx],'formula',5*(cn+1)/(2*cn+1))

# M55-1: horizontal cylinder, spring engaged after gap d; heat & kink in T-Q graph
import sympy as sp, numpy as np
P0,S,L,d,k,n,R,T0,x=sp.symbols('P0 S L d k n R T0 x',positive=True)
nR=P0*S*L/T0
T1=sp.simplify(P0*S*(L+d)/nR); print('T1',T1)
P2=P0+k*d/S; V2=S*(L+2*d)
dU=sp.Rational(3,2)*(P2*V2-P0*S*L)
W_int=P0*S*d+sp.integrate(P0+k*(x-d)/S,(x,d,2*d))*S
W_en=P0*S*2*d+k*d**2/2
print('W integral - energy',sp.simplify(W_int-W_en))
Q=sp.expand(dU+W_en); print('Q',Q)
# slopes dT/dQ before/after contact
# after contact: P=P0+k(x-d)/S, V=S(L+x); T=PV/nR; Q(x)=3/2 PV + ∫P dV
Px=P0+k*(x-d)/S; Vx=S*(L+x)
dTdx=sp.diff(Px*Vx,x)/nR; dQdx=sp.Rational(3,2)*sp.diff(Px*Vx,x)+Px*S
after=sp.simplify((dTdx/dQdx).subs(x,d)); before=1/(sp.Rational(5,2)*nR)
ratio=sp.simplify(after/before); print('ratio',ratio)
a=sp.symbols('a',positive=True)
target=5*(P0+a)/(5*P0+3*a)
print('check',sp.simplify(ratio-target.subs(a,k*(L+d)/S)))
print('k->oo',sp.limit(ratio,k,sp.oo))
# numeric: march heat with state switch
P0n,Sn,Ln,dn,kn,T0n=1e5,1e-2,0.20,0.05,2.0e4,300.
nRn=P0n*Sn*Ln/T0n
xs=np.linspace(0,2*dn,200001)
P=np.where(xs<dn,P0n,P0n+kn*(xs-dn)/Sn); Vv=Sn*(Ln+xs)
U=1.5*P*Vv; Wc=np.concatenate([[0],np.cumsum(0.5*(P[1:]+P[:-1])*np.diff(Vv))])
Qn=U-U[0]+Wc
print('Q numeric',Qn[-1],' formula',float(Q.subs({P0:P0n,S:Sn,L:Ln,d:dn,k:kn})))
T=P*Vv/nRn; i=np.searchsorted(xs,dn)
sb=(T[i-1]-T[i-50])/(Qn[i-1]-Qn[i-50]); sa=(T[i+50]-T[i+1])/(Qn[i+50]-Qn[i+1])
print('slope ratio numeric',sa/sb,' formula',float(ratio.subs({P0:P0n,S:Sn,L:Ln,d:dn,k:kn})))

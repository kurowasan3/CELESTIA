# N52-1 Millikan with polarity reversal: radius from difference, charge from sum
import numpy as np, sympy as sp
eta,rho,rho0,g,a,q,Ef,vu,vd=sp.symbols('eta rho rho0 g a q E v_u v_d',real=True)
W=sp.Rational(4,3)*sp.pi*a**3*(rho-rho0)*g
# polarity1: electric force up: q E - W = 6 pi eta a vu (vu up positive)
# polarity2: electric force down: q E + W = 6 pi eta a vd (vd down positive)
sol=sp.solve([sp.Eq(q*Ef-W,6*sp.pi*eta*a*vu),sp.Eq(q*Ef+W,6*sp.pi*eta*a*vd)],[q,rho],dict=True)[0]
print('q=',sp.simplify(sol[q]))
asol=sp.solve(sp.Eq(W,6*sp.pi*eta*a*(vd-vu)/2),a)
print('a=',asol)
# numeric data generation
etan,rhon,rho0n,gn,dn,Vn=1.82e-5,875.,1.20,9.80,6.00e-3,360.
an=9.0e-7; etrue=1.602e-19; En=Vn/dn
Wn=4/3*np.pi*an**3*(rhon-rho0n)*gn; k=6*np.pi*etan*an
data=[]
for n in [4,7,2,5]:
    qn=n*etrue
    vu_=(qn*En-Wn)/k; vd_=(qn*En+Wn)/k
    data.append((float('%.3g'%(vu_*1e5)),float('%.3g'%(vd_*1e5))))
print('data (1e-5 m/s, up-positive vu, down vd):',data)
# solve from rounded data
vus=np.array([d[0] for d in data])*1e-5; vds=np.array([d[1] for d in data])*1e-5
vg=(vds-vus)/2; print('vg each',vg)
vgm=vg.mean()
a_est=np.sqrt(9*etan*vgm/(2*(rhon-rho0n)*gn)); print('a est',a_est)
qs=6*np.pi*etan*a_est*(vus+vds)/(2*En); print('q',qs/1e-19)
diffs=np.diff(np.sort(qs)); print('sorted diffs',diffs/1e-19)
ns=np.round(qs/1.6e-19); print('n',ns)
print('e est',qs.sum()/ns.sum())
# alternative: each pair's own radius
a_i=np.sqrt(9*etan*vg/(2*(rhon-rho0n)*gn)); q_i=6*np.pi*etan*a_i*(vus+vds)/(2*En)
print('per-pair q',q_i/1e-19,' e',q_i.sum()/ns.sum())
# physical check: terminal-velocity relaxation time m/(6 pi eta a)
mn=4/3*np.pi*an**3*rhon; print('tau',mn/(6*np.pi*etan*an))
# solution-path numbers with averaged vg, and (4) balancing voltage for measurement 3
print('vg mean',vgm,' a',a_est,' coef',3*np.pi*etan*a_est/En)
v1,v2=vus[2],vds[2]; print("V' (4) =",Vn*(v2-v1)/(v1+v2))
# independent: true n=2 drop: W/(q/d)
print("V' true", Wn*dn/(2*etrue))

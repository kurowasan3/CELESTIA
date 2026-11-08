# M67-1: string (fundamental) drives closed air column; identify harmonics from tension ratio
import numpy as np
V=340.0; l=0.50; rho=5.0e-4; D=0.630; S1=80.0; S2=222.0
f1=np.sqrt(S1/rho)/(2*l); f2=np.sqrt(S2/rho)/(2*l)
r=f2/f1; print("f1",f1,"f2",f2,"ratio",r)
for m in range(1,6): print(m,(2*m+1)/(2*m-1))
# identify: f1 = (2m-1)V/(4(D+e)), f2=(2m+1)V/(4(D+e)) -> m=2
m=2
Le=(2*m-1)*V/(4*f1); e=Le-D; print("D+e",Le,"e",e)
# alternative independent: use f2 with 5th harmonic
print("D+e from f2 (exact 5/3 ratio)",5*V/(4*f1*5/3), " using measured S2:",5*V/(4*f2))
lam1=V/f1; print("lam1",lam1,"next depth",D+lam1/2)
lam2=4*Le/5; print("lam2",lam2,"V/f2(exact)",V/(f1*5/3))
nodes=[D-k*lam2/2 for k in range(0,3)]; print("nodes",nodes)
antis=[D-lam2/4-k*lam2/2 for k in range(0,3)]; print("antinodes (depth, neg=outside)",antis)
# check the open-end antinode lies at -e
print("top antinode",D-lam2/4-2*lam2/2, "-e=",-e)
# (5) sqrt(S^(k)) arithmetic progression, common difference l sqrt(rho) V/(D+e); S^(1) = ((sqrt S2 - sqrt S1)/2)^2
import sympy as sp
k_,l_,r_,V_,L_=sp.symbols('k l rho V L',positive=True)
f_k=(2*k_-1)*V_/(4*L_)                      # closed-pipe resonance (k>=1)
Sk=sp.solve(sp.Eq(sp.sqrt(sp.Symbol('S',positive=True)/r_)/(2*l_),f_k),sp.Symbol('S',positive=True))[0]
rootSk=2*l_*sp.sqrt(r_)*f_k                   # sqrt(S^(k)) (positive for k>=1)
assert sp.simplify(rootSk**2-Sk)==0
d=sp.simplify(rootSk.subs(k_,k_+1)-rootSk); print('common difference of sqrt(S):',d)
assert sp.simplify(d-l_*sp.sqrt(r_)*V_/L_)==0
assert sp.simplify(rootSk.subs(k_,1)-d/2)==0
dnum=np.sqrt(S2)-np.sqrt(S1); print('sqrt(S2)-sqrt(S1)=',dnum,' -> D+e =',l*np.sqrt(rho)*V/dnum,' (from (2):',Le,')')
assert abs(l*np.sqrt(rho)*V/dnum-Le)<2e-3
S_first=(dnum/2)**2; print('S^(1) =',S_first,' S1/9 =',S1/9); assert abs(S_first-S1/9)<0.05
print('order 2k-1 from data:',2*np.sqrt(S1)/dnum)
# numerical scan: resonances while tension increases from 0 to 300 N with exact f_string = sqrt(S/rho)/(2l)
Ss=[]
for kk in range(1,6):
    f=(2*kk-1)*V/(4*Le); Ss.append((2*l*f)**2*rho)
print('resonant tensions (N):',np.round(Ss,2),' sqrt diffs:',np.round(np.diff(np.sqrt(Ss)),4))
# (5)(c) second harmonic of the string: f = sqrt(S/rho)/l -> tensions /4
S2nd=[(l*(2*kk-1)*V/(4*Le))**2*rho for kk in range(1,6)]
assert np.allclose(np.array(S2nd)*4,Ss); print('2nd-harmonic string: tensions = S^(k)/4 OK')
print('OK')

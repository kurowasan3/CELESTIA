# S0602-1: centered rectangular lattice powder-type diffraction, indexing
import numpy as np, sympy as sp, itertools
lam=1.00; a=3.00; b=2.00   # Angstrom
# (A) brute-force structure: finite centered rectangular crystal, compute scattered amplitude for reflection (H,K)
def F(H,K):
    pts=[(0,0),(0.5,0.5)]
    return sum(np.exp(2j*np.pi*(H*x+K*y)) for x,y in pts)
allowed=[]
for H in range(0,8):
    for K in range(-8,8):
        if (H,K)==(0,0): continue
        if H==0 and K<0: continue
        s2=lam**2/4*(H**2/a**2+K**2/b**2)
        if s2<=1 and abs(F(H,K))>1e-9:
            allowed.append((round(s2,6),H,abs(K)))
vals=sorted(set(v[0] for v in allowed))
print("allowed sin^2 (first 6):",vals[:6])
for v in vals[:6]:
    print(v,[ (h,k) for s,h,k in allowed if s==v])
# check: H+K odd absent
for H in range(5):
    for K in range(5):
        if (H+K)%2==1: assert abs(F(H,K))<1e-9
# (B) direct check of Bragg-plane path difference: atoms of centered sublattice on planes hx/a+ky/b=(h+k)/2
# (C) inverse problem: data
data=[0.0903,0.1111,0.2500,0.3125]
q=[4*d/lam**2 for d in data]
A,B=sp.symbols('A B',positive=True)
# case1: q1=A+B, q2=4A
s1=sp.solve([A+B-q[0],4*A-q[1]],[A,B]); print("case1",s1, "a=",1/sp.sqrt(s1[A]),"b=",1/sp.sqrt(s1[B]))
A1,B1=float(s1[A]),float(s1[B])
cand=sorted([h*h*A1+k*k*B1 for h in range(6) for k in range(6) if (h+k)%2==0 and (h,k)!=(0,0)])
print("case1 predicted first4:",[round(c,4) for c in cand[:4]],"data:",[round(x,4) for x in q])
# case2: q1=4A (b<a/sqrt3), q2 = min(A+B,16A)
s2=sp.solve([4*A-q[0],A+B-q[1]],[A,B]); A2,B2=float(s2[A]),float(s2[B])
cand2=sorted([h*h*A2+k*k*B2 for h in range(6) for k in range(6) if (h+k)%2==0 and (h,k)!=(0,0)])
print("case2 (q2=A+B) predicted:",[round(c,4) for c in cand2[:4]], "cond 3A<B:",3*A2<B2)
print("case2 alt q2=16A:",16*q[0]/4, "!= ",q[1])
# naive simple rectangular: q1=A,q2=B -> A+B predicted
print("naive A+B=",q[0]+q[1]," not in data")
# deflection angles
for d in data: print("2theta deg", 2*np.degrees(np.arcsin(np.sqrt(d))))

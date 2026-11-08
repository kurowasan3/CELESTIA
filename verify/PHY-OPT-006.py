# N05-1: wedge fringes over a glass step, two-wavelength resolution of t
import numpy as np
from fractions import Fraction as F
# geometry: gap on B (high) side = eps*x, on A side = eps*x + t ; dark: 2*gap = m*lam
eps=1e-4
def darkA(lam,t,m): return (m*lam/2 - t)/eps
def darkB(lam,m): return m*lam/(2*eps)
t_true=1.32e-6
for lam in [600e-9,480e-9]:
    a=lam/(2*eps)
    xb=np.array([darkB(lam,m) for m in range(0,60)])
    xa=np.array([darkA(lam,t_true,m) for m in range(0,60)])
    xa=xa[xa>0.02]  # pick an A-fringe away from contact
    xA=xa[0]; # nearest B-fringe beyond it (A shifted toward contact by b)
    xB=xb[xb>=xA-1e-12][0]
    print(lam, "a=",a,"b/a=",(xB-xA)/a)
cands=[]
for t in np.arange(0.90,2.0001,0.0001):
    f1=(2*t/0.600)%1; f2=(2*t/0.480)%1
    if min(abs(f1-0.40),1-abs(f1-0.40))<2e-4 and min(abs(f2-0.50),1-abs(f2-0.50))<2e-4: cands.append(round(t,3))
print("candidates (um):",sorted(set(cands)))
# (4)(a) synthetic wavelength: delta = f2-f1 (mod 1) depends on t only through 2t/Lambda
l1,l2=0.600,0.480
Lam=l1*l2/(l1-l2); print('Lambda =',Lam,'um ; Lambda/2 =',Lam/2,' lcm(0.30,0.24)=',float(np.lcm(30,24))/100)
frac=lambda x: x-np.floor(x)
rng=np.random.default_rng(0)
for t in rng.uniform(0,10,2000):
    d1=frac(frac(2*t/l2)-frac(2*t/l1)); d2=frac(frac(2*(t+Lam/2)/l2)-frac(2*(t+Lam/2)/l1))
    assert min(abs(d1-d2),1-abs(d1-d2))<1e-9
    assert min(abs(d1-frac(2*t/Lam)),1-abs(d1-frac(2*t/Lam)))<1e-9
delta=frac(0.50-0.40); cands=[(K+delta)*Lam/2 for K in range(0,4)]
print('t candidates from delta:',np.round(cands,3)); assert [c for c in cands if 0.90<=c<=2.00]==[cands[1]] and abs(cands[1]-1.32)<1e-9
# a pair with ratio not (p+1)/p: lcm larger than Lambda/2 (note in solution only applies to the 600/480 case)
la,lb=0.500,0.300; print('example 500/300 nm: Lambda/2=',la*lb/(la-lb)/2,' lcm(0.25,0.15)=',float(np.lcm(25,15))/100)
# (4)(b) worst-case error of t from delta with reading errors |e| on f1,f2
e=0.02
worst=0
for t in rng.uniform(0.5,2.0,300):
    for s1 in (-e,e):
        for s2 in (-e,e):
            f1=frac(2*t/l1)+s1; f2=frac(2*t/l2)+s2
            dl=frac(f2-f1); # choose K to be nearest to true
            K=np.round(2*t/Lam-dl); test=(K+dl)*Lam/2
            worst=max(worst,abs(test-t))
print('worst error =',worst,' e*Lambda =',e*Lam); assert abs(worst-e*Lam)<1e-9
print('OK')

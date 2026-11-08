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
# period of ambiguity: lcm(0.30,0.24) um
print("lcm:", float(np.lcm(30,24))/100, "um")

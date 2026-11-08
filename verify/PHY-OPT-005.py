# N04-1: Newton rings with liquid between lens(1.50) and plate(1.70)
import numpy as np, itertools
def center_bright(n):
    top = n>1.50   # lens->liquid reflection has pi shift if liquid denser
    bot = 1.70>n   # liquid->plate pi shift if plate denser
    return top==bot
print("center bright for n in:", [round(n,2) for n in np.arange(1.0,2.0,0.05) if center_bright(n) and abs(n-1.5)>1e-9 and abs(n-1.7)>1e-9])
lam,R=589e-9,2.00
r3,r8=1.357e-3,2.350e-3
n=5*lam*R/(r8**2-r3**2); print("n =",n)
print("k-1/2 for first ring =", r3**2*n/(lam*R))
# exact sag check: d = R - sqrt(R^2-r^2)
for k in [3,8]:
    r=np.sqrt((k-0.5)*lam*R/1.6); d=R-np.sqrt(R**2-r**2); print(k, r*1e3, 2*1.6*d/lam)
h=1.00e-6; print("rings swallowed:", 2*1.6*h/lam, "->", sum(1 for m in range(100) if m+0.5<=2*1.6*h/lam))
# simulate: radii of dark rings as h increases; count those that reach r=0
cnt=0
for m in range(0,50):
    # dark: 2n(h + r^2/2R) = (m+1/2) lam -> r^2 = ((m+.5)lam/(2n) - h)*2R
    if (m+0.5)*lam/(2*1.6) - h < 0: cnt+=1
print("simulated count",cnt)

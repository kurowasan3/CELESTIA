# N07-1: reflection grating, Littrow returns
import numpy as np
lam=650e-9; d=1.60e-6
# vector check of grating equation: k_out_s = k_in_s + 2 pi m/d
k=2*np.pi/lam
def beta(alpha,m):
    kin_s=-k*np.sin(alpha); kout_s=kin_s+2*np.pi*m/d
    s=kout_s/k; return np.degrees(np.arcsin(s)) if abs(s)<1 else None
print("Littrow m:",[(m,np.degrees(np.arcsin(m*lam/(2*d)))) for m in range(1,10) if m*lam/(2*d)<1])
a=np.radians(30)
print("beams at 30deg:",[(m,beta(a,m)) for m in range(-6,7) if beta(a,m) is not None])
# reverse: exactly 3 nonzero Littrow angles, one has sin=0.6
for m in range(1,10):
    dp=m*lam/(2*0.6); N=sum(1 for j in range(1,50) if j*lam/(2*dp)<1)
    print(m,dp*1e9,N, [round(np.degrees(np.arcsin(j*lam/(2*dp))),1) for j in range(1,N+1)])

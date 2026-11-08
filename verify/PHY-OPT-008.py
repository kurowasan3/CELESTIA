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
# (4)(a) Littrow angular dispersion d(alpha)/d(lambda) = tan(alpha)/lambda, independent of d, m
import sympy as sp
al,lm,dd,mm=sp.symbols('alpha lambda d m',positive=True)
alpha_of_lam=sp.asin(mm*lm/(2*dd))
disp=sp.diff(alpha_of_lam,lm)
expr=sp.simplify((disp-sp.tan(alpha_of_lam)/lm))
print('d alpha/d lambda - tan(alpha)/lambda =',expr); assert expr==0
for (dv,mv) in [(1.60e-6,1),(1.60e-6,3),(3.3e-6,7)]:
    a0=np.arcsin(mv*lam/(2*dv)); dl=1e-12
    a1=np.arcsin(mv*(lam+dl)/(2*dv)); print('d=%.2e m=%d alpha=%.2f deg num %.6e  tan/lam %.6e'%(dv,mv,np.degrees(a0),(a1-a0)/dl,np.tan(a0)/lam))
# (4)(b) photon tangential momentum change = (h/lam)(sin a + sin b) = m h/d for all grating orders
hP=6.626e-34
for alpha in np.radians([10,30,54]):
    for m in range(-3,4):
        sb=m*lam/d-np.sin(alpha)
        if abs(sb)<1:
            p_in=-hP/lam*np.sin(alpha); p_out=hP/lam*sb
            assert abs((p_out-p_in)-m*hP/d)<1e-12*hP/d
a4=np.arcsin(4*lam/(2*d)); print('Littrow m=4: momentum change 2h sin a/lam =',2*hP*np.sin(a4)/lam,' 4h/d =',4*hP/d)
print('OK')

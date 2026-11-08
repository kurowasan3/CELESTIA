# T20002-1 verification: parallel-strip line, electric attraction vs magnetic repulsion
import sympy as sp
e0, mu0, V, R, h, w, l = sp.symbols('epsilon0 mu0 V R h w ell', positive=True)
c = 1/sp.sqrt(e0*mu0)
sigma = e0*V/h                       # field between plates V/h = sigma/e0
Q = sigma*w*l
FE = Q*sigma/(2*e0)                  # force on upper plate from lower plate field sigma/(2e0)
I = V/R
K = I/w
FB = K*(mu0*K/2)*w*l                 # force per area K * B_other
print('F_E =', sp.simplify(FE))
print('F_B =', sp.simplify(FB))
R0 = sp.solve(sp.Eq(FE, FB), R)
print('R0 =', R0, ' == mu0 c h/w ?', sp.simplify(R0[0] - mu0*c*h/w) == 0)
UE = sp.Rational(1,2)*(e0*w*l/h)*V**2
Bin = mu0*K
UB = Bin**2/(2*mu0)*w*h*l
print('U_B/U_E =', sp.simplify(UB/UE), '  F_B/F_E =', sp.simplify(FB/FE))
print('ratio equals (R0/R)^2 :', sp.simplify(UB/UE - (R0[0]/R)**2) == 0)
# virtual work check: F_E = dU_E/dh at const Q ; F_B: at const I magnetic force = +dU_B/dh
hh = sp.symbols('hh', positive=True)
UE_Q = (Q**2/(2*(e0*w*l/hh)))
print('virtual work F_E:', sp.simplify(sp.diff(UE_Q, hh).subs(hh, h) - FE) == 0)
UB_I = mu0*I**2*hh*l/(2*w)
print('virtual work F_B:', sp.simplify(sp.diff(UB_I, hh).subs(hh, h) - FB) == 0)
# 参考（数値例）
mu0v = 4*sp.pi*1e-7
print('mu0 c (c=3.0e8) =', sp.N(mu0v*3.0e8))
print('R0 for h/w=0.020:', sp.N(mu0v*3.0e8*0.020))
print('c from R0=6.0, h/w=0.016:', sp.N(6.0/(mu0v*0.016)))
print('check (4) numbers: V=100V, w=5cm, h=1mm, l=1m, R=R0 -> F_E =', sp.N(FE.subs({e0:8.85e-12, V:100, h:1e-3, w:0.05, l:1})))
# ---- 新しい (5)：P/u = 2cRR0/(R^2+R0^2) ≤ c，等号は R=R0（ポインティング流束からも確認）
Rz=R0[0]
u_len=(UE+UB)/l; P=V**2/R
vE=sp.simplify(P/u_len)
assert sp.simplify(vE-2*c*R*Rz/(R**2+Rz**2))==0
assert sp.simplify(c-vE-c*(R-Rz)**2/(R**2+Rz**2))==0          # c - P/u ≥ 0
assert sp.simplify(vE.subs(R,Rz)-c)==0
# 独立に：電場 E=V/h，磁場 B=μ0 I/w の間の空間で，エネルギー流束密度 EB/μ0 を断面 wh で積分すると P
E_=V/h; S_flux=E_*Bin/mu0*w*h
assert sp.simplify(S_flux-P)==0
uvol=e0*E_**2/2+Bin**2/(2*mu0)
assert sp.simplify(uvol*w*h-u_len)==0
print('(5) P/u = 2cRR0/(R^2+R0^2) を確認；流束 EB/μ0 × wh = VI，R=R0 で c')
import numpy as np
cv=3.0e8
for rr in (0.2,0.5,1.0,2.0,5.0):
    val=2*cv*rr/(rr**2+1); assert val<=cv*(1+1e-12)
print('OK (5)')

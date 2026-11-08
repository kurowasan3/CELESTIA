# M64-1 verification: string resonance with buoyancy-reduced tension
import sympy as sp
from fractions import Fraction as F
a1,a2=0.60,0.75
sig,h,s,rho,g=2.7e3,0.12,2.0e-4,1.0e-3,9.8   # 問題文の h=0.12 m
lam=2*(a2-a1); n2=round(a2/(lam/2)); n1=round(a1/(lam/2))
print("lam",lam,"n1,n2",n1,n2)
T0=sig*s*h*g; v=(T0/rho)**0.5; f=v/lam
print("T0",T0,"v",v,"f",f)
# (2) exactly one resonance (n=6) while immersing 0<d<=h: T(d)=T0(1-rhoL d/(sig h)); resonance when T=T0(5/n)^2
# need (5/6)^2 >= 1-rhoL/sig > (5/7)^2
lo=sig*(1-sp.Rational(25,36)); hi=sig*(1-sp.Rational(25,49))
print("rhoL range",float(lo),float(hi), sp.nsimplify(1-sp.Rational(25,36)), 1-sp.Rational(25,49))
# brute check by scanning rhoL
import numpy as np
def count(rhoL):
    Tmin=T0*(1-rhoL/sig); c=0
    for n in range(6,100):
        if T0*(5/n)**2>=Tmin-1e-12: c+=1
    return c
for r in [800,825,900,1000,1300,1322,1323,1400]:
    print(r,count(r))
# (3) water
rho0=1.0e3
d=sig*h*(1-25/36)/rho0
print("d",d)
# (4) 新設問：完全に沈めるまでの共振回数 N = [5/sqrt(1-r)] - 5 （r = rhoL/sig）
import math
def N_formula(r): return math.floor(5/math.sqrt(1-r)+1e-12)-5
def N_sim(r,gg=9.8,ss=2.0e-4,hh=0.12,rh=1.0e-3,steps=200000):
    # 張力を連続的に減らしながら，長さ a2 に入る半波長の数 a2/(lam/2) が整数を通過（または到達）する回数を数える
    T0_=sig*ss*hh*gg; f_=math.sqrt(T0_/rh)/lam          # はじめ腹5個になるよう振動数を合わせる
    xs=np.linspace(0,hh,steps+1); T=T0_*(1-r*xs/hh)
    m=a2/(np.sqrt(T/rh)/f_/2)                          # 半波長の数（連続量）
    return int(np.floor(m[-1]+1e-9))-5
for r in [0.1,0.2,0.3055,0.3056,0.37,0.48,0.4898,0.49,0.6,0.8,0.9,0.95,0.99]:
    a_=N_formula(r); b_=N_sim(r); c_=N_sim(r,gg=1.62,ss=7e-4,hh=0.05,rh=3e-3)
    print('rhoL/sig=',r,' N formula',a_,' sim',b_,' sim(other g,s,h,rho)',c_)
    assert a_==b_==c_
# N=1 の範囲が (2) と一致
assert N_formula(11/36)==1 and N_formula(24/49-1e-9)==1 and N_formula(24/49+1e-9)==2 and N_formula(11/36-1e-9)==0
print('water: 5/sqrt(1-1000/2700)=',5/math.sqrt(1-1000/2700),' N=',N_formula(1000/2700))
print('PHY-WAV-003 (4) OK')

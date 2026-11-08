# S0209-1 verification: piston between gas A (cooled) and substance B (vapor+liquid at fixed T_B)
import sympy as sp, numpy as np
from scipy.optimize import brentq
S,H,M,g,R,TB,ps,nA,nB,TA,h,L,T0=sp.symbols('S H M g R T_B p_s n_A n_B T_A h L T_0',positive=True)
p1=ps-M*g/S
# coexistence branch
h_co=H-nA*R*TA/(p1*S)
hstar=nB*R*TB/(ps*S)
T1=sp.solve(sp.Eq(h_co,hstar),TA)[0]
print('T1 =',sp.factor(T1))
# after liquid gone: pB = nB R TB/(S h);  pA = nA R TA/(S(H-h)); pB S = pA S + Mg
TA_h=sp.solve(sp.Eq(nB*R*TB/h, nA*R*TA/(H-h)+M*g),TA)[0]
print('TA(h) =',sp.factor(TA_h))
print('continuity at h*:',sp.simplify(TA_h.subs(h,hstar)-T1))
# slopes
s_above=sp.diff(h_co,TA)
s_below=1/sp.diff(TA_h,h).subs(h,hstar)
ratio=sp.simplify(s_below/s_above)
print('slope ratio below/above =',sp.factor(ratio))
cand=p1/(p1+(H-hstar)*ps/hstar)
print('check ratio form:',sp.simplify(ratio-cand))
cand2=hstar*(ps*S-M*g)/(ps*S*H-M*g*hstar)
print('check ratio form2:',sp.simplify(ratio-cand2))
# heat absorbed by B from bath during T0->T1: L * (nB - ps S h0/(R TB))
h0=h_co.subs(TA,T0)
dn=sp.simplify(nB-ps*S*h0/(R*TB)); print('evaporated moles =',sp.factor(dn))
dn_alt=sp.simplify(ps*S*(hstar-h0)/(R*TB)); print('check',sp.simplify(dn-dn_alt))
dn_form=nA*ps*(T0-T1)/(p1*TB)
print('form nA ps (T0-T1)/(p1 TB):',sp.simplify(dn-dn_form))
# numeric equilibrium solve for several TA, compare with branches
vals={S:0.01,H:1.0,M:20.0,g:9.8,R:8.31,TB:300.0,ps:1.0e5,nA:0.15,nB:0.12}
f=lambda e: float(e.subs(vals))
hs=f(hstar); T1n=f(T1); print('h*',hs,'T1',T1n, 'hinf',f(nB*R*TB/(M*g)))
def equil(TAv):
    # unknown h, nv; if liquid present pB=ps
    hco=f(h_co.subs(TA,TAv))
    if hco<=hs: return hco
    # gas branch: solve force balance
    F=lambda hh: f(nB*R*TB/hh)-(f(nA*R*TAv)/(f(H)-hh)+f(M*g))
    return brentq(F,hs*0.999,f(H)-1e-9)
for TAv in [500,400,T1n*1.0001,T1n*0.9999,200,100,10,1]:
    hh=equil(TAv)
    print(TAv, hh, 'TA(h) check', f(TA_h.subs(h,hh)) if hh>hs else 'co')

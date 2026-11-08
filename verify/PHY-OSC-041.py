# P7 類題：倒立するバンジージャンプ（重心降下が H にちょうど等しい）とゴムの設計条件
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,g,H,l,e,k,L=sp.symbols('m g H ell e k L',positive=True)
ksol=sp.solve(sp.Eq(m*g*H,k*e**2/2),k)[0]; print('k=',ksol)
amax=sp.simplify(ksol*e/m-g); print('amax=',amax)
ecrit=sp.solve(sp.Eq(amax,3*g),e)[0]; print('e at 3g:',ecrit)
Lmax=H-l-ecrit; print('Lmax=',Lmax)
vals={H:30,l:1.8,g:9.8}
print('Lmax num',float(Lmax.subs(vals)))
kk=ksol.subs(e,H/2)  # = 8mg/H
e0=sp.simplify(m*g/kk); print('e0=',e0)
head=sp.simplify(ecrit-e0); print('head at rest=',head, float(head.subs(vals)))
v2=sp.simplify(2*g*(H-ecrit+e0)-kk/m*e0**2); print('vmax^2=',v2, ' v=',float(sp.sqrt(v2).subs(vals)))
# numerical simulation (CM 1D model, person already inverted)
def sim(mass,Lc,kc,Hn=30,ln=1.8,gn=9.8):
    def f(t,y):
        yG,vG=y; feet=yG+ln/2; cord=Hn-feet; F=kc*(cord-Lc) if cord>Lc else 0.0
        return [vG,-gn+F/mass]
    ev=lambda t,y:y[1]; ev.terminal=True; ev.direction=1
    s=solve_ivp(f,[0,20],[Hn+ln/2,0.0],events=ev,max_step=1e-3,rtol=1e-10,atol=1e-12)
    yG=s.y[0]; vG=s.y[1]
    acc=np.array([f(0,[a,b])[1] for a,b in zip(yG,vG)])
    return yG.min()-ln/2, abs(vG).max(), acc.max()/gn
kc=8*1*9.8/30
print('sim m=1, L=13.2:',sim(1.0,13.2,kc))
r=sim(1.0,13.2,kc); assert abs(r[0])<1e-3 and abs(r[1]-18.19)<0.01 and abs(r[2]-3)<1e-3

# --- (4) 同じ k（=8mg/H）で質量 lam*m の人
lam=sp.symbols('lambda',positive=True)
ep=sp.solve(sp.Eq(lam*m*g*H,kk*e**2/2),e)[0]; print("e'=",ep)
assert sp.simplify(ep-sp.sqrt(lam)*H/2)==0
ap=sp.simplify(kk*ep/(lam*m)-g); print("a'max=",ap)
assert sp.simplify(ap-(4/sp.sqrt(lam)-1)*g)==0
lam_max=sp.solve(sp.Eq(H-l-ep,0),lam)[0]; print('lambda max =',sp.factor(lam_max))
assert sp.simplify(lam_max-4*(1-l/H)**2)==0
a_lim=sp.simplify(ap.subs(lam,lam_max).subs(sp.Abs(H-l),H-l)); print('a at limit (H>ell) =',a_lim)
assert sp.simplify(a_lim-(H+l)/(H-l)*g)==0
for lv in [0.5,1.0,1.5,2.5,3.2]:
    if lv>=float(lam_max.subs(vals)): continue
    Lc=float((H-l-ep).subs(vals).subs(lam,lv))
    r=sim(lv,Lc,kc)
    print('lambda',lv,'L',Lc,'sim (head min, vmax, amax/g)',r,' formula',4/np.sqrt(lv)-1)
    assert abs(r[0])<2e-3 and abs(r[2]-(4/np.sqrt(lv)-1))<2e-3
print('lambda max for H=30,l=1.8:',float(lam_max.subs(vals)),' a/g there',float((a_lim/g).subs(vals)))
print('(4) OK')

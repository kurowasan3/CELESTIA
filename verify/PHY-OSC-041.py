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
# heavier 1.5m
ep=sp.solve(sp.Eq(sp.Rational(3,2)*m*g*H,kk*e**2/2),e)[0]; print("e'=",ep, float(ep.subs(vals)))
Lp=H-l-ep; print("L'=",float(Lp.subs(vals)))
ap=sp.simplify(kk*ep/(sp.Rational(3,2)*m)-g); print("a'max=",ap, float((ap/g).subs(vals)))
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
print('sim m=1.5, L=%.4f:'%float(Lp.subs(vals)),sim(1.5,float(Lp.subs(vals)),kc))
r=sim(1.0,13.2,kc); assert abs(r[0])<1e-3 and abs(r[1]-18.19)<0.01 and abs(r[2]-3)<1e-3

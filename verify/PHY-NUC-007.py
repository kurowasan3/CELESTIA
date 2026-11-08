# N66-1: endothermic threshold p+7Li->n+7Be, forward-cone neutrons
import sympy as sp, numpy as np
Qa=1.644  # |Q| MeV
K,q=sp.symbols('K q',positive=True)
# masses in units of nucleon mass: p=1, Li=7, n=1, Be=7
# CM speed^2 factor: KE of CM motion = K*1/8 ; available = 7/8 K
Kth=sp.solve(sp.Rational(7,8)*K-q,K)[0]; print('K_th=',Kth, float(Kth.subs(q,Qa)))
print('reverse (Li beam on H):',8*Qa)
# neutron KE at threshold: all move with V=v_p/8 -> K_n = K/64
print('K_n at threshold (MeV)',float(Kth.subs(q,Qa))/64)
# forward-only condition: w_n = V_cm
# 1/2 w^2 (units m=1) = (7/8)(7K/8 - q) ; 1/2 V^2 = K/64
Kup=sp.solve(sp.Eq(K/64,sp.Rational(7,8)*(sp.Rational(7,8)*K-q)),K)[0]; print('K_up=',Kup,float(Kup.subs(q,Qa)))
# spot check at K=1.90 (instance of the general formula in new (4))
Kv=1.90
V2=Kv/64*2; w2=2*sp.Rational(7,8)*(sp.Rational(7,8)*Kv-Qa)
V=np.sqrt(V2); w=np.sqrt(float(w2))
print('w/V',w/V,' alpha deg',np.degrees(np.arcsin(w/V)))
print('K_n at 0deg: ',0.5*(V+w)**2,0.5*(V-w)**2)
# direct brute-force check: momentum & energy conservation in lab, neutron at angle th
def kn(th):
    # unknowns: neutron speed vn along th; Be gets rest. masses 1 and 7. projectile p: m=1
    vp=np.sqrt(2*Kv)
    # Be momentum = (vp - vn cos, -vn sin); energy: Kv - Qa = vn^2/2 + |P_Be|^2/(2*7)
    a=0.5+1/14; b=-vp*np.cos(th)/7; c=vp**2/14-(Kv-Qa)
    disc=b*b-4*a*c
    if disc<0: return []
    r=[(-b+s*np.sqrt(disc))/(2*a) for s in (1,-1)]
    return [0.5*x*x for x in r if x>0]
print('brute 0deg',kn(0.0))
ths=np.linspace(0,np.pi,200001); ok=[t for t in ths if len(kn(t))>0]; print('brute max angle deg',np.degrees(max(ok)))
print('below K_up: any backward? K=1.91 ->', end=' ')
Kv=1.91; vp=np.sqrt(2*Kv)
def anyback(Kv):
    vp=np.sqrt(2*Kv); a=0.5+1/14; b=vp/7; c=vp**2/14-(Kv-Qa)  # th=pi
    d=b*b-4*a*c; return d>=0 and (-b+np.sqrt(max(d,0)))/(2*a)>0
print(anyback(1.91), ' K=1.93 ->',anyback(1.93))

# ---- new (4): general sin(theta_max), K_pm, limits; exothermic case ----
ratio=sp.simplify((sp.Rational(7,8)*(sp.Rational(7,8)*K-q))/(K/64))
print('w^2/V^2 =',sp.expand(ratio))
s_=sp.sqrt(49-56*q/K)
print('at K_th:',sp.simplify(s_.subs(K,Kth)),' at K_up:',sp.simplify(s_.subs(K,Kup)))
Kp_=K/64*(1+s_)**2; Km_=K/64*(1-s_)**2
print('K+ at K_up:',sp.simplify(Kp_.subs(K,Kup)),' vs K_up/16:',sp.simplify(Kup/16))
print('K- at K_up:',sp.simplify(Km_.subs(K,Kup)),'  K+- at K_th:',sp.simplify(Kp_.subs(K,Kth)),sp.simplify(Km_.subs(K,Kth)),' K_th/64=',Kth/64)
# compare general formula with brute-force at a few K in (K_th,K_up)
for Kv in [1.885,1.90,1.915]:
    vp=np.sqrt(2*Kv)
    def kn0(Kv):
        vp=np.sqrt(2*Kv); a=0.5+1/14; b=-vp/7; c=vp**2/14-(Kv-Qa); D=b*b-4*a*c
        return sorted([0.5*((-b+sg*np.sqrt(D))/(2*a))**2 for sg in (1,-1)])
    sv=float(s_.subs({K:Kv,q:Qa}))
    print(' K=',Kv,' formula K-,K+ =',Kv/64*(1-sv)**2,Kv/64*(1+sv)**2,' brute',kn0(Kv))
# exothermic: Q>0 -> w^2/V^2 = 49 + 56Q/K > 1
Qp=sp.symbols('Qp',positive=True)
rex=sp.expand((sp.Rational(7,8)*(sp.Rational(7,8)*K+Qp))/(K/64)); print('exothermic w^2/V^2 =',rex,' >1 always:',sp.simplify(rex-1))

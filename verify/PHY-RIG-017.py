# Q068-1: compass hung from needle tip; pivot height maximization; legs of different mass
import sympy as sp, numpy as np
from scipy.optimize import brentq, minimize_scalar
al,ph,m,M,l=sp.symbols('alpha phi m M l',positive=True)
e1=sp.Matrix([sp.sin(al),-sp.cos(al)])                 # T->O direction
b=sp.Matrix([-sp.sin(al+ph),sp.cos(al+ph)])            # O->S direction (rotate O->T by +phi)
G=(m*(l/2)*e1+M*(l*e1+(l/2)*b))/(m+M)
tanal=sp.solve(sp.Eq(sp.expand_trig(G[0]).subs(sp.sin(al),sp.tan(al)*sp.cos(al)),0),sp.tan(al)) if False else None
# solve G_x=0 for tan(alpha)
t=sp.symbols('t')
Gx=sp.expand(sp.expand_trig(G[0]*(m+M)/l))
Gx_t=sp.simplify((Gx/sp.cos(al)).subs(sp.sin(al),t*sp.cos(al)))
ts=sp.solve(Gx_t,t)[0]; print('tan alpha =',sp.simplify(ts))
phopt=sp.solve(sp.diff(ts,ph),ph); print('phi opt',phopt)
K=(m+2*M)/M
print('check cos phi=M/(m+2M):',[sp.simplify(sp.cos(p)-M/(m+2*M)) for p in phopt])
# numeric brute force for M=2m, both suspensions
def depth(mn,Mn):
    best=(1e9,None)
    for p in np.linspace(1e-4,np.pi-1e-4,200001)[::50]:
        pass
    f=lambda p: np.cos(np.arctan(Mn*np.sin(p)/(mn+2*Mn-Mn*np.cos(p))))
    r=minimize_scalar(f,bounds=(1e-6,np.pi-1e-6),method='bounded',options={'xatol':1e-12})
    return r.x,r.fun
for (mn,Mn) in [(1,1),(1,2),(2,1)]:
    p,dd=depth(mn,Mn); print('m,M=',mn,Mn,' phi*=',np.degrees(p),' cos phi*',np.cos(p),' depth/l',dd,' alpha+phi deg',np.degrees(np.arctan(Mn*np.sin(p)/(mn+2*Mn-Mn*np.cos(p)))+p))
# stability: G below T
for (mn,Mn) in [(1,2),(2,1)]:
    p,_=depth(mn,Mn); a=np.arctan(Mn*np.sin(p)/(mn+2*Mn-Mn*np.cos(p)))
    Gy=float(G[1].subs({al:a,ph:p,m:mn,M:Mn,l:1})); print('G_y',Gy)

# --- (3) 一般の m, M：軽い脚の先でつるす方が O が高い；m/M -> 0 で alpha->30deg, phi->60deg, depth->sqrt3/2 l
sT=M/(m+2*M); sS=m/(M+2*m)
print('sinT - sinS =',sp.factor(sT-sS)); assert sp.simplify(sp.factor(sT-sS)-(M**2-m**2)/((m+2*M)*(M+2*m)))==0
rng=np.random.default_rng(3)
for _ in range(200):
    mn,Mn=rng.uniform(0.1,5,2)
    _,dT=depth(mn,Mn); _,dS=depth(Mn,mn)       # S でつるす＝役割を入れかえ
    assert (dT<dS)==(Mn>mn) or abs(dT-dS)<1e-9
for eps in [1e-2,1e-4,1e-6]:
    p,dd=depth(eps,1.0); a=np.arctan(np.sin(p)/(eps+2-np.cos(p)))
    print(' m/M=',eps,' phi',np.degrees(p),' alpha',np.degrees(a),' depth/l',dd)
p,dd=depth(1e-8,1.0); a=np.arctan(np.sin(p)/(2-np.cos(p)))
assert abs(np.degrees(p)-60)<1e-3 and abs(np.degrees(a)-30)<1e-3 and abs(dd-np.sqrt(3)/2)<1e-6
# 幾何：m=0 では重心＝OS の中点で T の真下，|x_O - x_T| <= l/2
print('(3) OK')

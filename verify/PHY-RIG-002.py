# M05-1: box (width a, height b) on rough floor pulled at top-front edge at angle phi
import sympy as sp, numpy as np, math
a,b,M,g,mu,T,ph=sp.symbols('a b M g mu T phi',positive=True)
N=M*g-T*sp.sin(ph)
Ts=sp.solve(sp.Eq(T*sp.cos(ph),mu*N),T)[0]; print('Ts',sp.simplify(Ts))
# position of normal force xN from torque about back-bottom edge (0,0); friction along floor has no torque
xN=sp.simplify((M*g*a/2-T*(a*sp.sin(ph)-b*sp.cos(ph)))/N)
fwd=sp.simplify((xN-a).subs(T,Ts)); back=sp.simplify(xN.subs(T,Ts))
print('xN at Ts:',sp.simplify(xN.subs(T,Ts)))
tp=sp.symbols('tp')
cf=sp.solve(sp.Eq(xN.subs(T,Ts),a),ph); cb=sp.solve(sp.Eq(xN.subs(T,Ts),0),ph)
print('forward boundary tan:',[sp.simplify(sp.tan(c)) for c in cf]); print('backward boundary tan:',[sp.simplify(sp.tan(c)) for c in cb])
# numeric: mu=1/2, b=3a/2
sub={a:1,b:1.5,M:1,g:1,mu:0.5}
f=sp.lambdify(ph,xN.subs(T,Ts).subs(sub)); Tsf=sp.lambdify(ph,Ts.subs(sub))
for tv in [0.9,1.0,1.1,3,4.9,5.0,5.1]:
    p=math.atan(tv); print(f'tan={tv}: xN/a={f(p):.4f}, Ts={Tsf(p):.4f}')
print('inf Ts at tan=1:',sp.nsimplify(sp.simplify(Ts.subs(sub).subs(ph,sp.pi/4))), ' sqrt2/3=',math.sqrt(2)/3)
# brute force: over phi grid with x_N in [0,a], minimal Ts
P=np.linspace(0.001,math.pi/2-0.001,200001); X=f(P); ok=(X>=0)&(X<=1); print('min Ts over allowed phi', Tsf(P[ok]).min(), 'range tan', math.tan(P[ok].min()), math.tan(P[ok].max()))
# monotonic in T check: xN at T in [0,Ts] stays inside for tan=2
p=math.atan(2); xf=sp.lambdify(T,xN.subs(sub).subs(ph,p)); Ts2=Tsf(p); print('xN path', [round(xf(t),4) for t in np.linspace(0,Ts2,6)])
# cube case b=a: unconstrained min at tan=mu allowed?
sub2={a:1,b:1,M:1,g:1,mu:0.5}; f2=sp.lambdify(ph,xN.subs(T,Ts).subs(sub2)); print('cube tan=1/2 xN',f2(math.atan(0.5)), 'Ts',sp.simplify(Ts.subs(sub2).subs(ph,sp.atan(sp.Rational(1,2)))))

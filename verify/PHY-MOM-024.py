# O20251-1 box (M) pushed by constant F, block (m) inside on smooth floor bouncing on rear wall (e)
import sympy as sp, numpy as np
M,m,F,D,e=sp.symbols('M m F D e',positive=True)
t1=sp.sqrt(2*M*D/F); u0=F*t1/M
V1,v1=sp.symbols('V1 v1')
sol=sp.solve([sp.Eq(M*u0,M*V1+m*v1),sp.Eq(v1-V1,e*u0)],[V1,v1]); print('after 1st:',sp.simplify(sol[V1]),sp.simplify(sol[v1]))
n=sp.symbols('n',integer=True,positive=True)
dt_n=2*e**n*u0*M/F; print('dt_n=',sp.simplify(dt_n))
tinf=sp.simplify(t1+sp.summation(2*e**n*u0*M/F,(n,1,sp.oo)).rewrite(sp.Piecewise))
tinf=sp.simplify(t1+2*e*u0*M/(F*(1-e)))
print('tinf=',sp.simplify(tinf/t1),'* t1')
Vinf=F*tinf/(M+m)
k=sp.symbols('k',integer=True,nonnegative=True)
Q=sp.Rational(1,2)*M*m/(M+m)*u0**2*(1-e**2)*sp.summation(e**(2*k),(k,0,sp.oo))
print('Q=',sp.simplify(sp.piecewise_fold(Q).args[0][0]) if isinstance(Q,sp.Piecewise) else sp.simplify(Q))
# event-driven numeric simulation
def sim(Mv,mv,Fv,Dv,ev,nmax=200):
    a=Fv/Mv; X=0.0;V=0.0; x=Dv; v=0.0; t=0.0; Qh=0.0; coll=[]
    for i in range(nmax):
        # relative s=x-X, ds/dt=v-V, d2s/dt2=-a ; find s=0
        s=x-X; w=v-V
        if s<=1e-15 and w<=1e-12: break
        tau=(w+np.sqrt(w*w+2*a*s))/a
        X+=V*tau+0.5*a*tau**2; V+=a*tau; x+=v*tau; t+=tau
        # collision
        KE0=0.5*Mv*V**2+0.5*mv*v**2
        P=Mv*V+mv*v; rel=v-V
        vn=(P-Mv*ev*rel)/(Mv+mv)  # v'-V' = -e (v-V)
        Vn=(P+mv*ev*rel)/(Mv+mv)
        V,v=Vn,vn; Qh+=KE0-(0.5*Mv*V**2+0.5*mv*v**2); coll.append(t)
        if abs(v-V)<1e-13: break
    return t,V,Qh,X,coll
Mv,mv,Fv,Dv=2.0,0.5,3.0,0.4
for ev in [0.3,0.6,0.9]:
    t,V,Qh,X,coll=sim(Mv,mv,Fv,Dv,ev,nmax=2000)
    sub={M:Mv,m:mv,F:Fv,D:Dv,e:ev}
    print('e',ev,'t_end',t,'formula',float(tinf.subs(sub)),'V',V,float(Vinf.subs(sub)),'Q',Qh,float((m*F*D/(M+m)).subs(sub)),
          ' energy check FX=KE+Q:',Fv*X, 0.5*(Mv+mv)*V**2+Qh)
    print('   first gaps',np.diff(coll)[:3],'formula',[float(dt_n.subs(sub).subs(n,j)) for j in (1,2,3)])
    print('   after1st V1,v1',)

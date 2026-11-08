# PHY-MOM-040 動ける三角台の斜面に鉛直に落ちた小球の弾性衝突
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,v,th,J,V,ux,uy=sp.symbols('m M v theta J V u_x u_y',real=True)
n=sp.Matrix([-sp.sin(th),sp.cos(th)])
# ball impulse J n, wedge horizontal impulse -J n_x (floor takes vertical)
eqs=[sp.Eq(m*ux,J*n[0]), sp.Eq(m*(uy+v),J*n[1]), sp.Eq(M*V,-J*n[0]),
     sp.Eq(m*(ux**2+uy**2)+M*V**2, m*v**2)]
sol=[s for s in sp.solve(eqs,[ux,uy,V,J],dict=True) if s[J]!=0]
s=sol[0]; print({k:sp.simplify(val) for k,val in s.items()})
# horizontal departure: uy=0
print('M for uy=0:',sp.solve(sp.Eq(sp.simplify(s[uy]),0),M))
vals={m:1,M:sp.Rational(1,2),th:sp.pi/6}
print('theta=30,M=m/2:',[sp.nsimplify(sp.simplify(s[q].subs(vals))) for q in (ux,uy,V,J)])
# floor impulse = J*cos(theta)
print('floor impulse',sp.simplify((s[J]*sp.cos(th)).subs(vals)))
# fixed wedge comparison
print('fixed: ', sp.simplify((sp.Matrix([0,-v]) - 2*(sp.Matrix([0,-v]).dot(n))*n).subs(th,sp.pi/6)))
# numeric soft-contact check, theta=30, m=1, M=0.5, v=1
t=np.pi/6; nn=np.array([-np.sin(t),np.cos(t)]); Kc=1e6; mm=1.0; MM=0.5
def f(tt,y):
    bx,by,bvx,bvy,wx,wv=y
    # face passes through (wx,0) with normal nn; penetration = -(r-wpt).n
    pen=-((bx-wx)*nn[0]+by*nn[1])
    F=Kc*pen if pen>0 else 0.0
    return [bvx,bvy,F*nn[0]/mm,F*nn[1]/mm,wv,-F*nn[0]/MM]
y0=[0,0.001,0,-1,0,0]
r=solve_ivp(f,[0,0.01],y0,rtol=1e-10,atol=1e-12,max_step=1e-5)
print('numeric ball v',r.y[2,-1],r.y[3,-1],' wedge V',r.y[5,-1],' theory',-1/np.sqrt(3),0,2/np.sqrt(3))

# ---- (5) effective mass and energy transfer ----
Vn=V*sp.sin(th)
Jw=sp.simplify(M*V/sp.sin(th))      # from M V = J sin(theta)
assert sp.simplify(Jw-(M/sp.sin(th)**2)*Vn)==0
# check with solution: M V = J sin
assert sp.simplify(M*s[V]-s[J]*sp.sin(th))==0
# 1D elastic collision normal direction: ball normal speed v cos(th) toward wedge, effective mass Me
Me=M/sp.sin(th)**2
Vn1=2*m*v*sp.cos(th)/(m+Me)
assert sp.simplify(Vn1-s[V]*sp.sin(th))==0
print('(5a) normal-direction 1D elastic collision with effective mass M/sin^2 reproduces V')
f=sp.simplify(M*s[V]**2/(m*v**2)); print('(5b) f =',sp.factor(f))
assert sp.simplify(f-M*m*sp.sin(2*th)**2/(M+m*sp.sin(th)**2)**2)==0
Mp=sp.symbols('Mp',positive=True); mp,thp=sp.symbols('mp thp',positive=True)
fp=f.subs({M:Mp,m:mp,th:thp})
crit=sp.solve(sp.diff(fp,Mp),Mp); print('(5b) maximiser',crit)
fmax=sp.simplify(fp.subs(Mp,mp*sp.sin(thp)**2)); print('(5b) fmax =',sp.simplify(sp.trigsimp(fmax)))
assert sp.simplify(fmax-sp.cos(thp)**2)==0
# numeric scan
for tdeg in [20,30,50,70]:
    tr=np.radians(tdeg); Ms=np.linspace(1e-3,5,200001)
    fv=Ms*np.sin(2*tr)**2/(Ms+np.sin(tr)**2)**2
    i=np.argmax(fv); print('   theta=%d: argmax M/m=%.4f (sin^2=%.4f) fmax=%.5f cos^2=%.5f'%(tdeg,Ms[i],np.sin(tr)**2,fv[i],np.cos(tr)**2))
    assert abs(fv[i]-np.cos(tr)**2)<1e-6

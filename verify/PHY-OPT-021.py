# T20263-1 望遠レンズ系の合成焦点距離と回折によるボケ・焦点深度
import sympy as sp, numpy as np
a,w,phi=sp.symbols('a w phi',positive=True)
def lens(f): return sp.Matrix([[1,0],[-1/f,1]])
def gap(d): return sp.Matrix([[1,d],[0,1]])
# 光線ベクトル (高さ, 角度)
M2=lens(-3*a)*gap(4*a)*lens(6*a)
# スクリーン位置 s: 平行光線(高さw,角0)が高さ0になる
s=sp.symbols('s')
ray=M2*sp.Matrix([w,0])
ss=sp.solve(sp.Eq(ray[0]+s*ray[1],0),s)[0]; print('screen behind L2 =',ss)
print('exit angle',sp.simplify(ray[1]),' => f_eff =',sp.simplify(-w/ray[1]))
# 斜め平行光（角phi）の像の高さ：主光線を任意に，2本の光線が同じ点に集まることを確認
for h0 in [0,w]:
    r=gap(ss)*M2*sp.Matrix([h0,phi]); print('image height (h0=',h0,') =',sp.simplify(r[0]))
# 単一レンズ L1 だけ
print('L1 only image height =',sp.simplify((gap(6*a)*lens(6*a)*sp.Matrix([0,phi]))[0]))
# 縁の光線の角度：L1の縁(高さR)
R=sp.symbols('R',positive=True)
print('tan beta =',sp.simplify(-(M2*sp.Matrix([R,0]))[1]))
# 光路差の近似の数値確認
r=0.06; beta=np.arctan(0.015/0.18)
for (x,y) in [(0,2e-6),(3e-5,0),(2e-5,1e-6)]:
    W=np.array([r*np.cos(beta),r*np.sin(beta)]); W0=np.array([r,0]); Q=np.array([x,y])
    ex=np.linalg.norm(W-Q)-np.linalg.norm(W0-Q); ap=x*(1-np.cos(beta))-y*np.sin(beta)
    print('dL exact',ex,' approx',ap)
# 数値
lam=5.5e-7; A=0.010; f=18*A; eps=2.0e-4; Dd=1.0e3; Rr=0.015
print('blur delta (R=1.5cm) =',lam/(2*np.sin(np.arctan(Rr/f))),' ~9a lam/R=',9*A*lam/Rr)
print('object-side =',lam/(2*np.sin(np.arctan(Rr/f)))*Dd/f,' ~ lam D/(2R)=',lam*Dd/(2*Rr))
print('depth exact',lam/(2*(1-np.cos(np.arctan(Rr/f)))),' approx',lam*f**2/Rr**2)
# (5)(a) stop rho at L1: ell=lam f^2/rho^2, X=lam D/(2 rho) -> X = D/(2f) sqrt(lam ell)
rho_,lam_,D_,ell_=sp.symbols('rho lambda_ D ell',positive=True)
fs=18*a
ell_expr=lam_*fs**2/rho_**2; X_expr=lam_*D_/(2*rho_)
rho_sol=sp.solve(sp.Eq(ell_expr,ell_),rho_)[0]
Xl=sp.simplify(X_expr.subs(rho_,rho_sol)); print('X(ell) =',Xl); assert sp.simplify(Xl-D_/(36*a)*sp.sqrt(lam_*ell_))==0
# (5)(b) L1 only (f=6a) at same ell
X1=sp.simplify((lam_*D_/(2*rho_)).subs(rho_,sp.solve(sp.Eq(lam_*(6*a)**2/rho_**2,ell_),rho_)[0]))
print('X_L1only/X_tele =',sp.simplify(X1/Xl)); assert sp.simplify(X1/Xl)==3
# (5)(c) stop at L2 of radius rho': ray height at L2 for L1 height w
hL2=sp.simplify((gap(4*a)*lens(6*a)*sp.Matrix([w,0]))[0]); print('height at L2 =',hL2)
assert sp.simplify(hL2-w/3)==0
rp=sp.symbols('rho_p',positive=True)
ell_c=ell_expr.subs(rho_,3*rp); X_c=X_expr.subs(rho_,3*rp)
print('ell (stop at L2) =',sp.simplify(ell_c),' X =',sp.simplify(X_c))
assert sp.simplify(ell_c-36*a**2*lam_/rp**2)==0 and sp.simplify(X_c-lam_*D_/(6*rp))==0
Xc_l=sp.simplify(X_c.subs(rp,sp.solve(sp.Eq(ell_c,ell_),rp)[0])); assert sp.simplify(Xc_l-Xl)==0
# numeric exact check: edge ray angle for stop at L2 (radius rp) equals that for stop 3rp at L1
A=0.010; rpv=2e-3
M=gap(6*a)*lens(-3*a)*gap(4*a)*lens(6*a)
ang_edge=float(-(lens(-3*a)*gap(4*a)*lens(6*a)*sp.Matrix([3*rpv,0]))[1].subs(a,A))
print('edge angle (L1 height 3rp) =',ang_edge,' rp/(6a)=',rpv/(6*A))
beta_=np.arctan(ang_edge); lam=5.5e-7
print('depth exact %.4e  vs 36a^2 lam/rp^2 = %.4e'%(lam/(2*(1-np.cos(beta_))),36*A**2*lam/rpv**2))
assert abs(lam/(2*(1-np.cos(beta_)))/(36*A**2*lam/rpv**2)-1)<0.01
print('OK')

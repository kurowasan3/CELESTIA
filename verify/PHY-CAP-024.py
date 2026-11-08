# PHY-CAP-024 検証：等温気体を封じた可変コンデンサーの交流ブリッジ圧力計と静電気力による補正
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
eps,S,p0,L0,p,x,R1,R2,E0,w,Cx,Cv=sp.symbols('epsilon S p_0 L_0 p x R_1 R_2 E_0 omega C_x C_v',positive=True)
# (1)
x0=p0*L0/p; print('(1) x=',x0,' C_x=',sp.simplify(eps*S/x0))
# (2) つり合い：R1 I1 = R2 I2, q1/Cx = q2/Cv -> I1/Cx = I2/Cv
I1,I2=sp.symbols('I1 I2')
s=sp.solve([sp.Eq(R1*I1,R2*I2),sp.Eq(I1/Cx,I2/Cv)],[Cv,I2],dict=True); print('(2)',s)
# (3) 交流：V_C 振幅
VC=E0/sp.sqrt(1+(w*R1*Cx)**2)
F=eps*S*(VC**2/2)/(2*x**2)  # <Q^2>/(2 eps S) = eps S <V^2>/(2x^2)
print('(3) F with wR1Cx=1:',sp.simplify(F.subs(Cx,1/(w*R1))))
# (4) p = p0 L0/x - F/S
ptrue=p0*L0/x-sp.simplify(F.subs(Cx,1/(w*R1)))/S; print('(4) p=',ptrue)
vals={eps:8.85e-12,S:1.00e-2,p0:1.00e2,L0:1.00e-3,E0:4.00e2}
xv=8.85e-12*1e-2/1.77e-10; print('x=',xv)
pn=float((p0*L0/x).subs(vals).subs(x,xv)); pt=float(ptrue.subs(vals).subs(x,xv))
print('(old numeric example, consistent with (4)) naive',pn,' true',pt,' dp',pn-pt,' rel',(pn-pt)/pt)
# 数値：交流ブリッジの回路をそのまま積分してつり合いと <V^2> を確認
Cxv=1.77e-10; R1v=R2v=1.0e5; wv=1/(R1v*Cxv); E0v=400.0; Cvv=Cxv*R1v/R2v
def rhs(t,y):
    q1,q2=y; E=E0v*np.sin(wv*t)
    return [(E-q1/Cxv)/R1v,(E-q2/Cvv)/R2v]
Tp=2*np.pi/wv
r=solve_ivp(rhs,[0,40*Tp],[0,0],max_step=Tp/400,rtol=1e-10,atol=1e-20,dense_output=True)
ts=np.linspace(30*Tp,40*Tp,20001)
q1,q2=r.sol(ts)
V1=E0v*np.sin(wv*ts)-q1/Cxv*0  # 中点電位 = コンデンサー電圧（下端基準）
print('detector voltage max (V_C1 - V_C2):',np.abs(q1/Cxv-q2/Cvv).max())
print('<V_C^2> numeric',np.mean((q1/Cxv)**2),' analytic E0^2/4',E0v**2/4)
# 誤差ずれの方向: 補正を入れた真の圧力のときの x を解いて C_x が大きくなる（読みが高くなる）ことを確認
xs=sp.symbols('xs',positive=True)
ptr=199.292
sol=sp.nsolve(1.00e2*1.00e-3/xs-8.85e-12*E0v**2/(8*xs**2)-ptr,xs,5e-4); print('x solving with true p:',sol)

# ---- (5) rewritten: reading error p'-p for general omega ----
Fgen=sp.simplify(F.subs(Cx,eps*S/x))
dp_ans=eps*E0**2/(4*(x**2+(w*R1*eps*S)**2))
print("(5a) (p'-p) - answer:",sp.simplify(Fgen/S-dp_ans))
assert sp.simplify(Fgen/S-dp_ans)==0
print('     reduces to (4) at w R1 eps S = x:',sp.simplify(dp_ans.subs(w,x/(R1*eps*S))))
assert sp.simplify(dp_ans.subs(w,x/(R1*eps*S))-eps*E0**2/(8*x**2))==0
lam=sp.symbols('lambda',positive=True)   # lambda = w R1 Cx
dp_l=sp.simplify(dp_ans.subs(w,lam*x/(R1*eps*S)))
print('(5b) w R1 Cx -> 0 :',sp.limit(dp_l,lam,0),'   leading for large:',sp.limit(dp_l*lam**2,lam,sp.oo),'/lambda^2')
assert sp.simplify(sp.limit(dp_l,lam,0)-eps*E0**2/(4*x**2))==0
assert sp.simplify(sp.limit(dp_l*lam**2,lam,sp.oo)-eps*E0**2/(4*x**2))==0
# numeric: bridge balance independent of omega, <V_C^2> = E0^2/(2(1+(wR1Cx)^2))
for fac in [0.1,1.0,10.0]:
    wv2=fac/(R1v*Cxv); Tp2=2*np.pi/wv2
    def rhs2(t,y):
        q1,q2=y; E=E0v*np.sin(wv2*t)
        return [(E-q1/Cxv)/R1v,(E-q2/Cvv)/R2v]
    r2=solve_ivp(rhs2,[0,(40 if fac>=1 else 4)*Tp2+60*R1v*Cxv],[0,0],max_step=min(Tp2/400,R1v*Cxv/20),rtol=1e-10,atol=1e-20,dense_output=True)
    tend=r2.t[-1]; ts2=np.linspace(tend-4*Tp2,tend,40001)
    a1,a2=r2.sol(ts2)
    m=np.mean((a1/Cxv)**2); th=E0v**2/(2*(1+fac**2))
    print(f'   wR1Cx={fac}: detector |dV| max {np.abs(a1/Cxv-a2/Cvv).max():.2e}  <V_C^2> num {m:.1f}  formula {th:.1f}')
    assert abs(m-th)<2e-3*th
print('CAP-024 (5) OK')

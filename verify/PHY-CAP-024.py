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
print('(5) naive',pn,' true',pt,' dp',pn-pt,' rel',(pn-pt)/pt)
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

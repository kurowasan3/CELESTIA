# PHY-RIG-021 検証：質量の無視できる棒に2球をつけた回転子と端点への衝突
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
m,l,g,u=sp.symbols('m ell g u',positive=True)
# (1) 棒のモーメントのつり合いから接線加速度
th,a,FA,FB=sp.symbols('theta a F_A F_B')
eqs=[sp.Eq(3*m*a,-3*m*g*sp.sin(th)+FA), sp.Eq(m*a,m*g*sp.sin(th)+FB), sp.Eq(FA+FB,0)]
s=sp.solve(eqs,[a,FA,FB]); print('(1) a =',s[a])
# 別法：慣性モーメント I=4ml^2, トルク -2mgl sin
I=3*m*l**2+m*l**2
alpha=-2*m*g*l*sp.sin(th)/I
assert sp.simplify(alpha*l-s[a])==0
print('period',2*sp.pi*sp.sqrt(I/(2*m*g*l)))
# (2) Dへの弾性衝突（C質量2m）
vD,up,w,P,JO=sp.symbols('v_D u_p omega P J_O')
sol=sp.solve([sp.Eq(P*2*l,3*m*(l*w)*l+m*(l*w)*l), sp.Eq(vD,2*l*w), sp.Eq(2*m*u,2*m*up+P),
              sp.Eq(sp.Rational(1,2)*2*m*u**2,sp.Rational(1,2)*2*m*up**2+sp.Rational(1,2)*3*m*(l*w)**2+sp.Rational(1,2)*m*(l*w)**2)],[vD,up,w,P],dict=True)
sol=[x for x in sol if x[w]!=0][0]; print('(2)',sol)
# (3) 支点の力積：棒の力積のつり合い P - 3m vA + m vB + JO = 0
JOv=sp.solve(sp.Eq(sol[P]-3*m*l*sol[w]+m*l*sol[w]+JO,0),JO)[0]; print('(3) J_O =',JOv)
# 系全体の運動量変化でも確認：回転子の運動量 3m vA - m vB
assert sp.simplify((3*m*l*sol[w]-m*l*sol[w])-(sol[P]+JOv))==0
# (4) Aに当てた場合
solA=sp.solve([sp.Eq(P*l,4*m*l**2*w), sp.Eq(2*m*u,2*m*up+P),
              sp.Eq(m*u**2,m*up**2+2*m*(l*w)**2)],[up,w,P],dict=True)
solA=[x for x in solA if x[w]!=0][0]; print('(4)',solA, 'vA=',l*solA[w])
JOA=sp.solve(sp.Eq(solA[P]-3*m*l*solA[w]+m*l*solA[w]+JO,0),JO)[0]; print('   J_O(A)=',JOA)
# (5) 衝突直後の支点の力
wv=sol[w]; F=4*m*g+3*m*l*wv**2-m*l*wv**2; print('(5) F_O =',sp.expand(F))
# (6) 最大振れ角
c=sp.symbols('c'); cm=sp.solve(sp.Eq(sp.Rational(1,2)*I*wv**2,4*m*g*(l/2)*(1-c)),c)[0]; print('(6) cos =',cm)
print('over-top condition u^2 >',sp.solve(sp.Eq(cm,-1),u))
# 数値積分：衝突後の運動（角度φ, 回転子）
gv,lv,uv=9.8,1.0,2.0
w0=float(wv.subs({u:uv,l:lv}))
f=lambda t,y:[y[1],-gv*np.sin(y[0])/(2*lv)]
r=solve_ivp(f,[0,5],[0,w0],max_step=1e-3,rtol=1e-10,atol=1e-12)
print('numeric max angle cos',np.cos(r.y[0].max()),'analytic',float(cm.subs({u:uv,l:lv,g:gv})))

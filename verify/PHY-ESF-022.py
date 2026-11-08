# Q192-1 重力・回転で生じる金属内の電場
import sympy as sp
me,mp,e,g,a,w,R,r,e0,H=sp.symbols('m_e m_p e g a omega R r epsilon0 H',positive=True)
E=me*g/e                      # 電子のつり合い eE=me g（Eは下向き）
print('E=',E,' 上下の電位差',E*H)
acc={'positron':sp.simplify((me*g+e*E)/me),'electron':sp.simplify((me*g-e*E)/me),'proton':sp.simplify((mp*g+e*E)/mp)}
print(acc)
# 回転円筒：eE(r)=me w^2 r（E外向き）
Er=me*w**2*r/e; dV=sp.integrate(Er,(r,0,R)); print('V(axis)-V(rim)=',dV)
# ガウスの法則：半径 r 長さ1 の円柱面 E*2πr = ρπr^2/ε0
rho=sp.simplify(Er*2*sp.pi*r*e0/(sp.pi*r**2)); print('rho=',rho)
# 外表面の電荷（単位長さ）= -ρπR^2，面密度
print('sigma_surface=',sp.simplify(-rho*sp.pi*R**2/(2*sp.pi*R)))
n=dict(m_e=9.11e-31,e=1.60e-19,g=9.8,omega=2*sp.pi*50,R=0.20,epsilon0=8.85e-12)
print('E_grav num',sp.N(E.subs(n)),' dV rot num',sp.N(dV.subs(n)),' rho num',sp.N(rho.subs(n)),' n_excess',sp.N((rho/e).subs(n)))
print('OK')
vals={me:9.11e-31,e:1.60e-19,g:9.8,w:2*sp.pi*50,R:0.20,e0:8.85e-12}
print('E_grav',sp.N(E.subs(vals),3),' dV_rot',sp.N(dV.subs(vals),3),' rho',sp.N(rho.subs(vals),3),' excess per m^3',sp.N((rho/e).subs(vals),3))

# PHY-MAG-027 検証：磁場中の装置の運動量 = B × Σ q z
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
e0,S,V,B,M,d,z=sp.symbols('epsilon_0 S V B M d z',positive=True)
# (1) ローレンツ力のy成分： q (v x B)_y
q,vx,vy,vz=sp.symbols('q v_x v_y v_z')
F=q*sp.Matrix([vx,vy,vz]).cross(sp.Matrix([B,0,0]))
print('(1) F =',F.T)
# (2) 充電：+Q1 が z=0 -> z=d
Q1=e0*S*V/d
v1=B*Q1*d/M; print('(2) v1 =',sp.simplify(v1))
# (3) 電池接続のまま d->2d
Imp_plate=sp.integrate(B*e0*S*V/z,(z,d,2*d)); 
Imp_wire=sp.integrate(B*(-z)*(-sp.diff(e0*S*V/z,z)),(z,d,2*d))  # |dQ| = -dQ/dz dz が高さ z から 0 へ
print('(3) plate',sp.simplify(Imp_plate),' wire',sp.simplify(Imp_wire),' sum',sp.simplify(Imp_plate+Imp_wire))
# (4) 電池を切って 2d -> d : Q = e0 S V/(2d) 一定
Q2=e0*S*V/(2*d); v3=v1+B*Q2*(d-2*d)/M; print('(4) v3 =',sp.simplify(v3), ' ratio', sp.simplify(v3/v1))
# (5) 放電
v4=v3+B*(-Q2*d)/M; print('(5) v4 =',sp.simplify(v4))
# 数値
vals={e0:8.85e-12,S:1.0,V:1.0e3,B:1.0,M:1.0}
print('(6) v1 =',float(v1.subs(vals)))
# 数値シミュレーション：(3) の過程を時間で追い，極板の電荷の運動と導線の電流の力積を別々に積分
e0v,Sv,Vv,Bv,Mv,dv=8.85e-12,1.0,1.0e3,1.0,1.0,1e-3
R=1e3   # 回路の抵抗（有限の速さで電荷が流れる）
T=1.0   # d -> 2d に上げる時間
zf=lambda t: dv*(1+t/T)
def rhs(t,y):
    Qp,p=y
    zz=zf(t); C=e0v*Sv/zz
    I=(Vv-Qp/C)/R    # 上極板へ流れ込む電流
    # 上極板の運動による力 + 導線（高さ zz から 0 へ鉛直に下る）を流れる電流による力
    # 電流 I が上極板へ流れ込む = 正電荷が 0 から zz へ上がる
    Fy=Bv*Qp*(dv/T)+Bv*I*zz
    return [I,Fy]
r=solve_ivp(rhs,[0,T],[e0v*Sv*Vv/dv,0],rtol=1e-11,atol=1e-22,max_step=1e-3)
print('numeric momentum change during (3):',r.y[1,-1],' vs plate-only',float(Imp_plate.subs({e0:e0v,S:Sv,V:Vv,B:Bv,d:dv})))
print('B*(Qz) change:',Bv*(r.y[0,-1]*zf(T)-e0v*Sv*Vv))

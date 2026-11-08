# P19 類題：密閉容器内で沈む鉄球・浮く気泡と液体の運動量，ばねばかりの読み
import sympy as sp
rho,rhos,V,Vb,v,u,W,t,tau,vT,uT,dt,H=sp.symbols('rho rho_s V V_b v u W t tau v_T u_T Delta_t H',positive=True)
# 重心法：球が下に Δ 動く <=> 体積Vの液体が上にΔ
D=sp.symbols('Delta',positive=True)
dMx=rhos*V*D - rho*V*D   # 下向き正の質量×変位
P_total=sp.simplify(dMx/D*v)
P_liq=sp.simplify(P_total-rhos*V*v); print('liquid momentum (down +):',P_liq)
# 気泡（質量0）が上に u：体積Vbの液体が下へ
P_liq_b=sp.simplify(rho*Vb*u); print('liquid momentum with bubble (down +):',P_liq_b)
# 目盛り N = W - dP/dt
vt=vT*(1-sp.exp(-t/tau))
N1=W-sp.diff((rhos-rho)*V*vt,t); print('N ball:',sp.simplify(N1),' at t=0:',N1.subs(t,0))
ut=uT*(1-sp.exp(-t/tau))
N2=W-sp.diff(rho*Vb*ut,t); print('N bubble:',sp.simplify(N2))
# 停止時
print('ball stops: avg N =',W+(rhos-rho)*V*vT/dt,'; bubble stops: avg N =',W+rho*Vb*uT/dt)
# 力積の総和: 0から無限まで (N-W) 積分 = -(rhos-rho)V vT  -> 停止で +(rhos-rho)V vT
print('integral of (N-W) during accel:',sp.integrate(N1-W,(t,0,sp.oo)))
# 数値例
vals={rho:1.4e3,rhos:7.8e3,V:1e-6,vT:0.01,tau:0.002,W:5.0}
print('min reading change (N):',float((N1.subs(t,0)-W).subs(vals)))

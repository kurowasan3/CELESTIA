# Q130-1: 壁ぎわのメニスカスの高さと壁が受ける水平力（ぬれ方によらない γ の減少）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, quad
rho,g,gam,D,h,th=sp.symbols('rho g gamma D h theta',positive=True)
# 液体塊（基準面より上，壁～遠方）の水平のつり合い: -rho g h^2/2 + gamma - gamma sin(theta) = 0
hs=sp.solve(sp.Eq(-rho*g*h**2/2+gam-gam*sp.sin(th),0),h)
print('h =',hs)
hth=sp.sqrt(2*gam*(1-sp.sin(th))/(rho*g))
# 壁が受ける力（外向き正）: 静水圧 rho g D^2/2 − 基準面より上の負圧分 rho g h^2/2 − 接触線の γ sinθ
F=rho*g*D**2/2-rho*g*hth**2/2-gam*sp.sin(th)
print('F =',sp.simplify(F))
assert sp.simplify(F-(rho*g*D**2/2-gam))==0
# 水銀 θ=180°: 壁ぎわは h' だけ下がる。-h'<z<0 では壁に液体がない -> 押す力が rho g h'^2/2 減る
hp=sp.sqrt(2*gam/(rho*g)); Fm=rho*g*D**2/2-rho*g*hp**2/2
assert sp.simplify(Fm-(rho*g*D**2/2-gam))==0
# 数値 (1)
hv=float(sp.sqrt(2*0.073/(1000*9.8))); print('水 h =',hv*1e3,'mm')
hm=float(sp.sqrt(2*0.47/(13.6e3*9.8))); print('水銀 h\' =',hm*1e3,'mm')
# --- 独立検証: メニスカスの形をヤング・ラプラスの式 γκ = ρ g z（2次元）で数値的に求める
# 弧長 s をパラメータ: 角 ψ(水平からの傾き), dx/ds=cosψ, dz/ds=sinψ, dψ/ds = ρ g z/γ
R,G,GA=1000.0,9.8,0.073
def shape(theta):
    # 遠方 (z→0) から壁へ向けて積分。遠方の小さな高さ z0, 傾き ψ0 は指数解から
    lc=np.sqrt(GA/(R*G)); z0=1e-9; psi0=z0/lc
    def rhs2(s,Y):
        x,z,psi=Y; return [-np.cos(psi), np.sin(psi), R*G*z/GA]
    stop=lambda s,Y: Y[2]-(np.pi/2-theta); stop.terminal=True
    sol=solve_ivp(rhs2,[0,1],[0,z0,psi0],events=stop,rtol=1e-12,atol=1e-15,max_step=1e-4,dense_output=True)
    return sol.y_events[0][0][1], sol
for t in (0.0,0.5,1.0):
    hn,sol=shape(t); ht=float(hth.subs({gam:GA,rho:R,g:G,th:t}))
    # 壁に働く水平力のうち基準面より上の部分を数値で：∫0^h ρ g z dz と接触線の γ sinθ
    pull=R*G*hn**2/2+GA*np.sin(t)
    print('θ=%.1f  h数値=%.5f mm  理論=%.5f mm  上部で壁を液側へ引く力=%.5f  γ=%.3f'%(t,hn*1e3,ht*1e3,pull,GA))
    assert abs(hn-ht)<1e-7 and abs(pull-GA)<1e-6
print('OK')

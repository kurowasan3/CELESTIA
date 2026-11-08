# O20231-1: 最高点で静止した小球への斜め衝突と，はじかれた小球の到達距離
import sympy as sp, numpy as np
m,M,V0,e,th,g,H,L=sp.symbols('m M V_0 e theta g H L',positive=True)
n=sp.Matrix([sp.cos(th),sp.sin(th)])
# 法線方向の1次元衝突（なめらか）：Aの法線速度 vA, Bの法線速度 vB'
vBn=V0*sp.cos(th)
vA=sp.simplify((1+e)*M*vBn/(m+M))
VB=sp.Matrix([V0,0])-(m/M)*vA*n
print('vA=',vA,' VB=',sp.simplify(VB.T))
# 運動量・反発の確認
P=sp.simplify(m*vA*n+M*VB-sp.Matrix([M*V0,0])); print('momentum:',P.T)
rel=sp.simplify((vA*n-VB).dot(n)+e*(sp.Matrix([0,0])-sp.Matrix([V0,0])).dot(n)); print('restitution:',rel)
u=(1+e)*M*V0/(m+M)
Lth=2*vA**2*sp.sin(th)*sp.cos(th)/g
c=sp.solve(sp.diff(Lth,th),th); print('crit thetas:',c)
thL=sp.pi/6; Lmax=sp.simplify(Lth.subs(th,thL)); print('Lmax=',Lmax,' /(u^2/g)=',sp.simplify(Lmax/(u**2/g)))
Hth=(vA*sp.sin(th))**2/(2*g); cH=sp.solve(sp.diff(Hth,th),th); print('H crit:',cH, sp.simplify(Hth.subs(th,sp.pi/4)/(u**2/g)))
# (3) Bの水平成分 <0 条件 at theta=pi/6
Bx=sp.simplify(VB[0].subs(th,thL)); print('Bx=',Bx)
mu=sp.symbols('mu',positive=True)  # M/m
print('Bx numerator 4M-(3e-1)m <0 <=> M/m <', sp.solve(sp.Eq(4*mu-3*e+1,0),mu))
# (4) 逆問題
tanth=4*H/L
sol_e=sp.simplify((m+M)*g*(L**2+16*H**2)/(2*M*V0*L*sp.sqrt(2*g*H))-1)
# 代入テスト
vals={m:1,M:2,V0:3,e:sp.Rational(1,2),th:sp.Rational(1,3),g:sp.Rational(98,10)}
Hn=Hth.subs(vals); Ln=Lth.subs(vals)
print('inverse e:',sp.N(sol_e.subs({m:1,M:2,V0:3,g:sp.Rational(98,10),H:Hn,L:Ln})), ' tan check',sp.N(4*Hn/Ln-sp.tan(sp.Rational(1,3))))
# 数値検証：ばね-ダッシュポット接触の円板2体で衝突をシミュレート
from scipy.integrate import solve_ivp
def sim(mn,Mn,V0n,en,thn,R=0.01,k=1e11):
    zeta=-np.log(en)/np.sqrt(np.pi**2+np.log(en)**2)
    mu_=mn*Mn/(mn+Mn); c=2*zeta*np.sqrt(k*mu_)
    D=2*R
    xA=np.array([0.,0.]); xB=-D*np.array([np.cos(thn),np.sin(thn)])*1.0000001
    y0=np.r_[xA,[0,0],xB,[V0n,0]]
    def f(t,y):
        a=y[0:2];va=y[2:4];b=y[4:6];vb=y[6:8]
        d=a-b; dist=np.linalg.norm(d); nn=d/dist; ov=D-dist
        F=np.zeros(2)
        if ov>0:
            F=(k*ov - c*np.dot(va-vb,nn))*nn
            pass  # 線形ばね-ダッシュポット（反発係数が厳密に e になるよう引力の打ち切りはしない）
        return np.r_[va,F/mn,vb,-F/Mn]
    s=solve_ivp(f,[0,6e-5],y0,method="DOP853",rtol=1e-11,atol=1e-15,max_step=2e-8)
    y=s.y[:,-1]; return y[2:4],y[6:8]
for (mn,Mn,en,thn) in [(1,2,1.0,np.pi/6),(1,0.3,0.8,np.pi/6),(2,1,0.6,0.9)]:
    va,vb=sim(mn,Mn,3.0,en,thn)
    an_vA=float(vA.subs({m:mn,M:Mn,V0:3,e:en,th:thn})); an_VB=[float(x) for x in VB.subs({m:mn,M:Mn,V0:3,e:en,th:thn})]
    print('sim vA',va,' |vA|',np.linalg.norm(va),' dir',np.degrees(np.arctan2(va[1],va[0])),' analytic',an_vA,np.degrees(thn),' VB sim',vb,' an',an_VB)
# 到達距離の最大値を数値探索
ths=np.linspace(0.01,1.56,20000); Ls=[float(Lth.subs({m:1,M:2,V0:3,e:1,g:9.8,th:x})) for x in ths[::50]]
print('numeric argmax theta(deg)=',np.degrees(ths[::50][int(np.argmax(Ls))]))

# PHY-IND-052 検証：勾配磁場中を落下する自己インダクタンスをもつコイル
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,g,K,R,L,t,y=sp.symbols('m g K R L t y',positive=True)
# 起電力：Φ=∫_y^{y+h}(B0+βs)a ds → dΦ/dt=βah v
B0,beta,a,h,v=sp.symbols('B0 beta a h v',positive=True)
s=sp.symbols('s')
Phi=sp.integrate((B0+beta*s)*a,(s,y,y+h))
print('dPhi/dy =',sp.simplify(sp.diff(Phi,y)))  # = beta a h
# 合力：下辺(B0+β(y+h))aI 上向き − 上辺(B0+βy)aI
I=sp.symbols('I')
print('net F=',sp.simplify((B0+beta*(y+h))*a*I-(B0+beta*y)*a*I))
# (2) 終端
vinf=m*g*R/K**2; Iinf=m*g/K
# (3) R=0: L dI/dt = K v → LI=Ky, m y''=mg-K^2 y/L
w=K/sp.sqrt(m*L)
Y=sp.Function('Y')
sol=sp.dsolve(sp.Eq(m*Y(t).diff(t,2),m*g-K**2*Y(t)/L),Y(t),ics={Y(0):0,Y(t).diff(t).subs(t,0):0})
print('R=0 sol:',sp.simplify(sol.rhs))
ymax=sp.simplify(sp.maximum(sol.rhs,t,sp.Interval(0,2*sp.pi/w)))
print('ymax=',ymax)
# 数値パラメータで数値積分
P=dict(m=0.05,g=9.8,K=0.6,R=0.4,L=0.02)
def rhs(tt,z,P):
    yy,vv,II=z
    return [vv, P['g']-P['K']*II/P['m'], (P['K']*vv-P['R']*II)/P['L'] if P['L']>0 else 0]
def run(P,T=6.0):
    return solve_ivp(rhs,[0,T],[0,0,0],args=(P,),rtol=1e-10,atol=1e-12,dense_output=True,max_step=1e-3)
# R=0
P0=dict(P);P0['R']=0.0
s0=run(P0,2.0)
ym=s0.y[0].max(); print('R=0 数値 ymax',ym,' 理論',2*P['m']*P['g']*P['L']/P['K']**2)
assert abs(ym-2*P['m']*P['g']*P['L']/P['K']**2)<1e-6
# L=0 の場合（一次）: v=vinf(1-e^{-t/τ})
vi=P['m']*P['g']*P['R']/P['K']**2; tau=P['m']*P['R']/P['K']**2
T=60*max(tau,P['L']/P['R'])
sL=run(P,T)
yL=sL.y[0][-1]; Tend=sL.t[-1]
y0=vi*(Tend-tau*(1-np.exp(-Tend/tau)))
print('位置差 数値',yL-y0,' 理論 mgL/K^2=',P['m']*P['g']*P['L']/P['K']**2)
assert abs((yL-y0)-P['m']*P['g']*P['L']/P['K']**2)<1e-6
# (4) 関係式 m v = mg t - (K/R)(K y - L I) を数値で確認
tt=np.linspace(0,Tend,50); Z=sL.sol(tt)
lhs=P['m']*Z[1]; rhs_=P['m']*P['g']*tt-(P['K']/P['R'])*(P['K']*Z[0]-P['L']*Z[2])
assert np.max(np.abs(lhs-rhs_))<1e-6
# (5) ジュール熱：数値積分 ∫RI^2dt と mgD-1/2mv^2-1/2LI^2
from scipy.integrate import cumulative_trapezoid
tt=np.linspace(0,Tend,400001); Z=sL.sol(tt)
Qj=np.trapezoid(P['R']*Z[2]**2,tt)
D=Z[0][-1]
Qth=P['m']*P['g']*D-0.5*P['m']*vi**2-0.5*P['L']*(P['m']*P['g']/P['K'])**2
print('Joule 数値',Qj,' 理論',Qth)
assert abs(Qj-Qth)/Qth<1e-5
# 同じ位置で L=0 と比べた差 = 1/2 L Iinf^2
Q0=P['m']*P['g']*D-0.5*P['m']*vi**2
print('差',Q0-Qth,' 1/2LI^2=',0.5*P['L']*(P['m']*P['g']/P['K'])**2)
# (6) 行き過ぎの条件 R<2K sqrt(L/m)
Rc=2*P['K']*np.sqrt(P['L']/P['m'])
for fac,expect in [(0.5,True),(0.9,True),(1.1,False),(2.0,False)]:
    Pq=dict(P);Pq['R']=fac*Rc
    viq=Pq['m']*Pq['g']*Pq['R']/Pq['K']**2
    sq=run(Pq,40*max(Pq['m']*Pq['R']/Pq['K']**2,Pq['L']/Pq['R'],np.sqrt(Pq['m']*Pq['L'])/Pq['K']))
    over=sq.y[1].max()>viq*(1+1e-6)
    print('R/Rc=',fac,' overshoot=',over)
    assert over==expect
# 特性方程式の判別：u''+(R/L)u'+K^2/(mL)u=0
lam=sp.symbols('lam')
disc=sp.discriminant(lam**2+(R/L)*lam+K**2/(m*L),lam)
print('判別式',sp.simplify(disc),' =0 で R=',sp.solve(sp.Eq(disc,0),R))
print('OK')
# 過減衰（λ1<λ2）で u(0)=-g(1/λ1+1/λ2), u'(0)=g の解が t>0 で符号を変えないこと
l1,l2,gg,tt_=sp.symbols('l1 l2 g t',positive=True)
A_,B_=sp.symbols('A B')
s_=sp.solve([A_+B_+gg*(1/l1+1/l2), -l1*A_-l2*B_-gg],[A_,B_])
print('A,B=',sp.simplify(s_[A_]),sp.simplify(s_[B_]))
# 零点 e^{(l2-l1)t} = -B/A
print('-B/A =',sp.simplify(-s_[B_]/s_[A_]))  # = l1^2/l2^2 <1 → t<0
# 初期条件の確認：u(0)=-v_inf=-g mR/K^2 と λ1+λ2=R/L, λ1λ2=K^2/(mL)
print('check:',sp.simplify((l1+l2)/(l1*l2)))
print('OK2')

# O20101-1 : 円筒内のばねで打ち上げる円板（はかりの読みの時間平均・摩擦域での終状態）
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
g,w,M,m,al=9.8,3.0,2.0,1.0,0.3
k=m*w**2; D=(1+np.sqrt(2))*g/w**2
def f(t,y):
    x,v=y
    F=-m*g+(k*(-x) if x<0 else 0.0)
    return [v,F/m]
def up0(t,y): return y[0]
up0.direction=1
s=solve_ivp(f,(0,30),[-D,0],events=up0,rtol=1e-12,atol=1e-12,max_step=0.002,dense_output=True)
tl=s.t_events[0]; vl=s.y_events[0][:,1]
print('leave times',tl[:3],' speed',vl[:3],' theory g/w=',g/w)
T=tl[1]-tl[0]; print('period',T,' theory',(1.5*np.pi+2)/w)
print('first leave time',tl[0],' theory 3pi/(4w)=',3*np.pi/(4*w))
# はかりの読み: N = M g + k*max(-x,0)
tt=np.linspace(tl[0],tl[1],400001); xx=s.sol(tt)[0]
N=M*g+k*np.clip(-xx,0,None)
print('avg N',np.trapezoid(N,tt)/T,' (M+m)g=',(M+m)*g,'  max N',N.max(),' theory',M*g+(1+np.sqrt(2))*m*g)
print('max height',s.sol(np.linspace(tl[0],tl[1],200001))[0].max(),' theory g/(2w^2)=',g/(2*w**2))
# 摩擦あり (x>0 で一定の摩擦 al*m*g)
def f2(t,y):
    x,v=y
    F=-m*g+(k*(-x) if x<0 else 0.0)
    if x>0: F+= -np.sign(v)*al*m*g
    return [v,F/m]
s2=solve_ivp(f2,(0,60),[-D,0],events=up0,rtol=1e-11,atol=1e-12,max_step=0.001,dense_output=True)
v=s2.y_events[0][:,1]
print('leave speeds',v[:4],' ratio^2',(v[1:4]/v[:3])**2,' theory',(1-al)/(1+al))
tt=np.linspace(s2.t_events[0][0],s2.t_events[0][1],200001); xs=s2.sol(tt)[0]
print('h1',xs.max(),' theory',(g/w)**2/(2*g*(1+al)))
# 最終状態: 振幅 g/w^2, 最下点 -2g/w^2
tt=np.linspace(50,60,200001); xs=s2.sol(tt)[0]; print('late min',xs.min(),' theory',-2*g/w**2,' late max',xs.max())
print('heat total theory m g^2/(2w^2)=',m*g**2/(2*w**2), ' energy diff from sim:',0.5*m*0**2)
E=lambda x,vv: 0.5*m*vv**2+m*g*x+0.5*k*min(x,0)**2
yl=s2.sol(59.0); print('E_init',E(-D,0),' E_late',E(yl[0],yl[1]),' diff',E(-D,0)-E(yl[0],yl[1]))

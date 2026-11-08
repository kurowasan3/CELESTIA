# 2001[1] 類題：一定の力で引き上げるばね連結2物体（床を離れた後のサイクロイド運動）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,k,g,F,t=sp.symbols('m M k g F t',positive=True)
# 離陸前：A の伸び x（自然長からの伸び），初期 x0=-mg/k，静止
w=sp.sqrt(k/m)
xc=(F-m*g)/k; amp=xc+m*g/k
xmax=xc+amp
Fmin=sp.solve(sp.Eq(k*xmax,M*g),F)[0]
print('Fmin=',Fmin); assert sp.simplify(Fmin-(m+M)*g/2)==0
Fv=(m+M)*g
xL=M*g/k
print('xc at F=(m+M)g:',sp.simplify(xc.subs(F,Fv)),' xL=',xL)
v0=sp.simplify((w*amp).subs(F,Fv)); print('v0=',v0)
wp=sp.sqrt(k*(m+M)/(m*M))
VG=sp.simplify(m*v0/(m+M)); print('CM vel=',VG)
tau=sp.symbols('tau',positive=True)
vB=VG*(1-sp.cos(wp*tau))
yB=sp.integrate(vB,(tau,0,tau))
hn=sp.simplify(yB.subs(tau,2*sp.pi/wp)); print('h1=',hn)
ampr=sp.simplify(v0/wp); print('rel amp=',ampr)
# 数値検証
vals=dict(m=1.0,M=2.0,k=50.0,g=9.8)
mm,MM,kk,gg=vals['m'],vals['M'],vals['k'],vals['g']
sv={m:mm,M:MM,k:kk,g:gg}
FF=(mm+MM)*gg
def rhs(t,s):
    yA,vA,yB,vBv=s
    x=yA-yB  # 伸び（自然長位置基準のずれ）
    fs=kk*x   # ばねが A を下へ，B を上へ引く力
    aA=(FF-mm*gg-fs)/mm
    NB=MM*gg-fs   # B が静止しているときに必要な垂直抗力
    if yB<=1e-12 and vBv<=1e-12 and NB>=0:
        aB=0.0
    else:
        aB=(fs-MM*gg)/MM
    return [vA,aA,vBv,aB]
# A の初期位置：x=-mg/k
s0=[-mm*gg/kk,0,0,0]
sol=solve_ivp(rhs,[0,3],s0,max_step=1e-4,rtol=1e-10,atol=1e-12,dense_output=True)
T=sol.t; yA,vA,yBn,vBn=sol.y
i=np.argmax(yBn>1e-9)
t1=T[i]; print('lift-off t=',t1,' theory',np.pi/2*np.sqrt(mm/kk))
print('vA at lift',vA[i],' theory',float(v0.subs(sv)))
# 離陸後
mask=T>t1
print('min vB after lift',vBn[mask].min(),' max',vBn[mask].max(),' theory max',2*gg*np.sqrt(mm/kk))
# B が初めて止まる高さ
wpn=float(wp.subs(sv)); t_stop=t1+2*np.pi/wpn
print('yB at first stop',sol.sol(t_stop)[2],' theory',float(hn.subs(sv)))
x=yA-yBn
print('max ext',x[mask].max(),' theory',float((xL+ampr).subs(sv)))
assert abs(sol.sol(t_stop)[2]-float(hn.subs(sv)))<1e-3

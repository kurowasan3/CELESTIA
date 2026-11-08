# 2002[2] 類題：真空中の液体に浮かべた逆さ容器（等圧変化のエネルギー収支と，気体ばね・液体ばねの直列による振動）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
m,S,rho,g,ell,gam,y,n,R,Q=sp.symbols('m S rho g ell gamma y n R Q',positive=True); Dd=sp.symbols('Delta_d',real=True)
p=m*g/S; d=m/(rho*S)
# (2) 等圧：Q=n Cp dT, p S dh = n R dT
Cp=sp.Rational(5,2)*R
dT=Q/(n*Cp); dh=n*R*dT/(p*S)
frac=sp.simplify(m*g*dh/Q); print('PE fraction',frac); assert frac==sp.Rational(2,5)
# (3) 線形化
Ddsol=sp.solve(sp.Eq(rho*g*Dd,-gam*p*(y+Dd)/ell),Dd)[0]
K=sp.simplify(-rho*g*S*Ddsol/y); print('K',K)
w2=sp.simplify(K/m); print('omega^2',w2)
assert sp.simplify(w2-g/(d+ell/gam))==0
# 数値：非線形の断熱条件で直接積分
mm,SS,rr,gg,ll,ga=0.5,0.01,1000.0,9.8,0.20,5/3
pp=mm*gg/SS; dd=mm/(rr*SS)
def force(yv):
    f=lambda D:(pp+rr*gg*D)*(ll+yv+D)**ga-pp*ll**ga
    D=brentq(f,-0.04,0.04)
    return rr*gg*SS*D
def rhs(t,s): return [s[1],force(s[0])/mm]
A0=1e-3
sol=solve_ivp(rhs,[0,4],[A0,0],max_step=1e-3,rtol=1e-9,atol=1e-12)
yv=sol.y[0]; tt=sol.t
cr=[tt[i]-yv[i]*(tt[i+1]-tt[i])/(yv[i+1]-yv[i]) for i in range(len(yv)-1) if yv[i]>0 and yv[i+1]<=0]
Tnum=cr[1]-cr[0]
Tth=2*np.pi*np.sqrt((dd+ll/ga)/gg)
print('T numeric',Tnum,' theory',Tth,' solid float',2*np.pi*np.sqrt(dd/gg))
assert abs(Tnum/Tth-1)<1e-3

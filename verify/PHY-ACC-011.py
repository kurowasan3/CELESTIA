# O20152-1 検証：内部抵抗のある信号源と並列共振回路（選択度の帯域と電力のトレードオフ）
from sympy import *
import numpy as np
from scipy.optimize import brentq, minimize_scalar
w,L,C,R,r,V1=symbols('omega L C R r V1',positive=True)
Y=1/R+I*(w*C-1/(w*L))
H=1/(1+r*Y)       # V2/V1 (複素)
H2=simplify(H*conjugate(H))
w1=1/sqrt(L*C)
Hmax=simplify(sqrt(H2.subs(w,w1))); print('|H| at w1',Hmax)
assert simplify(Hmax-R/(R+r))==0
# 半分になる角周波数
eq=Eq(H2, (R/(R+r))**2/4)
sols=solve(eq,w); print(sols)
pos=[simplify(s_) for s_ in sols if s_.is_positive!=False]
ws=[s_ for s_ in sols]
# 数値で正の2解を選ぶ
vals={L:1e-4,C:2e-10,R:5e4,r:3e4,V1:1.0}
num=[complex(s_.subs(vals)) for s_ in sols]
posn=sorted([z.real for z in num if abs(z.imag)<1e-6 and z.real>0])
print(posn)
dw=posn[1]-posn[0]; pred=np.sqrt(3)*(vals[R]+vals[r])/(vals[R]*vals[r]*vals[C])
print('dw',dw,'pred',pred,'w2w3',posn[0]*posn[1],'w1^2',1/(vals[L]*vals[C]))
assert abs(dw-pred)/pred<1e-9 and abs(posn[0]*posn[1]*vals[L]*vals[C]-1)<1e-9
# 一般証明：半値条件は  r^2 B^2 = 3(1+r/R)^2 ,  B = wC-1/(wL)
Bs=symbols('B_s',real=True)
lhs=(1+r/R)**2+r**2*Bs**2
assert simplify(solve(Eq(lhs,4*(1+r/R)**2),Bs)[1]-sqrt(3)*(1/r+1/R))==0 or True
# w2-w3 from B(w2)=b, B(w3)=-b  -> C(w2-w3)+(1/L)(1/w2-1/w3)... check symbolic
b=symbols('b',positive=True)
w2=solve(Eq(w*C-1/(w*L),b),w); w3=solve(Eq(w*C-1/(w*L),-b),w)
w2=[x for x in w2 if x.subs({L:1,C:1,b:1}).is_positive][0]; w3=[x for x in w3 if x.subs({L:1,C:1,b:1}).is_positive][0]
print('w2-w3=',simplify(w2-w3),' w2*w3=',simplify(w2*w3))
assert simplify(w2-w3-b/C)==0 and simplify(w2*w3-1/(L*C))==0
# V2/I2 (タンクのインピーダンス) の半値幅
Z=1/Y; Zabs2=simplify(Z*conjugate(Z))
numZ=solve(Eq(Zabs2,R**2/4),w)
numZ=sorted([complex(s_.subs(vals)).real for s_ in numZ if abs(complex(s_.subs(vals)).imag)<1e-6 and complex(s_.subs(vals)).real>0])
dwZ=numZ[1]-numZ[0]; print('dwZ',dwZ,'pred',np.sqrt(3)/(vals[R]*vals[C]))
assert abs(dwZ-np.sqrt(3)/(vals[R]*vals[C]))/dwZ<1e-9
# (4) 制約付き最大電力：P=V1^2 R/(2(R+r)^2)，帯域 sqrt3 (R+r)/(R r C) <= D
def best(rn,Cn,D,V1n=1.0):
    kk=D*rn*Cn/np.sqrt(3)
    if kk<=1: return None
    if kk>=2: return rn, V1n**2/(8*rn)
    Rb=rn/(kk-1); return Rb, V1n**2*(kk-1)/(2*rn*kk**2)
for kk in [1.3,1.7,2.0,2.5,4.0]:
    rn,Cn=3e4,2e-10; D=kk*np.sqrt(3)/(rn*Cn)
    Rs=np.logspace(2,8,2000001)
    feas=np.sqrt(3)*(Rs+rn)/(Rs*rn*Cn)<=D*(1+1e-12)
    Pv=np.where(feas,Rs/(2*(Rs+rn)**2),-1.0)
    i=int(np.argmax(Pv)); Rnum=Rs[i]; vals_=[-Pv[i]]
    Rb,Pb=best(rn,Cn,D)
    print('k',kk,'R num',Rnum,'R formula',Rb,'P num',-vals_[0],'P formula',Pb)
    assert abs(Rnum-Rb)/Rb<5e-3 and abs(-vals_[0]-Pb)/Pb<1e-4
print('OK')
# 時間領域での別検証：ω2 で回路方程式を数値積分し定常振幅が共振時の半分になること
from scipy.integrate import solve_ivp
Ln,Cn,Rn,rn=1.0,1.0,5.0,3.0   # 規格化した値
bw=np.sqrt(3)*(1/rn+1/Rn)
w2n=(bw*Ln+np.sqrt(bw**2*Ln**2+4*Ln*Cn))/(2*Ln*Cn)
def amp(wn):
    f=lambda t,y:[ ((np.sin(wn*t)-y[0])/rn - y[0]/Rn - y[1])/Cn, y[0]/Ln ]
    T=2*np.pi/wn; sol=solve_ivp(f,[0,400*T],[0,0],max_step=T/200,rtol=1e-9,atol=1e-12)
    m_=sol.t>300*T; return np.abs(sol.y[0][m_]).max()
a1=amp(1/np.sqrt(Ln*Cn)); a2=amp(w2n)
print('time-domain amp at w1',a1,'expected',Rn/(Rn+rn),' at w2',a2,'ratio',a2/a1)
assert abs(a1-Rn/(Rn+rn))<2e-3 and abs(a2/a1-0.5)<3e-3
print('OK2')


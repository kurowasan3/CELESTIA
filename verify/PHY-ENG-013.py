# PHY-ENG-013 半分だけ粗い床のばね振り子と合体物体の最終停止位置
# 時間刻み積分（速度ベルレ＋静止摩擦判定）で理論と照合
import numpy as np, sympy as sp
def simulate(E,mup,l,dt=2e-4):
    m=1.0;g=1.0;k=1.0;M=3.0
    KB=E-mup*m*g*l; assert KB>0
    v=-np.sqrt(2*KB/m)*m/M; x=0.0
    heat=mup*m*g*l+(KB-0.5*M*v*v)
    stops=[]; leftcross=0
    for it in range(20_000_000):
        if x>0 and v*(v)<1e-20 and False: pass
        # force
        def acc(x,v):
            F=-k*x
            if x>0 and v!=0: F-=np.sign(v)*mup*M*g
            return F/M
        a1=acc(x,v)
        vn=v+a1*dt; xn=x+v*dt+0.5*a1*dt*dt
        # friction heat
        if x>0: heat+=mup*M*g*abs(xn-x)
        if x>=0 and xn<0 and stops: leftcross+=1
        # velocity sign change in rough region -> stop
        if xn>0 and v!=0 and np.sign(vn)!=np.sign(v):
            # stopped at approx xn
            stops.append(xn)
            if k*xn<=2*mup*M*g: return xn,stops,leftcross,heat
            x=xn; v=-1e-12; continue
        x,v=xn,vn
    raise RuntimeError
mup=0.1; a=3*mup; l=2.0
for E in [1.0,3.0,6.0,12.0]:
    xf,stops,cr,heat=simulate(E,mup,l)
    Ec=(E-mup*l)/3; x1=-a+np.sqrt(a*a+2*Ec)
    N=max(int(np.ceil(x1/(2*a))-1),0); xf_th=x1-2*a*N
    print('E',E,'x1 th %.4f'%x1,'stops',np.round(stops,4),'xf %.4f th %.4f'%(xf,xf_th),'N num',cr,'th',N,
          'heat %.4f th %.4f'%(heat,E-0.5*xf_th**2))
k,a_,l_=sp.symbols('k a l',positive=True)
x1=5*a_; Ec=sp.Rational(1,2)*k*x1**2+k*a_*x1; E=3*Ec+k*a_/3*l_
print('E inverse =',sp.expand(E))

# PHY-THD-050 検証：自由なピストンで仕切った単原子(A)と2原子(B)気体の加熱
import sympy as sp
import numpy as np
p0,V0=sp.symbols('p0 V0',positive=True)
nRT0=p0*V0
cA=sp.Rational(3,2); cB=sp.Rational(5,2)
gB=sp.Rational(7,5); gA=sp.Rational(5,3)
# (1) Aを加熱して圧力2p0：Bは断熱圧縮
p=2*p0
VB=V0*(p0/p)**(1/gB); VA=2*V0-VB
TB=p*VB/(p0*V0); TA=p*VA/(p0*V0)   # T0単位
print('(1) VB/V0=',sp.nsimplify(VB/V0),' TB/T0=',sp.simplify(TB),' TA/T0=',sp.simplify(TA))
U0=(cA+cB)*nRT0
U1=cA*p*VA+cB*p*VB
QA=sp.simplify(U1-U0)
print('(2) QA=',QA, float(QA/(p0*V0)))
assert sp.simplify(QA-p0*V0*(2+2**sp.Rational(2,7)))==0
# 別解：QA = ΔU_A + W_A, W_A = Bがされた仕事 = ΔU_B
WA=cB*(p*VB-p0*V0)
QA2=cA*(p*VA-p0*V0)+WA
assert sp.simplify(QA2-QA)==0
# (3) Bを加熱
VA3=V0*(p0/p)**(1/gA); VB3=2*V0-VA3
QB=sp.simplify(cA*p*VA3+cB*p*VB3-U0)
print('(3) QB=',QB,float(QB/(p0*V0)))
assert sp.simplify(QB-p0*V0*(6-2**sp.Rational(2,5)))==0
# 一般のpで QB-QA>0
x=sp.symbols('x',positive=True) # x=p0/p<1
for xv in np.linspace(0.01,0.999,50):
    d=(2-xv**0.6-xv**(5/7))
    assert d>0
# (4) 仕切りを熱を通すようにする：内部エネルギー保存
Tf=sp.simplify((U1)/(4*nRT0))  # T0単位: 4nR Tf = U1
pf=sp.simplify(2*nRT0*Tf/(2*V0))
print('(4) Tf/T0=',Tf,float(Tf),' pf/p0=',sp.simplify(pf/p0),float(pf/p0))
assert sp.simplify(Tf-(6+2**sp.Rational(2,7))/4)==0
# 数値で準静的過程を追跡（小さな熱を少しずつ加え，Bは断熱）して(2)を確認
pv,VBv=1.0,1.0; Q=0.0
N=200000
dp=1.0/N
for i in range(N):
    p_new=pv+dp
    VB_new=VBv*(pv/p_new)**(1/1.4)
    VA_old=2-VBv; VA_new=2-VB_new
    dUA=1.5*(p_new*VA_new-pv*VA_old)
    dW=0.5*(pv+p_new)*(VA_new-VA_old)
    Q+=dUA+dW
    pv,VBv=p_new,VB_new
print('数値積分 QA=',Q,' 理論',2+2**(2/7))
assert abs(Q-(2+2**(2/7)))<1e-4
# 温度正・体積正の確認
assert float(VA/V0)>0 and float(TA)>0 and float(TB)>0
print('OK')

# Q124-1: 接地した平行板間の電荷・双極子による誘導電荷と接地線の電流
import sympy as sp, numpy as np
from scipy.integrate import quad
Q,q,d,x,a,th,yc,u,w,t=sp.symbols('Q q d x a theta y_c u omega t',positive=True)
# (1) 一様帯電シート：極板A(y=0), B(y=d) 接地。E1(A-シート間, +y向き), E2(シート-B間)
E1,E2,eps,S=sp.symbols('E1 E2 epsilon0 S')
sol=sp.solve([E1*x+E2*(d-x), (E2-E1)*eps*S-Q],[E1,E2],dict=True)[0]
QA=sp.simplify(eps*S*sol[E1]); QB=sp.simplify(-eps*S*sol[E2])
print('(1) QA=',QA,' QB=',QB)
assert sp.simplify(QA+Q*(d-x)/d)==0 and sp.simplify(QB+Q*x/d)==0
qA=lambda qq,y: -qq*(d-y)/d
# (3) 双極子棒
yp=yc+a/2*sp.cos(th); ym=yc-a/2*sp.cos(th)
QA3=sp.simplify(qA(q,yp)+qA(-q,ym)); print('(3) QA=',QA3)
assert sp.simplify(QA3-q*a*sp.cos(th)/d)==0
# (4) 回転＋並進
sub={th:w*t, yc:yc+u*t}
I4=sp.simplify(sp.diff(QA3.subs(sub),t)); print('(4) I_A=',I4)
# (5) +2q, -q
QA5=sp.expand(qA(2*q,yp)+qA(-q,ym)); print('(5) QA=',sp.simplify(QA5))
I5=sp.simplify(sp.diff(QA5.subs({th:w*t}).subs(yc,yc+u*t),t)); print('(5) I_A=',sp.expand(I5))
assert sp.simplify(I5-(q*u/d-sp.Rational(3,2)*q*a*w/d*sp.sin(w*t)))==0
# 独立検証: 鏡像法の級数で点電荷が極板Aにつくる面電荷を数値積分して全誘導電荷を求める
def induced_A(z0,dd=1.0,qq=1.0,N=400):
    # 像電荷: +q at 2nd+z0, -q at 2nd-z0 (n整数)
    def sigma(r):
        Ey=0.0
        for n in range(-N,N+1):
            for s,zz in ((1,2*n*dd+z0),(-1,2*n*dd-z0)):
                R2=r*r+zz*zz
                Ey+= s*qq*(0-zz)/R2**1.5   # y成分 (k=1)
        return Ey/(4*np.pi)*2*np.pi*r   # sigma=eps0*E_y, sigma = eps0*E_y(上向き正)…k=1/(4pi eps0) 単位で eps0=1/(4pi)
    val,err=quad(sigma,0,60,limit=400)
    return val
for z0 in (0.2,0.5,0.7):
    num=induced_A(z0); print('数値 QA/q (x=%.1f):'%z0, round(num,4),' 理論:',-(1-z0))
    assert abs(num+(1-z0))<1.5e-2
print('OK')

# PHY-THD-053: cycle = polytropic expansion pV^n=const (V0->8V0), isobaric compression, isochoric heating; monatomic
import sympy as sp, numpy as np
n,a,R,T0=sp.symbols('n a R T0',positive=True)
g=sp.Rational(5,3); Cv=sp.Rational(3,2)*R; Cp=sp.Rational(5,2)*R
# (1) molar heat along pV^n: dQ = Cv dT + p dV ; pV=RT, pV^n=K -> p dV = -R dT/(n-1)
C=sp.simplify(Cv - R/(n-1))
print('(1) C_n =',C, '  zero at n=',sp.solve(sp.Eq(C,0),n))
# direct numeric check of C_n by integrating along process
def Cnum(nv):
    K=1.0; V=np.linspace(1,1.001,2001); p=K/V**nv; T=p*V
    W=np.trapezoid(p,V); dU=1.5*(T[-1]-T[0]); return (W+dU)/(T[-1]-T[0])
for nv in (1.3,2.0,2.5):
    print('  n=',nv,' C num',Cnum(nv),' formula',float(C.subs({R:1,n:nv})))
def cycle(nv,av):
    # R=1,T0=1 -> p0V0=1
    if nv==1: W1=np.log(av)
    else: W1=(1-av**(1-nv))/(nv-1)
    T3=av**(1-nv); Q1=W1+1.5*(T3-1)
    p3=av**(-nv); T4=p3
    W2=-p3*(av-1); Q2=2.5*(T4-T3)
    Q3=1.5*(1-T4)
    assert abs(Q1+Q2+Q3-(W1+W2))<1e-12
    Qin=sum(q for q in (Q1,Q2,Q3) if q>0)
    return (W1+W2)/Qin
print('(2) e(n=1)=',cycle(1,8),' formula',(3*np.log(2)-7/8)/(3*np.log(2)+21/16))
print('    e(n=5/3)=',cycle(5/3,8),' 58/93=',58/93)
print('    e(n=2)=',cycle(2,8),' 14/27=',14/27)
# exact
print('    exact n=5/3:',sp.nsimplify(1-g*7/(sp.Integer(8)**g-1)),'  n=2:',sp.Rational(49,64)/(sp.Rational(3,2)*(1-sp.Rational(1,64))))
# (3) limits a->oo
for nv in (1.0,1.3,1.6,5/3,2.0,3.0):
    print('(3) n=',nv,' e(a=1e8)=',cycle(nv,1e8),'  pred', (1 if nv<=5/3 else 2/(3*(nv-1))))
# symbolic limit for n>gamma
e_big=((1-a**(1-n))/(n-1)-a**(-n)*(a-1))/(sp.Rational(3,2)*(1-a**(-n)))
print('   limit n=2:',sp.limit(e_big.subs(n,2),a,sp.oo),' n=3:',sp.limit(e_big.subs(n,3),a,sp.oo))
# check that the max over n (a=8) is at n=5/3
ns=np.linspace(1,3,20001); es=[cycle(x,8) for x in ns]
print('argmax n =',ns[int(np.argmax(es))])
# sign of Q1 along process uniform: for n in (1,5/3) heat absorbed, n>5/3 released
print('Q1 sign n=1.5:',np.sign((1-8**(-0.5))/0.5+1.5*(8**-0.5-1)),' n=2:',np.sign((1-8**-1)/1+1.5*(8**-1-1)))

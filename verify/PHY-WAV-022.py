# PHY-WAV-022 2本の基線による位相差の方向探知
import numpy as np, itertools
d1,d2=4,5   # in wavelengths
frac=lambda x: x-np.floor(x)
s=0.30
pB,pC=frac(d1*s),frac(d2*s); print('readings',pB,pC)
print('(1) candidates',[(pB+k)/4 for k in range(4)],[np.degrees(np.arcsin((pB+k)/4)) for k in range(4)])
print('(2) coarse',frac(pC-pB))
def decode(rb,rc):
    sc=frac(rc-rb)
    c=np.array([(rb+k)/d1 for k in range(-1,d1+1)])
    return c[np.argmin(abs(c-sc))], sc
# (3) worst-case search: for delta just below / above 1/18
def worst_ok(delta,N=2001):
    ok=True
    for st in np.linspace(0.15,0.85,N):
        for eb,ec in itertools.product([-delta,delta],[-delta,delta]):
            rb=frac(d1*st+eb); rc=frac(d2*st+ec)
            sf,_=decode(rb,rc)
            if abs(sf-st)>delta/d1+1e-9: ok=False
    return ok
print('delta=0.055 ok?',worst_ok(0.055),' delta=0.057 ok?',worst_ok(0.057), '1/18=',1/18)
# (5) 新設問：一般の N（d1=N λ, d2=(N+1) λ）での条件 δ < 1/(2(2N+1))，精度 δ/N
def decodeN(rb,rc,N):
    sc=frac(rc-rb)
    c=np.array([(rb+k)/N for k in range(-2,N+2)])
    return c[np.argmin(abs(c-sc))]
def worst_okN(delta,N,M=1401):
    worst=0.0
    for st in np.linspace(0.15,0.85,M):
        for eb,ec in itertools.product([-delta,0,delta],[-delta,0,delta]):
            rb=frac(N*st+eb); rc=frac((N+1)*st+ec)
            sf=decodeN(rb,rc,N)
            err=abs(sf-st)
            if err>delta/N+1e-9: return False,err
            worst=max(worst,err)
    return True,worst
for N in [3,4,6,10]:
    dc=1/(2*(2*N+1))
    okb,wb=worst_okN(dc*0.99,N); oka,_=worst_okN(dc*1.02,N)
    print(f'N={N}: delta_c=1/(2(2N+1))={dc:.5f}  ok just below? {okb} (max err {wb:.5f}, delta/N={dc*0.99/N:.5f})  ok just above? {oka}')
    assert okb and not oka
    assert 2*dc<0.15   # N>=3 では小数部分の折り返しが起きない
print('N=4 ->',1/(2*9),'= 1/18')
print('max N for delta=0.02:',(1-2*0.02)/(4*0.02),' best error ~ 4 delta^2 =',4*0.02**2)

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
# (4)
lo,hi=(1.2-0.02)/4,(1.2+0.02)/4
print('(4) sin range',lo,hi,' deg',np.degrees(np.arcsin(lo)),np.degrees(np.arcsin(hi)))
print('coarse only deg',np.degrees(np.arcsin(0.26)),np.degrees(np.arcsin(0.34)))

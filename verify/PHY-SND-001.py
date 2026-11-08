# M66-1: decode closed/open tubes from resonance list; coincidence parity conditions
from fractions import Fraction as F
from itertools import product
V=340
obs=[(125,1),(375,2),(625,1),(750,1)]
def series(kind,f1,fmax=2000):
    if kind=='closed': return [f1*(2*m-1) for m in range(1,200) if f1*(2*m-1)<=fmax]
    return [f1*n for n in range(1,200) if f1*n<=fmax]
sols=[]
# brute force: candidate fundamentals on a 0.5 Hz grid up to 800
cands=[F(k,2) for k in range(1,1601)]
for kx,ky in [('closed','open'),('open','closed')]:
    for fx in cands:
        sx=series(kx,fx,760)
        if not sx or sx[0]!=125: continue
        for fy in cands:
            sy=series(ky,fy,760)
            allf=sorted(set(sx)|set(sy))[:4]
            rec=[(int(f) if f==int(f) else f,(f in sx)+(f in sy)) for f in allf]
            if rec==obs: sols.append((kx,fx,ky,fy))
print("solutions (X kind, fX, Y kind, fY):",sols)
for s in sols:
    a=F(V,4)/s[1] if s[0]=='closed' else F(V,2)/s[1]
    b=F(V,4)/s[3] if s[2]=='closed' else F(V,2)/s[3]
    print("a=",a,float(a),"b=",b,float(b))
sx=series('closed',125,5000); sy=series('open',375,5000)
print("common",sorted(set(sx)&set(sy)))
syc=series('closed',F(375,2),5000); print("common after closing Y",sorted(set(sx)&set(syc)))
# parity check, generic
import math
def closed_closed(p,q):  # a/b=p/q, fundamentals V/4a, V/4b: (2m-1)/p == (2n-1)/q
    return any((2*m-1)*q==(2*n-1)*p for m in range(1,60) for n in range(1,60))
def closed_open(p,q):    # closed a, open b, b/a=p/q: (2m-1)/(4a)=n/(2b) -> (2m-1) p = 2 n q
    return any((2*m-1)*p==2*n*q for m in range(1,60) for n in range(1,60))
ok=True
for p in range(1,12):
    for q in range(1,12):
        if math.gcd(p,q)!=1: continue
        ok&= closed_closed(p,q)==(p%2==1 and q%2==1)
        ok&= closed_open(p,q)==(p%2==0)
print("parity rules hold:",ok)

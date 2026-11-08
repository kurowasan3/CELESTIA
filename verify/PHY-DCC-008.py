# N33-1: ideal diodes. E1=12V: T node (+12). X: T-10ohm-X, X-R-G. E2=6V: U node(+6): U-10ohm-Y, Y-10ohm-G.
# D1: X->O, D2: Y->O, load O-10ohm-G.  brute force over diode states, check consistency.
import itertools, sympy as sp
def solve(Rv):
    Rv=sp.nsimplify(Rv); res=[]
    for s1,s2 in itertools.product([0,1],repeat=2):
        VX,VY,VO,i1,i2=sp.symbols('VX VY VO i1 i2')
        eqs=[sp.Eq((12-VX)/10, VX/Rv + i1), sp.Eq((6-VY)/10, VY/10 + i2), sp.Eq(i1+i2, VO/10)]
        eqs.append(sp.Eq(VX,VO) if s1 else sp.Eq(i1,0))
        eqs.append(sp.Eq(VY,VO) if s2 else sp.Eq(i2,0))
        so=sp.solve(eqs,[VX,VY,VO,i1,i2],dict=True)[0]
        ok=True
        if s1: ok&= so[i1]>=0
        else: ok&= so[VX]<=so[VO]
        if s2: ok&= so[i2]>=0
        else: ok&= so[VY]<=so[VO]
        if ok: res.append(((s1,s2),so[VO],so[i1],so[i2]))
    return res
for Rv in [1,2,3,4,5,8,30]:
    print(Rv,solve(Rv))
# symbolic boundaries
R=sp.symbols('R',positive=True)
VthX=12*R/(R+10); print('D1 on threshold (VthX=2):',sp.solve(sp.Eq(VthX,2),R))
print('D2 off threshold (X alone gives VO=3):',sp.solve(sp.Eq(VthX*10/(10*R/(R+10)+10),3),R))
# both-on: VO(R) and i2(R) formula
VO=sp.symbols('VO'); vo=sp.solve(sp.Eq((12-VO)/10-VO/R+(6-VO)/10-VO/10, VO/10),VO)[0]
print('both on VO(R)=',sp.factor(vo),' i2=',sp.factor((6-vo)/10-vo/10),' i1=',sp.factor((12-vo)/10-vo/R))

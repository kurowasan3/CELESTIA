# Q089-1: 向きの異なる2つの小さな電気双極子：トルクは逆向き・等大でないが，力の偶力を含めた合計は0
import sympy as sp
k,q,dl,L=sp.symbols('k q delta L',positive=True)
p=q*dl
# C: 原点中心，x軸向き（+q が x=+δ/2）．D: (L,0) 中心，y軸向き（+q が y=+δ/2）
C=[(q,sp.Matrix([dl/2,0])),(-q,sp.Matrix([-dl/2,0]))]
D=[(q,sp.Matrix([L,dl/2])),(-q,sp.Matrix([L,-dl/2]))]
def force(Q1,r1,Q2,r2):  # 2 が 1 に及ぼす力
    v=r1-r2; return k*Q1*Q2*v/(v.dot(v))**sp.Rational(3,2)
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
FC=sp.zeros(2,1); tC=0; FD=sp.zeros(2,1); tD=0; tot=0
cC=sp.Matrix([0,0]); cD=sp.Matrix([L,0]); mid=sp.Matrix([L/2,0])
for Q1,r1 in C:
    for Q2,r2 in D:
        f=force(Q1,r1,Q2,r2); FC+=f; tC+=cross(r1-cC,f); tot+=cross(r1-mid,f)
for Q1,r1 in D:
    for Q2,r2 in C:
        f=force(Q1,r1,Q2,r2); FD+=f; tD+=cross(r1-cD,f); tot+=cross(r1-mid,f)
ser=lambda e: sp.simplify(sp.series(e,dl,0,3).removeO())
print('F_C =',[ser(c) for c in FC]); print('F_D =',[ser(c) for c in FD])
print('tau_C (about own center) =',ser(tC)); print('tau_D =',ser(tD))
print('total torque about midpoint (exact) =',sp.simplify(tot))
assert sp.simplify(ser(FC[1])+3*k*q**2*dl**2/L**4)==0
assert sp.simplify(ser(tC)+k*q**2*dl**2/L**3)==0 and sp.simplify(ser(tD)+2*k*q**2*dl**2/L**3)==0
assert sp.simplify(tot)==0
# couple of forces about midpoint
couple=cross(cC-mid,FC)+cross(cD-mid,FD); print('couple =',ser(couple))
# 双極子の電場の公式の確認：軸上 2kp/r^3，垂直二等分線上 kp/r^3（逆向き）
r=sp.symbols('r',positive=True)
Eax=k*q/(r-dl/2)**2-k*q/(r+dl/2)**2; print('E_axis ~',sp.simplify(sp.series(Eax,dl,0,2).removeO()))
Eperp=-2*k*q/(r**2+dl**2/4)*(dl/2)/sp.sqrt(r**2+dl**2/4); print('E_perp ~',sp.simplify(sp.series(Eperp,dl,0,2).removeO()))
print('OK')

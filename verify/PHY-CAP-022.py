# O20162-1  Y-Delta equivalence of capacitors, bridge, precharged isolated node
import sympy as sp
C1,C2,C3,pA,pB,pD,C,V,q0=sp.symbols('C1 C2 C3 phiA phiB phiD C V q0',positive=True)
pO=sp.symbols('phiO')
# (1) Y with neutral center
sol=sp.solve(sp.Eq(C1*(pO-pA)+C2*(pO-pB)+C3*(pO-pD),0),pO)[0]
QA=sp.simplify(C1*(pA-sol))
S=C1+C2+C3
QA_delta=C1*C2/S*(pA-pB)+C1*C3/S*(pA-pD)
assert sp.simplify(QA-QA_delta)==0
print("phiO=",sol,"  Delta: C_AB=C1C2/S etc. OK")
# with center charge q0
sol2=sp.solve(sp.Eq(C1*(pO-pA)+C2*(pO-pB)+C3*(pO-pD),q0),pO)[0]
QA2=sp.simplify(C1*(pA-sol2))
print("QA with q0 - delta form =",sp.simplify(QA2-QA_delta))
# (2)(3) bridge: P-X C, X-Q 2C, P-Y 2C, Y-Q C, X-Y C ; phiP=V, phiQ=0
x,y=sp.symbols('x y')
def bridge(q):
    eqs=[sp.Eq(C*(x-V)+C*(x-y)+2*C*x,q), sp.Eq(2*C*(y-V)+C*(y-x)+C*y,0)]
    s=sp.solve(eqs,[x,y])
    QP=sp.simplify(C*(V-s[x])+2*C*(V-s[y]))
    return s,QP
s,QP=bridge(0)
print("phiX,phiY=",s[x],s[y]," QP=",QP," Ceq=",sp.simplify(QP/V))
# Y-Delta route: star at X: P(C),Q(2C),Y(C)
S_=4*C; CPQ=C*2*C/S_; CPY=C*C/S_; CQY=2*C*C/S_
Ceq2=sp.simplify(CPQ+ 1/(1/(2*C+CPY)+1/(C+CQY)))
print("Ceq via Y-Delta =",Ceq2)
assert sp.simplify(Ceq2-sp.Rational(7,5)*C)==0
print("bridge cap voltage",sp.simplify(s[y]-s[x]),"charge",sp.simplify(C*(s[y]-s[x])))
# energy check : sum 1/2 C V^2 == 1/2 QP V
U=sp.Rational(1,2)*(C*(V-s[x])**2+2*C*s[x]**2+2*C*(V-s[y])**2+C*s[y]**2+C*(s[y]-s[x])**2)
print("U=",sp.simplify(U)," 1/2 QP V=",sp.simplify(QP*V/2))
# (4) precharged X node q0 such that phiX=phiY
s4,QP4=bridge(q0)
q0sol=sp.solve(sp.Eq(s4[x],s4[y]),q0)[0]
print("q0 for zero bridge voltage:",q0sol,"  phi=",sp.simplify(s4[x].subs(q0,q0sol)),"  QP=",sp.simplify(QP4.subs(q0,q0sol)))
# charge flowing in at Q terminal (out of Q plates): charge on Q-side plates = -(2C*phiX + C*phiY)
QQ=sp.simplify(-(2*C*(s4[x]-0)+C*(s4[y]-0))).subs(q0,q0sol)
print("charge on Q-side plates:",sp.simplify(QQ))
# general: QP as function of q0
print("QP(q0)=",sp.expand(QP4))
# energy in (4) and work of source
U4=sp.Rational(1,2)*(C*(V-s4[x])**2+2*C*s4[x]**2+2*C*(V-s4[y])**2+C*s4[y]**2+C*(s4[y]-s4[x])**2)
U4i=sp.simplify(U4.subs(q0,q0sol))
print("U4 final =",U4i)

# N30-1: potentiometer with unknown contact resistance at end A: V(l)=k(l+u)
import sympy as sp
k,u,Ex,r=sp.symbols('k u E_x r',positive=True)
Es=sp.Rational(102,100); ls=sp.Rational(30,100); l1=sp.Rational(46,100); l2=sp.Rational(36,100); l3=sp.Rational(26,100)
R=4; Rp=sp.Rational(3,2)
sol=sp.solve([k*(ls+u)-Es, k*(l1+u)-Ex, k*(l2+u)-Ex*R/(R+r), k*(l3+u)-Ex*Rp/(Rp+r)],[k,u,Ex,r],dict=True)
print(sol)
# naive (u=0)
print('naive Ex',sp.N(Es*l1/ls),'naive r from R',sp.N(R*(l1-l2)/l2),'from Rp',sp.N(Rp*(l1-l3)/l3))
s=sol[0]
# forward check with true values
kk,uu,EE,rr=3,sp.Rational(4,100),sp.Rational(3,2),1
print([sp.Rational(1)*Es/kk-uu, EE/kk-uu, EE*R/(R+rr)/kk-uu, EE*Rp/(Rp+rr)/kk-uu])

# ---- (4) rewritten: systematic errors of the naive analysis ----
ls_,l1_,l_,u_,Es_,R_,k_=sp.symbols('l_s l_1 l u E_s R k',positive=True)
Ex_true=Es_*(l1_+u_)/(ls_+u_); Ex_naive=Es_*l1_/ls_
ratio=sp.simplify(Ex_naive/Ex_true); print('(4a) Ex\'/Ex =',ratio)
assert sp.simplify(ratio-l1_*(ls_+u_)/(ls_*(l1_+u_)))==0
print('     ratio-1 =',sp.factor(ratio-1))
assert sp.simplify(sp.factor(ratio-1)-u_*(l1_-ls_)/(ls_*(l1_+u_)))==0
rn=R_*(l1_-l_)/l_; rt=R_*(l1_-l_)/(l_+u_)
print("(4b) r'/r =",sp.simplify(rn/rt))
assert sp.simplify(rn/rt-(1+u_/l_))==0
# data: true r=1, u=4cm
for ll,RR in [(sp.Rational(36,100),4),(sp.Rational(26,100),sp.Rational(3,2))]:
    print('   l=',ll,' naive r',sp.N(RR*(l1-ll)/ll,4),'  (1+u/l)*r_true',sp.N((1+s[u]/ll)*s[r],4))
    assert abs(sp.N(RR*(l1-ll)/ll-(1+s[u]/ll)*s[r]))<1e-12
print('   naive Ex / true Ex',sp.N(Es*l1/ls/s[Ex],4),' formula',sp.N(l1*(ls+s[u])/(ls*(l1+s[u])),4))
print('DCC-005 (4) OK')

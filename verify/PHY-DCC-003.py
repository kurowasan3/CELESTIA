# N28-1: determine R, rA, rV from two connections with ideal battery
import sympy as sp
E,R,rA,rV=sp.symbols('E R r_A r_V',positive=True)
# Fig1: voltmeter across R, ammeter in series with (R||V)
I1=E/(rA+R*rV/(R+rV)); V1=E-rA*I1
# Fig2: ammeter in series with R, voltmeter across (R+A) = battery terminals
I2=E/(R+rA); V2=E
val={E:6,R:60,rA:sp.Rational(5,2),rV:228}
print('readings',[sp.nsimplify(z.subs(val)) for z in (I1,V1,I2,V2)])
v1,i1,v2,i2=sp.symbols('V1 I1 V2 I2',positive=True)
rA_e=(v2-v1)/i1; R_e=v2/i2-rA_e; rV_e=v1/(i1-v1/R_e)
for name,ex,tr in [('rA',rA_e,rA),('R',R_e,R),('rV',rV_e,rV)]:
    print(name,sp.simplify(ex.subs({v1:V1,i1:I1,v2:V2,i2:I2})-tr))
num={v1:sp.Rational(57,10),i1:sp.Rational(12,100),v2:6,i2:sp.Rational(96,1000)}
print('numeric rA,R,rV',[sp.nsimplify(z.subs(num)) for z in (rA_e,R_e,rV_e)],'naive',num[v1]/num[i1],num[v2]/num[i2])
# criterion: |err1| = R/(R+rV) vs err2 = rA/R
Rs=sp.symbols('R',positive=True)
print(sp.solve(sp.Eq(Rs/(Rs+rV),rA/Rs),Rs))
print('errors',(R/(R+rV)).subs(val),(rA/R).subs(val), 'Rc',sp.N(((rA+sp.sqrt(rA**2+4*rA*rV))/2).subs(val)))

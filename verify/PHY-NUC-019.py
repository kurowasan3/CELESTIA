# O20183-1 Bohr model with additional attractive force beta/r^3
import sympy as sp
m,k,e,h,beta,r,v,n,eps=sp.symbols('m k e h beta r v n epsilon',positive=True)
# circular motion + quantization
sol=sp.solve([sp.Eq(m*v**2/r,k*e**2/r**2+beta/r**3),sp.Eq(m*v*r,n*h/(2*sp.pi))],[r,v],dict=True)
rn=[s[r] for s in sol if s[r].is_positive!=False][0]; print("r_n=",sp.simplify(rn))
a0=h**2/(4*sp.pi**2*m*k*e**2); E0=2*sp.pi**2*m*k**2*e**4/h**2
epsdef=4*sp.pi**2*m*beta/h**2
print("check r_n=a0(n^2-eps):",sp.simplify(rn-a0*(n**2-epsdef)))
# potential of beta/r^3 attraction: U=-beta/(2r^2); energy
U=-k*e**2/r-beta/(2*r**2)
vn=n*h/(2*sp.pi*m*rn)
En=sp.simplify(sp.Rational(1,2)*m*vn**2+U.subs(r,rn)); print("E_n=",En)
print("check E_n=-E0/(n^2-eps):",sp.simplify(En+E0/(n**2-epsdef)))
print("check E=-ke^2/(2r):",sp.simplify(En+k*e**2/(2*rn)))
# force from U
print("F=-dU/dr :",sp.simplify(-sp.diff(-beta/(2*r**2),r)))
# line ratio lambda(3->2)/lambda(2->1)
x=sp.symbols('x')
E=lambda N:-1/(N**2-x)
rho=(E(2)-E(1))/(E(3)-E(2))
ser=sp.series(rho,x,0,2).removeO(); print("rho series:",sp.simplify(ser), " c=",sp.simplify(sp.diff(rho,x).subs(x,0)/rho.subs(x,0)))
# relative shifts of each line
for (a,b) in [(2,1),(3,2),(3,1),(4,3)]:
    dE=E(a)-E(b); print(f"{a}->{b} relative change of dE: ",sp.simplify(sp.diff(dE,x).subs(x,0)/dE.subs(x,0)))
# order of magnitude used in the new (5) remarks: beta=G'Mm, G'=G R ; eps = 4pi^2 m G R M m / h^2
G=6.67e-11; M=1.67e-27; me=9.11e-31; hh=6.63e-34; R=2e-4
epsn=4*3.1416**2*G*R*M*me**2/hh**2; print("eps numeric (R=2e-4):",epsn)
# critical beta
print("beta_c=",sp.solve(sp.Eq(epsdef,1),beta))
# numeric check of ratio with eps=1e-3
ep=1e-3
Ev=lambda N:-1/(N**2-ep)
rho_ex=(Ev(2)-Ev(1))/(Ev(3)-Ev(2)); print("rho exact",rho_ex," approx",27/5*(1+8/9*ep))
print("eps with given rounded constants:",4*3.14159**2*6.7e-11*2.0e-4*1.7e-27*(9.1e-31)**2/(6.6e-34)**2)

# ---- new (5): eps as force ratio, factorization, detection condition ----
print("eps = beta/(required centripetal force * r^3) at n=1:",sp.simplify(epsdef-beta/(h**2/(4*sp.pi**2*m))))
Gs,Ms,Rs,eta=sp.symbols('G M R eta',positive=True)
eps_h=epsdef.subs(beta,Gs*Rs*Ms*m)
print("eps - (GMm/ke^2)(R/a0) =",sp.simplify(eps_h-(Gs*Ms*m/(k*e**2))*(Rs/a0)))
Rcond=sp.solve(sp.Eq(sp.Rational(8,9)*eps_h,eta),Rs)[0]
print("R threshold - 9eta/8 * ke^2/(GMm) * a0 =",sp.simplify(Rcond-sp.Rational(9,8)*eta*k*e**2/(Gs*Ms*m)*a0))
ke2=8.99e9*(1.602e-19)**2; ratio=G*M*me/ke2; a0n=5.29e-11
print("GMm/ke^2 =",ratio,"  R needed for eta=1e-10 [in a0]:",9/8*1e-10/ratio)
print("check eps(R=2e-4) via factorization:",ratio*R/a0n)

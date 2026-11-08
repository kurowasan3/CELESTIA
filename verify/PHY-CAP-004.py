# N19-1: dielectric liquid rising between vertical plates
import sympy as sp
e,er,w,d,L,V,rho,g,h,Q=sp.symbols('epsilon0 epsilon_r w d L V rho g h Q',positive=True)
C=e*w*(L+(er-1)*h)/d
# constant V: upward force = d/dh(1/2 C V^2)
F=sp.diff(C*V**2/2,h); print('F_V',F)
hV=sp.solve(sp.Eq(F,rho*g*w*d*h),h)[0]; print('hV',hV)
# energy bookkeeping at hV
Wb=V*(C.subs(h,hV)-C.subs(h,0))*V
dU=sp.Rational(1,2)*(C.subs(h,hV)-C.subs(h,0))*V**2
Ug=rho*g*w*d*hV**2/2
print('dU/Wb',sp.simplify(dU/Wb),'Ug/Wb',sp.simplify(Ug/Wb),'heat/Wb',sp.simplify((Wb-dU-Ug)/Wb))
# constant Q: Q = C(0) V0 ; force = d/dh(-Q^2/(2C)) * (-1)... upward force = Q^2/(2C^2) dC/dh
V0=sp.symbols('V0',positive=True)
Qv=C.subs(h,0)*V0
FQ=-sp.diff(Qv**2/(2*C),h)
print('FQ',sp.simplify(FQ))
eqQ=sp.simplify(FQ-rho*g*w*d*h)
H0=sp.symbols('H0',positive=True)  # H0 = hV evaluated at V0
H0expr=hV.subs(V,V0)
# claim: h (1+(er-1)h/L)^2 = H0
print('check claim', sp.simplify(eqQ.subs(V0, sp.sqrt(H0*2*rho*g*d**2/((er-1)*e))) * (L+(er-1)*h)**2 - (rho*g*w*d)*(H0*L**2 - h*(L+(er-1)*h)**2)) )
# numeric example: er=5? choose (er-1)h1/L = 1/2  -> H0 = h1*(9/4)
er_v=3; L_v=1.0; h1=L_v/(2*(er_v-1))
print('H0/h1', (1+(er_v-1)*h1/L_v)**2)
# stability of constant-Q equilibrium: d/dh(FQ - weight) <0
print('dF/dh', sp.simplify(sp.diff(FQ-rho*g*w*d*h,h)))

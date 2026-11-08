# J095-1: slit camera with point source at finite distance a: geometric + diffraction blur, resolution limit
import sympy as sp, numpy as np
a,f,lam,d,s=sp.symbols('a f lambda d s',positive=True)
geo=d*(a+f)/a   # rays from point source through slit edges, similar triangles
dif=2*f*lam/d   # first zeros at +-lam/d
D=geo+dif
dstar=sp.solve(sp.diff(D,d),d)[0]; Dmin=sp.simplify(D.subs(d,dstar))
print('d* =',sp.simplify(dstar),' Dmin =',Dmin)
F=a*f/(a+f); print('d*^2 - 2 lam F =',sp.simplify(dstar**2-2*lam*F))
# check geometric width by explicit rays: source at (−a, 0), slit edges at (0,±d/2), screen x=f
yedge=lambda ye: ye+ (ye-0)/a*f
print('geo width check', sp.simplify(yedge(d/2)-yedge(-d/2)-geo))
# image separation of two sources separated s: s f / a ; resolution criterion s f/a >= Dmin
smin=sp.simplify(Dmin*a/f); print('s_min =',smin)
print('limit f->oo', sp.limit(smin,f,sp.oo))
print('d/df s_min sign:', sp.simplify(sp.diff(smin,f)))
# first diffraction zero: pair rays d/2 apart, path diff d/2 sin th = lam/2 -> sin th = lam/d
vals={lam:5.0e-7,a:10.0,f:0.20}
print('d* (mm)=',float(dstar.subs(vals))*1e3,' Dmin(mm)=',float(Dmin.subs(vals))*1e3,' smin(mm)=',float(smin.subs(vals))*1e3)
print('f->oo smin (mm)=',float((2*sp.sqrt(2*lam*a)).subs(vals))*1e3)
# numeric minimization check
dd=np.linspace(1e-5,2e-3,200001); DD=dd*(10.2/10)+2*0.2*5e-7/dd; k=DD.argmin(); print('numeric d*',dd[k]*1e3,'Dmin',DD[k]*1e3)
# if f doubled -> 0.40
vals2={lam:5.0e-7,a:10.0,f:0.40}; print('f=0.4: smin(mm)=',float(smin.subs(vals2))*1e3)

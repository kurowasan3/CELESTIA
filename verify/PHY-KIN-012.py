# O20181-1 car clearing a wall at distance d, height h: ramp redirect vs vertical impulse
import sympy as sp, numpy as np
g,h,d,a,v,u,t,T=sp.symbols('g h d a v u t T',positive=True)
# (1) 45 deg ramp: y(d)=d - g d^2/v^2 >= h
v45=sp.solve(sp.Eq(d-g*d**2/v**2,h),v)[0]; s1=sp.simplify(v45**2/(2*a)); print("s1=",s1)
# (2) free angle: discriminant
T_=sp.symbols('tau')
expr=d*T_ - g*d**2*(1+T_**2)/(2*v**2) - h
disc=sp.discriminant(sp.expand(expr*2*v**2),T_)
vmin2=sp.solve(sp.Eq(sp.simplify(disc),0),v**2); print("vmin^2 candidates:",vmin2)
vm2=g*(h+sp.sqrt(h**2+d**2)); tan_opt=sp.simplify(v**2/(g*d)).subs(v**2,vm2); print("tan opt=",tan_opt)
# numeric brute force
gv,hv,dv=9.8,2.0,6.0
best=min((np.sqrt(gv*dv**2/(2*np.cos(th)**2*(dv*np.tan(th)-hv))),th) for th in np.linspace(np.arctan(hv/dv)+1e-4,np.pi/2-1e-4,200000))
print("numeric vmin^2",best[0]**2," formula",gv*(hv+np.sqrt(hv**2+dv**2))," tan",np.tan(best[1]),(hv+np.sqrt(hv**2+dv**2))/dv)
print("s1 > s2 ?", gv*dv**2/(dv-hv)/2, gv*(hv+np.sqrt(hv**2+dv**2))/2)
# (3) vertical impulse giving vertical speed u, horizontal vx
ts=sp.solve(sp.Eq(u*t-g*t**2/2,h),t); print("times:",ts)
# vx range d/t2 .. d/t1 ; s = vx^2/(2a)
t1=(u-sp.sqrt(u**2-2*g*h))/g; t2=(u+sp.sqrt(u**2-2*g*h))/g
smin=sp.simplify((d/t2)**2/(2*a)); smax=sp.simplify((d/t1)**2/(2*a)); print("smin=",smin,"smax=",smax)
# (4) u = 2 sqrt(gh)
U=2*sp.sqrt(g*h)
print("u=2sqrt(gh): smin=",sp.simplify(smin.subs(u,U)),"  smax=",sp.simplify(smax.subs(u,U)), " ratio",sp.simplify(smax.subs(u,U)/smin.subs(u,U)))
# simulate check: for vx in range clears, outside fails
uv=2*np.sqrt(gv*hv)
t1v=(uv-np.sqrt(uv**2-2*gv*hv))/gv; t2v=(uv+np.sqrt(uv**2-2*gv*hv))/gv
for vx in [dv/t2v*0.98,dv/t2v*1.02,dv/t1v*0.98,dv/t1v*1.02]:
    tt=dv/vx; y=uv*tt-gv*tt**2/2; print("vx",vx,"height at wall",y,"clears" if y>=hv else "fails")


# ---- (5) total energy (per unit mass) E = a s + u^2/2 minimised over u and s ----
D=sp.sqrt(h**2+d**2)
s2=g*(h+D)/(2*a)
Uopt=(h+D)*sp.sqrt(g/(2*D)); sopt=g*d**2/(4*a*D)
# trajectory with vx^2=2 a sopt, vertical u=Uopt passes exactly through (d,h)
vx=sp.sqrt(2*a*sopt); tw=d/vx
print('(5) height at wall - h:',sp.simplify(Uopt*tw-g*tw**2/2-h))
assert sp.simplify(Uopt*tw-g*tw**2/2-h)==0
print('(5) E(opt) - a s2:',sp.simplify(a*sopt+Uopt**2/2-a*s2)); assert sp.simplify(a*sopt+Uopt**2/2-a*s2)==0
assert sp.simplify((Uopt**2-2*g*h)-g*(D-h)**2/(2*D))==0
# numeric brute force over u: min of a*smin(u)+u^2/2 (smin = lowest run-up for that u)
smin_f=sp.lambdify((u,),(smin+u**2/(2*a)).subs({g:gv,h:hv,d:dv,a:1.0}))
us=np.linspace(np.sqrt(2*gv*hv)+1e-9,40,400001); Es=smin_f(us)
i=np.argmin(Es)
print('(5) numeric min E=%.6f at u=%.5f ; formula E=%.6f u=%.5f'%(Es[i],us[i],gv*(hv+np.hypot(hv,dv))/2,float(Uopt.subs({g:gv,h:hv,d:dv}))))
assert abs(Es[i]-gv*(hv+np.hypot(hv,dv))/2)<1e-6
# also: E over the whole allowed s-range for each u is >= a s2 (any launch reaching wall needs V^2>=g(h+D))
smax_f=sp.lambdify((u,),(smax+u**2/(2*a)).subs({g:gv,h:hv,d:dv,a:1.0}))
assert np.all(smax_f(us[::100])>=gv*(hv+np.hypot(hv,dv))/2-1e-9)

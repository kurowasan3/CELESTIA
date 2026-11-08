# M49-1: slowly expanding mirrored sphere (v r = const, theta const) and collisionless piston expansion
import sympy as sp, numpy as np
v,u,th,r,T=sp.symbols('v u theta r T',positive=True)
# (1) collision with wall receding at u along normal: normal comp v cos(th) -> -(v cos th - 2u)
vp2=(v*sp.sin(th))**2+(v*sp.cos(th)-2*u)**2
dv=sp.series(sp.sqrt(vp2),u,0,2).removeO()-v
print('dv first order',sp.simplify(dv))
# (2) rate
rate=sp.simplify(dv*v/(2*r*sp.cos(th))); print('dv/dt',rate)
# ODE dv/dr = rate/u -> v r = const
f=sp.Function('f'); sol=sp.dsolve(sp.Eq(f(r).diff(r),-f(r)/r)); print(sol)
# numeric simulation in sphere, event driven
def sim(r0=1.0,u0=1e-4,v0=1.0,steps=None,rfinal=2.0):
    rng=np.random.default_rng(3)
    pos=np.array([0.3,0.1,-0.2]); d=rng.normal(size=3); d/=np.linalg.norm(d); vel=v0*d
    t=0.0
    while r0+u0*t<rfinal:
        # solve |pos+vel s| = r0+u0(t+s)
        a=vel@vel-u0*u0; R0=r0+u0*t
        b=2*(pos@vel-R0*u0); c=pos@pos-R0*R0
        s=(-b+np.sqrt(b*b-4*a*c))/(2*a)
        pos=pos+vel*s; t+=s
        n=pos/np.linalg.norm(pos); vn=vel@n
        vel=vel-2*(vn-u0)*n   # reflect relative to moving wall
        pos=pos*(1-1e-12)
    R=r0+u0*t
    sinth=np.linalg.norm(np.cross(pos,vel))/(np.linalg.norm(pos)*np.linalg.norm(vel))
    return R,np.linalg.norm(vel),sinth
# initial sin(theta) at first hit? record after first hit instead
R,vf,sf=sim()
print('sphere: R',R,'v*R',vf*R)
# theta check: run short
def sim_theta(rfinal):
    return sim(rfinal=rfinal)[2]
print('sin theta at r=1.001, 1.5, 2:',sim_theta(1.001),sim_theta(1.5),sim_theta(2.0))
# cylinder: collisionless, piston at x=Lp(t), only vx changes
rng=np.random.default_rng(5); Nm=4000
vxs=rng.normal(size=Nm); vys=rng.normal(size=Nm); vzs=rng.normal(size=Nm)
E0=np.mean(vxs**2+vys**2+vzs**2); P0=np.mean(vxs**2)/1.0
u0=2e-3; vx_f=[]
for vx in vxs:
    x=rng.random(); w=vx; t=0.0
    while True:
        Lp=1.0+u0*t
        if w>0:
            s=(Lp-x)/(w-u0)
            if w<=u0: break
            t+=s; x=1.0+u0*t
            if x>=2.0: break
            w=-(w-2*u0)
        else:
            s=x/(-w); t+=s; x=0.0
            if 1.0+u0*t>=2.0: break
            w=-w
    vx_f.append(abs(w))
vx_f=np.array(vx_f)
E1=np.mean(vx_f**2+vys**2+vzs**2)
print('cyl E1/E0',E1/E0,'(expect 3/4)','P1/P0',np.mean(vx_f**2)/2.0/P0,'(expect 1/8)')
print('quasi-static adiabatic T ratio',2**(-2/3),' P ratio',2**(-5/3),' isotropized P',0.75/2)

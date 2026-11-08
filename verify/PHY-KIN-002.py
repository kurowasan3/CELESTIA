# M03-1: bouncing up a smooth incline
import sympy as sp, math
g,v0,al,th,e,t=sp.symbols('g v0 alpha theta e t',positive=True)
t1=2*v0*sp.sin(al)/(g*sp.cos(th)); T=sp.summation(t1*sp.Symbol('e')**sp.Symbol('k'),(sp.Symbol('k'),0,sp.oo))
T=t1/(1-e)
ttop=v0*sp.cos(al)/(g*sp.sin(th)); xmax=sp.simplify(v0*sp.cos(al)*ttop-g*sp.sin(th)*ttop**2/2)
print('xmax',xmax)
print('T<=ttop  <=>', sp.simplify(sp.solve(sp.Eq(T,ttop),e)[0]), '(e at boundary) -> tan a tan th = (1-e)/2')
print('T<=2ttop boundary e:', sp.simplify(sp.solve(sp.Eq(T,2*ttop),e)[0]))
# exact event-based simulation
def sim(v0v,alv,thv,ev,gv=9.8):
    gx=-gv*math.sin(thv); gy=-gv*math.cos(thv)
    tt=0; x=0; vx=v0v*math.cos(alv); vy=v0v*math.sin(alv)
    xmax=0; log=[]
    while vy>1e-12:
        dt=-2*vy/gy
        # max x during this hop
        if vx>0 and vx/(-gx)<dt: xmax=max(xmax,x+vx**2/(-2*gx))
        x+=vx*dt+0.5*gx*dt*dt; vx+=gx*dt; tt+=dt; xmax=max(xmax,x)
        vy=ev*vy; log.append((tt,x))
        if len(log)>2000: break
    if vx>0: xmax=max(xmax,x+vx**2/(-2*gx))
    return tt,xmax,vx,x,log
for (alv,thv,ev) in [(0.2,0.5,0.5),(0.6,0.5,0.5),(0.35,0.6,0.3)]:
    v0v=5; gv=9.8
    tt,xm,vx,x,log=sim(v0v,alv,thv,ev)
    Tp=2*v0v*math.sin(alv)/(gv*math.cos(thv))/(1-ev); tt_top=v0v*math.cos(alv)/(gv*math.sin(thv))
    print(f'al={alv},th={thv},e={ev}: T sim {tt:.6f} pred {Tp:.6f}; xmax sim {xm:.6f} pred {v0v**2*math.cos(alv)**2/(2*gv*math.sin(thv)):.6f}; tan a tan th={math.tan(alv)*math.tan(thv):.3f} (1-e)/2={(1-ev)/2}, 1-e={1-ev}; stops before top? {tt<tt_top}; sliding at return? {tt<2*tt_top}')
    if tt<2*tt_top:
        # sliding from tt: speed at x=0
        gx=-gv*math.sin(thv); vret=-math.sqrt(vx**2+2*(-gx)*x) if True else None
        print('   speed at O on return:',abs(vret),' pred v0 cos a=',v0v*math.cos(alv),' energy lost/m:',0.5*v0v**2-0.5*vret**2,' pred',0.5*(v0v*math.sin(alv))**2)

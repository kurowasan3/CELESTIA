# 2001[2] 類題：一方の音源が動く2音源の干渉模様（節の移動速度とうなり）
import sympy as sp, numpy as np
f,V,u,L,x,t,te,w=sp.symbols('f V u L x t t_e w',positive=True)
# S1: x_s=-L+u t。到達条件 x - x_s(te) = V (t-te)
tes=sp.solve(sp.Eq(x-(-L+u*te),V*(t-te)),te)[0]
Phi1=sp.simplify(2*sp.pi*f*tes); print('Phi1=',Phi1)
Phi2=2*sp.pi*f*(t-(L-x)/V)
D=sp.expand(Phi1-Phi2)
a=sp.simplify(sp.diff(D,t)); b=sp.simplify(sp.diff(D,x))
vn=sp.simplify(-a/b); print('node speed',vn)
assert sp.simplify(vn-u*V/(2*V-u))==0
spacing=sp.simplify(2*sp.pi/(-b)); print('spacing',spacing)
fb=sp.simplify(a/(2*sp.pi)); print('beat at fixed mic',fb)
assert sp.simplify(fb-vn/spacing)==0
# 動くマイク：ドップラーで両方の振動数が等しくなる w
f1=f*V/(V-u)*(V-w)/V; f2=f*(V+w)/V
ws=sp.solve(sp.Eq(f1,f2),w); print('w no beat',ws)
assert sp.simplify(ws[0]-vn)==0
# 数値：包絡線の零点を追跡
fn,Vn,un,Ln=200.0,340.0,40.0,50.0
xs=np.linspace(-10,10,200001)
def env(tt):
    P1=2*np.pi*fn*(Vn*tt-xs-Ln)/(Vn-un); P2=2*np.pi*fn*(tt-(Ln-xs)/Vn)
    return np.abs(np.cos((P1-P2)/2))
def node_near(tt,x0):
    e=env(tt); sel=np.abs(xs-x0)<0.6
    return xs[sel][np.argmin(e[sel])]
t0=0.3
n0=node_near(t0,0.0); n1=node_near(t0+0.002,n0+0.002*un/2)
print('numeric node speed',(n1-n0)/0.002,' theory',un*Vn/(2*Vn-un))
assert abs((n1-n0)/0.002-un*Vn/(2*Vn-un))<0.5

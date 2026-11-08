# N02-1: two plane waves in a channel, trace wavelength along wall
import sympy as sp, numpy as np
lam,th,f,D,Lam=sp.symbols('lambda theta f D Lambda',positive=True)
# superpose cos(k(x cos th + y sin th) - wt) + cos(k(x cos th - y sin th) - wt) = 2 cos(k y sin th) cos(k x cos th - wt)
x,y,t=sp.symbols('x y t',real=True); k=2*sp.pi/lam; w=2*sp.pi*f
s=sp.cos(k*(x*sp.cos(th)+y*sp.sin(th))-w*t)+sp.cos(k*(x*sp.cos(th)-y*sp.sin(th))-w*t)
print(sp.simplify(sp.expand_trig(s - 2*sp.cos(k*y*sp.sin(th))*sp.cos(k*x*sp.cos(th)-w*t))))
# wall antinodes at y=0,D with one node: k D sin th = pi
sol=sp.solve([sp.Eq(D, lam/(2*sp.sin(th))), sp.Eq(Lam, lam/sp.cos(th))],[lam,th],dict=True)
print(sol)
Dv,fv,Lv=0.30,2.0,0.80
lv=2*Dv*Lv/np.sqrt(4*Dv**2+Lv**2); v=fv*lv
thv=np.arcsin(lv/(2*Dv))
print("lambda",lv,"v",v,"Vp",v/np.cos(thv),"fc",v/(2*Dv), "check Lam",lv/np.cos(thv))
# numeric check of pattern: along wall y=0, crest positions spacing
xs=np.linspace(0,5,500001); 
def field(xx,yy,tt): 
    kk=2*np.pi/lv; return np.cos(kk*(xx*np.cos(thv)+yy*np.sin(thv))-2*np.pi*fv*tt)+np.cos(kk*(xx*np.cos(thv)-yy*np.sin(thv))-2*np.pi*fv*tt)
z=field(xs,0,0); pk=xs[1:-1][(z[1:-1]>z[:-2])&(z[1:-1]>z[2:])]; print("crest spacing",np.diff(pk)[:3])
print("amp at y=D/2:",np.max(np.abs(field(xs,Dv/2,0.1))), " at y=D:",np.max(np.abs(field(xs,Dv,0.1))))

# (4) 新設問：遮断振動数 f_c=v/(2D)，エネルギーの x 方向の速さ U=v cos(th)，U V = v^2
fs_,vs_,Ds_=sp.symbols('f v D',positive=True)
fc_=vs_/(2*Ds_)
sth=fc_/fs_                          # sin(theta)=lambda/(2D)=f_c/f
U_=vs_*sp.sqrt(1-sth**2); V_=vs_/sp.sqrt(1-sth**2)
assert sp.simplify(U_*V_-vs_**2)==0
print('U at f->fc:',sp.limit(U_,fs_,fc_,'+'),'  V at f->fc:',sp.limit(V_,fs_,fc_,'+'))
# 数値：浅水波の線形理論（水深 h=v^2/g）でエネルギー流束を直接計算し U=v cos(th) を確かめる
#   各平面波 η_i と流速 u=(g/v)(η_1 n_1+η_2 n_2)，エネルギー密度 E=½ρgη²+½ρh|u|²，流束 F=ρ g h η u
#   ρg で割ると E=½[η²+|η_1 n_1+η_2 n_2|²]，F=v η (η_1 n_1+η_2 n_2)。x=0 の断面で水路幅・1周期の平均をとる
vn,Dn=0.96,0.30
for fn in [3.0,2.0,1.7,1.61]:
    lam_=vn/fn; th_=np.arcsin(lam_/(2*Dn)); k_=2*np.pi/lam_
    ky=k_*np.sin(th_)
    yy=np.linspace(0,Dn,4001); tt=np.linspace(0,1/fn,400,endpoint=False)
    Y,Tt=np.meshgrid(yy,tt)
    e1=np.cos(ky*Y-2*np.pi*fn*Tt); e2=np.cos(-ky*Y-2*np.pi*fn*Tt)
    eta=e1+e2
    flux=np.mean(vn*eta*(e1+e2)*np.cos(th_))
    dens=np.mean(0.5*(eta**2+((e1+e2)*np.cos(th_))**2+((e1-e2)*np.sin(th_))**2))
    Unum=flux/dens; Uth=vn*np.cos(th_)
    print('f=',fn,' U numeric',Unum,' v cos(th)',Uth,' V',vn/np.cos(th_),' U*V',Unum*vn/np.cos(th_),' v^2',vn**2)
    assert abs(Unum-Uth)<1e-3*vn
print('PHY-WAV-007 (4) OK')

# T20133-1 : 2層の板を伝わる縦波・横波の先頭波（屈折波）による走時と逆問題
import numpy as np, sympy as sp
from scipy.optimize import minimize_scalar
VA,VB,k,h=3.0,5.0,sp.Rational(4,3),12.0   # mm/us, mm
VBs=VB/float(k)
print('VB/k =',VBs,'> VA')
# 先頭波: フェルマーの原理で経路 (斜めに下りる→境界を速さcで→斜めに上る) の時間最小を数値で
def t_head(x,c):
    T=lambda u: 2*np.hypot(h,u)/VA+(x-2*u)/c
    r=minimize_scalar(T,bounds=(0,x/2),method='bounded'); return r.fun,r.x
for c,name in [(VB,'P'),(VBs,'S')]:
    for x in [60,100]:
        tn,u=t_head(x,c); th=x/c+2*h*np.sqrt(1/VA**2-1/c**2)
        print(name,'x=',x,' num %.5f th %.5f  sin(alpha)=%.4f VA/c=%.4f'%(tn,th,u/np.hypot(h,u),VA/c))
# 切片と最小距離
for c,name in [(VB,'P'),(VBs,'S')]:
    sa=VA/c; tau=2*h*np.sqrt(1/VA**2-1/c**2); xmin=2*h*sa/np.sqrt(1-sa**2)
    print(name,' intercept %.4f us  xmin %.4f mm'%(tau,xmin))
# 逆問題: 傾き1/5.0, 切片6.4 から h
hh=sp.symbols('h',positive=True)
print('h from data:',sp.solve(sp.Eq(2*hh*sp.sqrt(sp.Rational(1,9)-sp.Rational(1,25)),sp.Rational(32,5)),hh))
# 初動（最初に届く波）
xs=np.linspace(0.1,300,30000)
ta=xs/VA
tb=np.where(xs>=18,xs/VB+6.4,np.inf); tc=np.where(xs>=32,xs/VBs+4.8,np.inf)
first=np.argmin(np.vstack([ta,tb,tc]),axis=0)
ch=np.where(np.diff(first)!=0)[0]
print('first-arrival switches at x=',xs[ch],' kinds',first[ch],first[ch+1],' c ever first?',(first==2).any())
# 反射縦波 2√(h²+(x/2)²)/VA は直接波より常に遅い
print('reflected never first:',np.all(2*np.hypot(h,xs/2)/VA>ta))
# 横波の屈折角の確認（境界に沿った波面の交線の速さ）
a=sp.symbols('alpha'); print('sin phi\' =',sp.simplify((VBs)*sp.sin(a)/VA),' sin theta\'=',sp.sin(a)/k)

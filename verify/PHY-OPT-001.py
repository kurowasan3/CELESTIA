# M73-1 検証：スラブ型光導波路の受光角と到達時間の広がり（入射側媒質で場合分け）
import sympy as sp, numpy as np
n0,n1,n2,l,c,al=sp.symbols('n0 n1 n2 l c alpha',positive=True)
# 入射角alpha -> 内部屈折角beta: n0 sin a = n1 sin b; 側面入射角 90-b > 臨界角
sinb=n0*sp.sin(al)/n1
cond=sp.simplify(n1**2*(1-sinb**2)-n2**2)   # cos^2 b > (n2/n1)^2
print('TIR condition: n1^2 cos^2b - n2^2 =',sp.expand(cond),'>0  <=> n0^2 sin^2 a < n1^2-n2^2')
t=lambda cb: n1*l/(c*cb)
dtA=sp.simplify(t(n2/n1)-t(1)); dtB=sp.simplify(t(sp.sqrt(n1**2-n0**2)/n1)-t(1))
print('case n0>=sqrt(n1^2-n2^2):',dtA); print('case n0<sqrt:',dtB)
# 数値：alpha を走査して導波される光線の最大到達時間を直接探す
def spread(n0v,n1v,n2v,lv=1e3,cv=3e8):
    a=np.linspace(0,np.pi/2*0.999999,2000001)
    sb=n0v*np.sin(a)/n1v; ok=sb<1
    cb=np.sqrt(1-sb[ok]**2)
    guided=cb>n2v/n1v   # 側面入射角の正弦 cos b > n2/n1
    tt=n1v*lv/(cv*cb[guided])
    return tt.max()-tt.min()
for n0v in [1.00,1.33]:
    print(n0v,'sqrt(n1^2-n2^2)=',np.sqrt(1.6**2-1.2**2),'num dt=',spread(n0v,1.6,1.2))
print('formula air:',1.6*1e3/3e8*(1.6/np.sqrt(1.6**2-1)-1),' water:',1.6*1e3/3e8*(1.6/1.2-1))
# (4) 受光角
print('alpha_max water deg:',np.degrees(np.arcsin(np.sqrt(1.12)/1.33)))
# 比
print('ratio water/air:',(1.6/1.2-1)/(1.6/np.sqrt(1.56)-1))

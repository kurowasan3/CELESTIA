# T20142-1 : 一部が影になった太陽電池の直列接続とバイパスダイオード（理想化モデル）
import numpy as np
from scipy.optimize import brentq
s,P,V0=1.0,1.0,1.0
r=1e-4; eps=1e-4          # r: V>V0 側の傾き（十分小），eps: V<V0 側のわずかな傾き（正則化）
def Vcell(I,Pk):
    Isc=s*Pk
    return V0+r*(Isc-I) if I<=Isc else V0-(I-Isc)/eps
def Vunit(I,Pk,diode):
    V=Vcell(I,Pk); return max(V,0.0) if diode else V
def solve(R,alpha,diode):
    F=lambda I: Vunit(I,P,diode)+Vunit(I,alpha*P,diode)-I*R
    I=brentq(F,0,5); return I,Vunit(I,P,diode),Vunit(I,alpha*P,diode)
def I_theory(R,alpha,diode):
    if not diode:
        return alpha*s*P if R<=2*V0/(alpha*s*P) else 2*V0/R
    if R<=V0/(s*P): return s*P
    if R<V0/(alpha*s*P): return V0/R
    if R<=2*V0/(alpha*s*P): return alpha*s*P
    return 2*V0/R
alpha=0.3; err=0
for diode in [False,True]:
    for R in np.linspace(0.05,12,240):
        I,V1,V2=solve(R,alpha,diode); err=max(err,abs(I-I_theory(R,alpha,diode)))
print('max |I-I_theory| =',err)
I,V1,V2=solve(1e-6,alpha,False); print('no diode R=0: I=%.4f V1=%.4f V2=%.4f  power absorbed by cell2=%.4f (th %.4f)'%(I,V1,V2,-I*V2,alpha*s*P*V0))
for alpha in [0.3,0.5,0.7]:
    Rs=np.linspace(0.01,20,40000)
    Pw=[solve(R,alpha,True)[0]**2*R for R in Rs]; k=int(np.argmax(Pw))
    print('alpha',alpha,' Pmax num %.4f th %.4f at R=%.3f'%(Pw[k],max(s*P*V0,2*alpha*s*P*V0),Rs[k]))
# コンデンサー充電（ダイオードあり）: dQ/dt = I(V_C)
alpha=0.3; C=1.0
def Ich(Vc):
    G=lambda I: Vunit(I,P,True)+Vunit(I,alpha*P,True)-Vc
    if G(0)<=0: return 0.0
    return brentq(G,0,5)
t,Q,dt=0.0,0.0,1e-4; marks=[]
while Q/C<2*V0*0.9999 and t<20:
    Q+=Ich(Q/C)*dt; t+=dt
    if not marks and Q/C>=V0*0.9999: marks.append(t)
print('time to V0: %.3f (th %.3f)  time to 2V0: %.3f (th %.3f)'%(marks[0],C*V0/(s*P),t,C*V0/(s*P)+C*V0/(alpha*s*P)))

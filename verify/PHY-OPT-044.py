# PHY-OPT-044 検証：重ねた2枚の回折格子と並べた2枚の回折格子（レンズの焦点面＝フーリエ面）
import numpy as np, sympy as sp
N=1024
x=np.arange(N); X,Y=np.meshgrid(x,x,indexing='xy')
kA=np.array([64,0]); kB=np.array([63,11])   # 空間周波数（ビン単位）
def grating(k,binary=True):
    ph=2*np.pi*(k[0]*X+k[1]*Y)/N
    return (np.cos(ph)>np.cos(np.pi*0.4)).astype(float) if binary else 0.5*(1+np.cos(ph))
def peaks(t,thr=1e-3):
    F=np.abs(np.fft.fft2(t))**2; F/=F.max()
    idx=np.argwhere(F>thr)
    res=[]
    for iy,ix in idx:
        kx=ix if ix<N//2 else ix-N; ky=iy if iy<N//2 else iy-N
        res.append((kx,ky,F[iy,ix]))
    return res
A=grating(kA); B=grating(kB)
stack=A*B
pk=peaks(stack)
S={(a,b) for a,b,_ in pk}
diff=tuple(kA-kB)
print('重ねた場合 a-b のピーク',diff, diff in S or (-diff[0],-diff[1]) in S)
assert diff in S
# 並べた場合（上半分A，下半分B）
side=np.where(Y<N//2,A,B)
F=np.abs(np.fft.fft2(side))**2; F/=F.max()
v=F[diff[1]%N,diff[0]%N]
vA=F[0,kA[0]]
print('並べた場合 a-b の強さ',v,' 1次(a)の強さ',vA)
assert v<1e-3*vA
# 重ねた場合の格子点 m a + n b がすべてピーク候補であること（強いもの）
for a,b,val in pk:
    if val>1e-2:
        # 解 m kA + n kB = (a,b) が整数
        M=np.array([[kA[0],kB[0]],[kA[1],kB[1]]],float)
        mn=np.linalg.solve(M,[a,b])
        assert np.allclose(mn,np.round(mn),atol=1e-6),(a,b,mn)
print('強いピークはすべて m a + n b 上')
# 等周期・角α：最近接の新しい明点の距離 2 s sin(α/2) と方向
s,al=sp.symbols('s alpha',positive=True)
a=sp.Matrix([s,0]); b=sp.Matrix([s*sp.cos(al),s*sp.sin(al)])
dvec=a-b
print('|a-b|=',sp.simplify(sp.sqrt(sp.trigsimp(dvec.dot(dvec)))))
# 60°未満なら |a-b|<s
for alv in np.linspace(0.01,np.pi/3-0.01,30):
    assert 2*np.sin(alv/2)<1
# モアレ周期 D=d/(2sin(α/2))：2次元で重ねた透過率を低域通過し，残る縞の周期と向きを測る
N2=2048
x2=np.arange(N2); X2,Y2=np.meshgrid(x2,x2,indexing='xy')
kA2=np.array([128,0]); kB2=np.array([127,16])
d=N2/np.linalg.norm(kA2); alv=np.arctan2(16,127)
tA=(np.cos(2*np.pi*(kA2[0]*X2+kA2[1]*Y2)/N2)>np.cos(np.pi*0.4)).astype(float)
tB=(np.cos(2*np.pi*(kB2[0]*X2+kB2[1]*Y2)/N2)>np.cos(np.pi*0.4)).astype(float)
F=np.fft.fft2(tA*tB)
fx=np.fft.fftfreq(N2)*N2; FX,FY=np.meshgrid(fx,fx,indexing='xy')
mask=(np.hypot(FX,FY)<0.5*np.linalg.norm(kA2))&(np.hypot(FX,FY)>0)
Fm=np.abs(F)*mask
iy,ix=np.unravel_index(np.argmax(Fm),Fm.shape)
kdom=np.array([FX[iy,ix],FY[iy,ix]])
D_num=N2/np.linalg.norm(kdom)
D_th=d/(2*np.sin(alv/2))
print('モアレ周期 数値',D_num,' 理論',D_th,' (周期のわずかな差を含む)')
assert abs(D_num-D_th)/D_th<5e-3
# 縞の向き：縞（等位相線）は kdom に垂直．x軸となす角
ang=np.degrees(np.arctan2(kdom[1],kdom[0]))+90
ang=(ang+90)%180-90
print('縞がx軸となす角',ang,' 理論 α/2=',np.degrees(alv/2))
assert abs(abs(ang)-np.degrees(alv/2))<0.3
# (5) 平行で周期が d と d(1+ε)：最近接 fλ ε/(d(1+ε)), モアレ周期 d(1+ε)/ε
eps=0.05
print('ノギス型: 周期',d*(1+eps)/eps)
print('OK')

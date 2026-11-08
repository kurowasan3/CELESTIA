# O20193-1 : thickness resonance of a free-free plate driven at depth h
# Independent check: lumped mass-spring chain (N+1 nodes, free-free), harmonic force at the node at depth h,
# light damping, tridiagonal steady-state solve; response amplitude of the bottom-surface node vs frequency.
import numpy as np
from scipy.linalg import solve_banded
d=0.020; VL=6000.0
def peaks(V,hfrac,fmin,fmax,force_harm=1,N=240):
    dz=d/N; rho=1.0; ks=rho*V**2/dz
    mm=np.full(N+1,rho*dz); mm[0]/=2; mm[-1]/=2
    diag=np.full(N+1,2*ks); diag[0]=ks; diag[-1]=ks
    off=np.full(N,-ks)
    F=np.zeros(N+1,complex); F[int(round(hfrac*N))]=1.0
    fs=np.arange(fmin,fmax,0.5e3); amp=np.empty(len(fs))
    for j,f in enumerate(fs):
        w=2*np.pi*f*force_harm
        ab=np.zeros((3,N+1),complex)
        ab[0,1:]=off; ab[2,:-1]=off
        ab[1]=diag-w**2*mm+1j*w*mm*2*np.pi*500.0
        amp[j]=abs(solve_banded((1,1),ab,F)[-1])
    med=np.median(amp)
    pk=[fs[i] for i in range(1,len(fs)-1) if amp[i]>amp[i-1] and amp[i]>=amp[i+1] and amp[i]>8*med]
    return np.round(np.array(pk)/1e3,1)
print('h=d/2 , VL :',peaks(VL,0.5,50e3,1100e3))
print('h=d/4 , VL :',peaks(VL,0.25,100e3,1100e3))
print('h=3d/4, VL :',peaks(VL,0.75,100e3,1100e3))
print('h=d/4 , VL, force at 2f :',peaks(VL,0.25,100e3,1100e3,force_harm=2))
f0=VL/(2*d); print('VL/(2d)=',f0/1e3,'kHz ; allowed k:',[k for k in range(1,16) if abs(np.cos(k*np.pi/4))>1e-9])
print('VL=2d*150kHz =',2*d*150e3)

# (4) 新設問：h と d-h で，k 番目の共振での下面の変位の位相は (-1)^k 倍
def bottom_response(V,hfrac,f,N=240):
    dz=d/N; rho=1.0; ks=rho*V**2/dz
    mm=np.full(N+1,rho*dz); mm[0]/=2; mm[-1]/=2
    diag=np.full(N+1,2*ks); diag[0]=ks; diag[-1]=ks
    off=np.full(N,-ks)
    F=np.zeros(N+1,complex); F[int(round(hfrac*N))]=1.0
    w=2*np.pi*f
    ab=np.zeros((3,N+1),complex); ab[0,1:]=off; ab[2,:-1]=off
    ab[1]=diag-w**2*mm+1j*w*mm*2*np.pi*500.0
    return solve_banded((1,1),ab,F)[-1]
# 離散モデルの共振振動数（固有値）で比較
for k in [1,3,4,5,7]:
    N=240; dz=d/N; ks=VL**2/dz**2
    lam_k=2*np.sqrt(ks)*np.sin(k*np.pi/(2*N))/(2*np.pi)   # 両端自由の鎖の固有振動数
    r1=bottom_response(VL,0.25,lam_k); r3=bottom_response(VL,0.75,lam_k)
    ratio=r3/r1
    print(f'k={k} f={lam_k/1e3:.1f} kHz  resp(3d/4)/resp(d/4) = {ratio.real:+.4f}{ratio.imag:+.4f}i  expected {(-1)**k:+d}')
    assert abs(ratio-(-1)**k)<0.05
# 解析：cos(kπ(d-z)/d) = (-1)^k cos(kπz/d)
zz=np.linspace(0,d,101)
for k in range(1,9):
    assert np.allclose(np.cos(k*np.pi*(d-zz)/d),(-1)**k*np.cos(k*np.pi*zz/d))
print('PHY-WAV-016 (4) OK')

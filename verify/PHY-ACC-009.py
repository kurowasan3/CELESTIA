# J137-1 検証：3電圧計法による抵抗をもつコイルの測定
import sympy as sp, numpy as np
R,VR,Vc,V,f=45,90,100,170,50
I=VR/R; cphi=(V**2-VR**2-Vc**2)/(2*VR*Vc); sphi=np.sqrt(1-cphi**2)
r=Vc*cphi/I; X=Vc*sphi/I; L=X/(2*np.pi*f); Pc=(V**2-VR**2-Vc**2)/(2*R)
print('I',I,'cos',cphi,'r',r,'X',X,'L',L,'Pc',Pc,'I^2 r',I**2*r)
# 複素数で順方向に確認
Z=complex(R+r,X); Ic=V/abs(Z); print('|Z|',abs(Z),'I',Ic,'VR',Ic*R,'Vc',Ic*abs(complex(r,X)))
# 力率と電源の電力
print('power factor',(R+r)/abs(Z),'Psource',V*Ic*(R+r)/abs(Z), 'check',Ic**2*(R+r))
# (4) 直列にCを入れて電流最大
C=1/((2*np.pi*f)*X); Imax=V/(R+r); print('C',C,'Imax',Imax,'Vcoil',Imax*abs(complex(r,X)),'VC',Imax*X,'V(coil+C)',Imax*r)
# 位相差：電源電圧に対するコイル電圧の進み（C挿入後）
print('coil phase lead deg',np.degrees(np.arctan2(X,r)))
# 周波数を変えずに，コイルの代わりに同じ読みを与えるもの（容量性）:3つの読みは同じ
Z2=complex(R+r,-X); print('capacitive same readings:',V/abs(Z2)*R, V/abs(Z2)*abs(complex(r,-X)))
# (5) 周波数 100Hz で電源 170V：VR の読み
for sgn,name in [(2,'coil'),(0.5,'cap')]:
    Xn=X*sgn; Zn=abs(complex(R+r,Xn)); print(name,'VR at 100Hz',V/Zn*R)

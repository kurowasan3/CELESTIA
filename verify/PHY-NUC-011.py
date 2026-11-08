# P85-1: モーズリーの法則（遮蔽）と電子捕獲後のX線・ニュートリノ
import numpy as np
from scipy.optimize import brentq
ER=13.6
E22=4.50e3; E40=15.51e3
# E = k (Z-σ)^2 として2点から σ, k
r=np.sqrt(E40/E22); s=(40-r*22)/(1-r); k=E22/(22-s)**2
print('sigma',s,'k',k,'3/4 ER',0.75*ER)
# σ=1, k=3/4 ER と仮定したときの再現
print([0.75*ER*(Z-1)**2 for Z in (22,40)])
# 6.40 keV の X線
Zx=1+np.sqrt(6400/(0.75*ER)); print('Z emitter',Zx)
Zx2=s+np.sqrt(6400/k); print('Z emitter (fitted)',Zx2)
# ニュートリノ
Q=836.0; Eex=136.5; BK=7.1
Enu=Q-Eex-BK; print('Enu keV',Enu)
Mc2=57*931.5e3  # keV
print('recoil eV',Enu**2/(2*Mc2)*1e3)
# K捕獲の条件 Q-E*>BK ; 最小 Q
print('Qmin for K capture',Eex+BK)
# 遮蔽なしだと誤認する Z
print('Z w/o screening',np.sqrt(6400/(0.75*ER)))

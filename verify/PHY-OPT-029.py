# J090-1: grating with white light 400-700nm: order overlap, inverse determination of d, oblique incidence complete orders
import numpy as np, math
l1,l2=400.,700.
# overlap of order m (red end) and m+1 (violet end)
for m in range(1,5): print('m',m,' overlap',m*l2>(m+1)*l1, ' boundary eq', m*l2==(m+1)*l1)
# inverse: at sin(theta0)=0.6 only 600 and 400 present -> L=d*0.6 ; find d values (multiples) giving exactly these two wavelengths
for j in range(1,5):
    L=1200*j; ks=[k for k in range(1,100) if l1<=L/k<=l2]; print('d',L/0.6,' wavelengths',[L/k for k in ks])
# check with d=2000 that other wavelengths absent, scanning orders
d=2000.
# complete orders, normal incidence
comp=[m for m in range(-10,11) if m!=0 and abs(m*l2/d)<1]; print('normal complete',comp)
alpha=math.radians(30); s=math.sin(alpha)
comp=[m for m in range(-10,11) if m!=0 and abs(s+m*l2/d)<1 and abs(s+m*l1/d)<1]; print('oblique complete',comp,len(comp))
# partial orders at oblique
part=[m for m in range(-10,11) if m!=0 and (abs(s+m*l1/d)<1)!=(abs(s+m*l2/d)<1)]; print('partial',part)
for m in part:
    # wavelength range visible
    lo=max(l1, (-1-s)*d/m if m<0 else l1); hi=min(l2,(1-s)*d/m if m>0 else (-1-s)*d/m)
    print(' m',m,' visible up to', (1-s)*d/m if m>0 else (-1-s)*d/m)
# overlap under oblique incidence: direction sin=s+m*l/d ; overlap m & m+1 iff same as normal
print('normal partial',[m for m in range(-10,11) if m!=0 and (abs(m*l1/d)<1)!=(abs(m*l2/d)<1)])

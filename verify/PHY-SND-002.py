# M67-1: string (fundamental) drives closed air column; identify harmonics from tension ratio
import numpy as np
V=340.0; l=0.50; rho=5.0e-4; D=0.630; S1=80.0; S2=222.0
f1=np.sqrt(S1/rho)/(2*l); f2=np.sqrt(S2/rho)/(2*l)
r=f2/f1; print("f1",f1,"f2",f2,"ratio",r)
for m in range(1,6): print(m,(2*m+1)/(2*m-1))
# identify: f1 = (2m-1)V/(4(D+e)), f2=(2m+1)V/(4(D+e)) -> m=2
m=2
Le=(2*m-1)*V/(4*f1); e=Le-D; print("D+e",Le,"e",e)
# alternative independent: use f2 with 5th harmonic
print("D+e from f2 (exact 5/3 ratio)",5*V/(4*f1*5/3), " using measured S2:",5*V/(4*f2))
lam1=V/f1; print("lam1",lam1,"next depth",D+lam1/2)
lam2=4*Le/5; print("lam2",lam2,"V/f2(exact)",V/(f1*5/3))
nodes=[D-k*lam2/2 for k in range(0,3)]; print("nodes",nodes)
antis=[D-lam2/4-k*lam2/2 for k in range(0,3)]; print("antinodes (depth, neg=outside)",antis)
# check the open-end antinode lies at -e
print("top antinode",D-lam2/4-2*lam2/2, "-e=",-e)

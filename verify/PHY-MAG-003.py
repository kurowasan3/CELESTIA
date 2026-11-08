# N50-1: Hall element (sign, n, surface charge) vs electromagnetic flowmeter
import numpy as np
e=1.6e-19
a,b,c=2.0e-3,5.0e-3,2.0e-4   # width(x), length(y), thickness(z)
I,B,VH=4.0e-3,0.50,25e-3
eps=1.0e-10
# sign: current +y, B +z. For carrier q with velocity v (parallel to current if q>0, antiparallel if q<0)
for q in [+1,-1]:
    v=np.array([0,1,0])*np.sign(q)       # drift direction
    F=q*np.cross(v,[0,0,1])
    print('carrier sign',q,'force dir',F,'-> accumulates on', '+x face' if F[0]>0 else '-x face',
          '; +x face potential', 'higher' if q>0 else 'lower')
n=I*B/(e*c*VH); vd=VH/(B*a)
print('n =',n,' drift v =',vd,' check I = e n a c v =',e*n*a*c*vd)
E=VH/a; sigma=eps*E; N=sigma/e
print('E =',E,'sigma =',sigma,'excess per m^2 =',N,' n*a =',n*a,' ratio =',N/(n*a))
# flowmeter
D,h,Bf,Vf=2.0e-2,1.0e-2,0.20,2.0e-3
u=Vf/(Bf*D); print('u =',u,' flow rate =',u*D*h)
for q in [+1,-1]:
    F=q*np.cross([0,1,0],[0,0,1]); print('ion sign',q,'force',F)

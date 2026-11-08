import numpy as np
e=1.6e-19; I=2.0e-3; B=0.50; h=0.50e-3; w=2.0e-3; s=0.20e-3
Vp,Vm=0.30e-3,2.10e-3
V0=(Vp+Vm)/2; VH=(Vm-Vp)/2
n=I*B/(e*h*VH); rho=V0*w*h/(I*s); mu=1/(n*e*rho)
Bstar=V0/(VH/B)
print('V0',V0,'VH',VH,'n',n,'rho',rho,'mu',mu,'B*',Bstar,' s/(w mu)',s/(w*mu),' tan at 0.5T',mu*B,np.degrees(np.arctan(mu*B)))
# independent forward model: electrons, drift velocity, fields, potential difference P(0,0)-Q(s,w)
vd=I/(n*e*w*h)           # electrons move -x
def reading(Bz):
    Ex=rho*I/(w*h)        # along +x (current direction)
    # force balance on electron in y: (-e)(Ey + (v x B)_y)=0 ; v=(-vd,0,0), B=(0,0,Bz): (v x B)_y = -(-vd*Bz)= vd*Bz
    Ey=-vd*Bz
    phi=lambda x,y: -Ex*x-Ey*y
    return phi(0,0)-phi(s,w)
print('forward model V(+B),V(-B),V(B*)',reading(B),reading(-B),reading(Bstar))

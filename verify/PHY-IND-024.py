# T20222-1: coil cart with battery opposing ideal diode -> current only when BLv>V; energy split
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
# symbolic: exact identity dE_b = V I dt = -(V/BL) m dv  -> E_b = m u (v0-u); heat = KE loss - E_b
m,v0,u=sp.symbols('m v0 u',positive=True)
Eb=m*u*(v0-u); H=sp.Rational(1,2)*m*(v0**2-u**2)-Eb
print('heat =',sp.factor(H))
print('maximize Eb over u:',sp.solve(sp.diff(Eb,u),u),' Eb max=',sp.simplify(Eb.subs(u,v0/2)),' heat=',sp.simplify(H.subs(u,v0/2)))
# small-change per strip: impulse = BL*(BLv-V)/R * d/v
B,L,R,d,V,v=sp.symbols('B L R d V v',positive=True)
print('Delta v per strip =',sp.simplify(B*L*(B*L*v-V)/R*(d/v)/m))
# numeric: many strips, exact ODE (only front side conducts)
Bn,Ln,Rn,dn,mn,Vn,v0n=0.5,0.2,0.05,0.02,0.1,0.2,4.0   # u=V/(BL)=2.0
un=Vn/(Bn*Ln); gap=0.5   # strip period (>L)
def rhs(t,y):
    x,vv,Eb,Q=y
    xs=(x+Ln/2)%gap   # front side position relative to strip start
    I=0.0
    if xs<dn:
        I=max(Bn*Ln*vv-Vn,0)/Rn
    return [vv,-Bn*Ln*I/mn,Vn*I,I*I*Rn]
sol=solve_ivp(rhs,[0,400],[0,v0n,0,0],max_step=0.0005,rtol=1e-9,atol=1e-12)
vf,Ebf,Qf=sol.y[1,-1],sol.y[2,-1],sol.y[3,-1]
print('final v',vf,' u',un)
print('E_b num',Ebf,' formula m u (v0-v_f)',mn*un*(v0n-vf), ' limit m u(v0-u)',mn*un*(v0n-un))
print('heat num',Qf,' KEloss - Eb',0.5*mn*(v0n**2-vf**2)-mn*un*(v0n-vf),' limit',0.5*mn*(v0n-un)**2)
# check per-strip small change formula on first strip
k=Bn**2*Ln**2*dn/(mn*Rn); print('predicted first-strip dv',k*(1-un/v0n))
i=np.argmax(sol.y[0]>gap-Ln/2+dn+1e-3); print('numeric first-strip dv',v0n-sol.y[1,i])

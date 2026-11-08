# N64-1: alpha decay with excited daughter levels; recoil correction
import sympy as sp
u=931.5
mTh,mRa,mHe=228.028741,224.020211,4.002603
Q=(mTh-mRa-mHe)*u; print('Q=',Q)
# split from rest (masses ~ mass numbers) - symbolic check
ma,mY,Qs,va,vY=sp.symbols('m_a m_Y Q v_a v_Y',positive=True)
s=sp.solve([sp.Eq(ma*va,mY*vY),sp.Eq(ma*va**2/2+mY*vY**2/2,Qs)],[va,vY],dict=True)[0]
Ka=sp.simplify(ma*s[va]**2/2); print('K_alpha =',Ka)
K0=Q*224/228; print('predicted main alpha KE',K0)
# with exact masses
print('exact-mass ratio',Q*mRa/(mRa+mHe))
Ks=[5.423,5.340,5.177]
Qeff=[k*228/224 for k in Ks]; print('Q_eff',Qeff)
E1=(Ks[0]-Ks[1])*228/224; E2=(Ks[0]-Ks[2])*228/224
print('E1*,E2* (MeV)',E1,E2,' E2-E1',E2-E1)
print('recoil KE of daughter main',Ks[0]*4/224,' total Q from K0',Qeff[0])
# gamma recoil correction for 84 keV from Ra-224 (Mc^2 = 224*931.5)
Mc2=224*u; print('gamma recoil shift (keV)',1e3*E1**2/(2*Mc2))
# naive error
print('naive ΔK',Ks[0]-Ks[1],' rel err',((Ks[0]-Ks[1])-E1)/E1)

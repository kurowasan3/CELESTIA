# Q166-1 地磁気中を水平飛行する航空機の誘導起電力（ベクトルで直接計算）
import numpy as np
v=250.0; w=60.0; h=12.0; L=50.0
def field(lat_deg,B0):
    lam=np.radians(lat_deg)
    # x: east, y: north, z: up. dipole: horizontal north B0 cos(lam), vertical 2B0 sin(lam) downward in north hemisphere
    return np.array([0.0, B0*np.cos(lam), -2*B0*np.sin(lam)])
def emfs(lat,B0,heading_deg):
    psi=np.radians(heading_deg)     # heading measured from north toward east
    fwd=np.array([np.sin(psi),np.cos(psi),0.0]); up=np.array([0,0,1.0]); left=np.cross(up,fwd)
    F=np.cross(v*fwd,field(lat,B0))  # force per unit charge on positive charges = v x B ; E_eq inside = -(v x B)
    # potential difference phi(left)-phi(right): positive charges pushed toward +F side -> that side higher potential
    return dict(wing_left_minus_right=np.dot(F,left)*w, top_minus_bottom=np.dot(F,up)*h, nose_minus_tail=np.dot(F,fwd)*L)
B0=3.0e-5
for hd in [0,45,90,180,270]:
    print('lat30 heading',hd, {k:round(x,4) for k,x in emfs(30,B0,hd).items()})
# part: southern hemisphere, heading east: right wing higher by 0.40 V, top higher by 0.10 V -> lat, B0
tanl=h*0.40/(2*w*0.10); lam=np.degrees(np.arctan(tanl)); B0s=0.10/(v*h*np.cos(np.radians(lam)))
print('tan lam=',tanl,' lam=',lam,' B0=',B0s)
print('check:',emfs(-lam,B0s,90))
print('dip angle tanI=2tan lam ->',np.degrees(np.arctan(2*tanl)))

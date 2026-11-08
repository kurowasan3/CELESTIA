# Q193-1 正方形の頂点の陽電子2個・陽子2個：放す順番と最終速度（質量比による時間スケール分離）
import numpy as np
from scipy.integrate import solve_ivp
k=e=a=1.0; m=1.0; M=1836.15
q=np.ones(4); mass=np.array([m,M,m,M])   # 0:A(e+),1:B(p),2:C(e+),3:D(p)
pos0=np.array([[0,a],[a,a],[a,0],[0,0]],float)
def rhs(t,s,free):
    x=s[:8].reshape(4,2); v=s[8:].reshape(4,2); acc=np.zeros((4,2))
    for i in range(4):
        if not free[i]: continue
        for j in range(4):
            if i!=j:
                d=x[i]-x[j]; acc[i]+=k*q[i]*q[j]*d/np.linalg.norm(d)**3/mass[i]
    vv=v*np.array(free)[:,None]
    return np.concatenate([vv.ravel(),acc.ravel()])
def U(x):
    s=0
    for i in range(4):
        for j in range(i+1,4): s+=k/np.linalg.norm(x[i]-x[j])
    return s
U0=U(pos0); print('U0',U0,' theory',(4+1/np.sqrt(2)))
th_e=np.sqrt(2*(2+np.sqrt(2)/4)/m); th_p=np.sqrt(2*(np.sqrt(2)/4)/M)
# シナリオ1：同時に放す（長時間積分）
s=solve_ivp(rhs,[0,4e5],np.concatenate([pos0.ravel(),np.zeros(8)]),args=([1,1,1,1],),rtol=1e-10,atol=1e-12,method='DOP853')
x=s.y[:8,-1].reshape(4,2); v=s.y[8:,-1].reshape(4,2)
KE=0.5*mass*np.sum(v**2,1); print('sim1 KE e+',KE[0],KE[2],' p',KE[1],KE[3],' remaining U',U(x))
ve=np.linalg.norm(v[0]); vp=np.linalg.norm(v[1])
print(' v_e/theory',ve/th_e,' v_p/theory',vp/th_p, ' sqrt(m/M)=',np.sqrt(m/M))
# シナリオ2：陽電子を固定して陽子だけ放す → 十分遠くへ → 陽電子を放す
s2=solve_ivp(rhs,[0,4e5],np.concatenate([pos0.ravel(),np.zeros(8)]),args=([0,1,0,1],),rtol=1e-10,atol=1e-12,method='DOP853')
st=s2.y[:,-1].copy()
s3=solve_ivp(rhs,[0,4e3],st,args=([1,1,1,1],),rtol=1e-10,atol=1e-12,method='DOP853')
v=s3.y[8:,-1].reshape(4,2); x=s3.y[:8,-1].reshape(4,2); KE=0.5*mass*np.sum(v**2,1)
print('sim2 KE e+',KE[0],' p',KE[1],' remaining U',U(x))
print(' theory: e+ each',np.sqrt(2)/4,' p each',2+np.sqrt(2)/4)

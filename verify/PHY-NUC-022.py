# S0604-1 hydrogen-like ion series overlap, counting, identification
import sympy as sp, math
m,n=sp.symbols('m n',positive=True)
# (1) overlap condition: longest wl of series m+1 is inf? No: series m+1 occupies 1/lam in (Z^2R(1/(m+1)^2-1/(m+2)^2), Z^2R/(m+1)^2]
# series m occupies 1/lam in [Z^2R(1/m^2-1/(m+1)^2), Z^2R/m^2). Overlap iff 1/(m+1)^2 > 1/m^2-1/(m+1)^2
for M in range(1,8):
    ov = 1/(M+1)**2 > 1/M**2-1/(M+1)**2
    print(M, ov)
print(sp.solve(sp.Eq(2*m**2,(m+1)**2),m))
# (2) count lines of series m with wavelength longer than limit of series m+1
def count(M):
    return [N for N in range(M+1,200) if 1/M**2-1/N**2 < 1/(M+1)**2]
for M in [3,4,5,6]:
    c=count(M); print(M,c,len(c), M*(M+1)/math.sqrt(2*M+1))
# (3) data
l1=91.18 # nm, hydrogen Lyman limit hc/13.6eV
hc=1239.84; print("hc/13.6=",hc/13.6)
lim1=hc/13.6
lp,lq=51.3,91.2
print("ratio",lq/lp, 16/9)
Z2=9*lim1/lp; print("Z^2 (m=3)",Z2)
Z2b=36*lim1/lp; print("Z^2 if (6,8)",Z2b, " series7 limit",49*lim1/64)
# all series limits for Z=4
for M in range(1,7): print("limit",M, M*M*lim1/16)
# (4) identify line 80.1 nm
Z=4
for M in range(1,6):
    for N in range(M+1,40):
        lam=1/(Z*Z/lim1*(1/M**2-1/N**2))
        if abs(lam-80.1)<1.0: print("match",M,N,lam)
# ionization energy from n=2 of Z=4
print("E_ion(n=2)",16*13.6/4)
print("==== final version uses Z=6, limits 40.5 & 63.3 nm, line 72.9 nm")
lp,lq=40.5,63.3
print("ratio",lq/lp,25/16)
print("Z^2 (m=4)",16*lim1/lp,"Z^2 if (8,10)",64*lim1/lp,"series9 limit if Z=12",81*lim1/144)
Z=6
for M in range(1,8): print("limit",M, M*M*lim1/Z**2)
for M in range(1,9):
    for N in range(M+1,400):
        lam=1/(Z*Z/lim1*(1/M**2-1/N**2))
        if abs(lam-72.9)<0.3: print("match",M,N,lam)
# n continuous for series 5 at 72.9
x=0.04-lim1/(36*72.9); print("series5 needs 1/n^2=",x,"n=",1/math.sqrt(x))
print("series4 lines:",[round(1/(Z*Z/lim1*(1/16-1/N**2)),2) for N in range(5,9)])
for M in range(1,9):
    for N in range(M+1,400):
        lam=1/(Z*Z/lim1*(1/M**2-1/N**2))
        if abs(lam-52.1)<0.3: print("match52",M,N,lam)

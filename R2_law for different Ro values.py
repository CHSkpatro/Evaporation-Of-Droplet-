import numpy as np
import matplotlib.pyplot as mpt
r1=.8
r2=.85
r3=.90
r4=.95
r5=1
t=np.arange(0,1500,10)
k=2*3.63312*10**-4
mpt.subplot(1,1,1)
mpt.plot(t,.64- k*t,'r',label=r1)
mpt.plot(t,(r2**2)-k*t,'g',label=r2)
mpt.plot(t,(r3**2)-k*t,'b',label=r3)
mpt.plot(t,(r4**2)-k*t,'m',label=r4)
mpt.plot(t,(r5**2)-k*t,'c',label=r5)
mpt.ylim(bottom=0)
mpt.xlim(left=0,right=1500)
mpt.legend()
mpt.grid(True)
mpt.xlabel('Time t in seconds')
mpt.ylabel('$R^2$(t) in $mm^2$')
mpt.title('$R^2$(t) Vs t for different Ro(mm) values , for H=0.4')
mpt.show()

import numpy as np
import matplotlib.pyplot as mpt
r=float(input('The value of initial radius in mm is='))
t=np.arange(0,5000,10)
H1=0.2
H2=0.4
H3=0.6
H4=0.8
K1=2*(1-H1)*6.06*10**-4
K2=2*(1-H2)*6.06*10**-4
K3=2*(1-H3)*6.06*10**-4
K4=2*(1-H4)*6.06*10**-4
mpt.figure(figsize=(8,6))
mpt.subplot(1,1,1)
mpt.plot(t,r**2- K1*t,'r',label='H=0.2')
mpt.plot(t,(r**2)-K2*t,'g',label='H=0.4')
mpt.plot(t,(r**2)-K3*t,'b',label='H=0.6')
mpt.plot(t,(r**2)-K4*t,'m',label='H=0.8')

mpt.ylim(bottom=0)
ax = mpt.gca()
ax.set_xlim([0, 5000])
mpt.legend()
mpt.grid(True)
mpt.xlabel('Time t in seconds')
mpt.ylabel('$R^2$(t) in $mm^2$')
mpt.title('$R^2$(t) Vs t for different H values , for Ro=%.2f mm' %r)
mpt.show()

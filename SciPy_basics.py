from scipy import cluster
#help(cluster)
#import scipy
#scipy.info(cluster)
#scipy.source(cluster)
######################################  SPECIAL FUNCTIONS :EXPONENTIAL AND TRIGNOMETRIC ################################
from scipy import special
y=special.exp2(45) #only exp10 and exp2 are present 
print(y)
z=special.sindg(49)
print(z)
a=special.cosdg(67)
print(a)
############ Integration Functions ##############
#### quad function calculates the integration of a function with respect to a single variable where as the dblquad function calculates 
####              double integral of the function with respect to two variables ######
from scipy import integrate # in integrate the quad feature is available
i=integrate.quad(lambda x:special.exp10(x),0,2)# in python lamn=bda can be used for a f(x)
print(i)
e=lambda x,y:x*y**2
f=lambda x: 1
g=lambda x:-1
f=integrate.dblquad(e,0,2,gfun=f,hfun=g)# in double integral here limits of y are defined by the help of functions of x, f=lambda(x=1),g=lambda(x=-1)
print(f)#the result we will see that contain two values 1st is the approximate result and second is the approximate error in the result

#################   FOURIER TRANSFORMATION ########
### fft and ifft for inverse fourrier transform ####
# the fft pack 
from scipy.fftpack import fft,ifft
import numpy as np
list=[1,2,23,4]
y=fft(np.array(list))
print(y)
z=ifft(np.array(list))
print(z)
############# Linear Algebra #########
#### ATLAS lAPACK and BLAS libraries #######
### calculating inverse of a matrix###inv pack
import scipy.linalg as lin
a=np.array([[1,2],[3,4]])
b=lin.inv(a)
print(b)
############# Interpolation Functions ####
####SCIPY.INTERPOLATE , CREATING A SET OF NEW DATA POINTS WITHIN KNOWN DATA POINTS ###WE WANT TO PLOT IT SO ALSO IMPORT MATPLOTLIB.PYPLOT
import matplotlib.pyplot as mpt
import scipy.interpolate as int
x=np.arange(5,30)
y=np.exp(x/3)
f=int.interp1d(x,y)
x1=np.arange(6,12)
y1=f(x1)
mpt.plot(x,y,'o',x1,y1,'--')
mpt.show()






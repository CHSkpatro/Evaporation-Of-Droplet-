import numpy as np
import matplotlib.pyplot as mpt
import pandas as pd
r_min=float(input('The minimun radius value in mm: '))
r_max=float(input('The maximum radius value in mm: '))
N=int(input('No. of divisions of r axis : '))
del_r=(r_max-r_min)/N
Z_min=float(input('The minimun z value in mm: '))
Z_max=float(input('The maximum z value in mm: '))
M=round((Z_max-Z_min)/(del_r))
print(M)
C1=float(input('Concentration value for r=r_min : '))
C2=float(input('Concentration value for Z=Z_min : '))
C3=float(input('Concentration value for r=r_max : '))
C4=float(input('Concentration value for Z=Z_max : '))
##print(del_r)
##print((Z_max-Z_min)/M)
r=np.linspace(r_min,r_max,N+1)
##print(r)
z=np.linspace(Z_min,Z_max,M+1)
###  A*C=  B ####
A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
###filling up diagonal elements###
for i in range ((N-1)*(M-1)):
    A[i,i]=round(-(4 +(del_r/r[(i%(N-1)) + 1])),1)
    if i< (((N-1)*(M-1))-1):  
        if i==0 or (i+1)%(N-1)!=0:
          A[i,i+1]=round(1+(del_r/r[(i%(N-1)) + 1]),1)###super diagonal elements###but there are still slective rows in which the super diagonal element is 0
        if i==0 or (i+1)%(N-1)!=0:
          A[i+1,i]=1###sub diagonal elements ####but also there are selective rows in which the sub diagonal element is 0 
for i in range((N-1),((N-1)*(M-1))):
    A[i,i-(N-1)]=1  ###filling the second sub diagonal elements###
    A[i-(N-1),i] =1 ###filling up the 2nd super diagonal elements###
df = pd.DataFrame(A, index=range(1, (N - 1) * (M - 1) + 1), columns=range(1, (N - 1) * (M - 1) + 1))
df.to_csv('matrix_A_cyl_PPT.csv', index=True)
print("Matrix A cyl PPT has been saved to 'matrix_A.csv'.")  
import numpy as np
import matplotlib.pyplot as plt
##A.B=C ###  A is  a square matrix of NXN , B and C are column vectors of N rows , To find is B while  A and C are knowns 
A=np.array([  [1,2,0,0],
              [2,1,2,0], 
              [0,2,1,2],
              [0,0,2,1]  ])
A_inv=np.linalg.inv(A)
C=np.zeros((4,1))
C[0][0]=1
C[3][0]=0.4
##print(C)
B=np.dot(A_inv,C)
print(B) 
x=np.linspace(0,1,4)
plt.subplot(1,2,1)
plt.scatter(x,B.flatten(), marker='o', color='b', label='Matrix B values')
plt.xlabel('Index')
plt.ylabel('Matrix B Values')
plt.title('Scatter Plot of Matrix B Values')
plt.grid(True)
plt.legend()
plt.subplot(1,2,2)
plt.plot(x,B.flatten(),'-')
plt.xlabel('Index')
plt.ylabel('Matrix B Values')
plt.title('curve joining  of Matrix B Values')
plt.grid(True)
plt.legend()
plt.show()
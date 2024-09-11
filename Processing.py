import numpy as np
import matplotlib.pyplot as mpt
from PIL import Image 
pic=Image.open('Photo.jpg')
pic.show()
picarray=np.asarray(pic)
print(picarray.shape)
#mpt.imshow(picarray)
#mpt.show()
A=picarray.copy()
picR=A[:,:,0]
print(picR)


mpt.imshow(picR,cmap='gray')
mpt.show()
picG=A[:,:,1]
print(picG)
mpt.imshow(picG,cmap='gray')
mpt.show()

picB=A[:,:,2]
print(picB)
mpt.imshow(picB,cmap='gray')
mpt.show()

A[:,:,2]=0
B=A.copy()
mpt.imshow(B)
mpt.show()

A[:,:,1]=0
C=A.copy()
mpt.imshow(C)
mpt.show()

A[:,:,0]=0
mpt.imshow(A)
mpt.show()
import numpy as np
import matplotlib.pyplot as mpt
import cv2

img=cv2.imread('image.jpg')#the image is saved as numpy.ndarray type directly
print(type(img))  #if type printed as none there is some erroer in directory name
print(img.shape)
mpt.imshow(img)
mpt.show()  #the image being shown will be bluish and will be different from the original imaage as open cv and matplotlib read the color channels differently

fix_img=cv2.cvtColor(img,cv2.COLOR_BGR2RGB) #Matplotlib reads color chanels as RGB and open CV as BGR, so we converted the color channel from BGR to RGB as at last we need in RGB 
mpt.imshow(fix_img)
mpt.show()
#now we want to convert into grayscale image
img_gray=cv2.imread('image.jpg',cv2.IMREAD_GRAYSCALE)
mpt.imshow(img_gray) 
mpt.show() #but it still wont show gray image even though the image is converted into 2D array which can be checked from shape() for gray image 
mpt.imshow(img_gray,cmap='gray') 
mpt.show()#now gray image is shown 
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
#RESIZING AND FLIPPING OF IMAGES#
#LET WE WANT TO RESIZE THE IMG_GRAY#
print(img_gray.shape) #(956,1300)956 is height and 1300 is width#
new_img=cv2.resize(img_gray,(300,500))#here it means the 956 pixel is converted into 500 pixel and 1300 pixel is converted into 1300 pixel is converted into 300#
mpt.imshow(new_img)
mpt.show()
#now let we want to resize the image into a ratio#
w_ratio=.67
h_ratio=.45#means  we want to convert the image width to 67% and image height to 45%#
rat_img=cv2.resize(img_gray,(0,0),img_gray,w_ratio,h_ratio)
mpt.imshow(rat_img)
mpt.show()
print(rat_img.shape)
#Flipping images#
flip_image=cv2.flip(rat_img,0)#the image is flipped about horizontal axis#
mpt.imshow(flip_image)
mpt.show()

flip_img1=cv2.flip(rat_img,1)#the image is flipped about middle vertical axis$
mpt.imshow(flip_img1)
mpt.show()

flip_img=cv2.flip(rat_img,-1)#the image is flipped 1st about horizonatl and then about vertical axis#
mpt.imshow(flip_img)
mpt.show()
#Saving the image#
cv2.imwrite('New_image.jpg',img_gray)
#How to enlarge the displayment of an image without resizing it#
fig_x=mpt.figure(figsize=(10,8))#10 inches and 8 inches canvas size on which the image will be displayed
ax=fig_x.add_subplot(111)
ax.imshow(fix_img)#the image name which we want to be enlarged or downsized#
mpt.show()
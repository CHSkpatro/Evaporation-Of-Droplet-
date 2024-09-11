import cv2
img=cv2.imread('image.jpg')
while True:
    cv2.imshow('Beautiful',img)
    if cv2.waitKey(10000) & 0xFF==27:#means if image is open for 10000millisecond and escape key(27) is pressed then the while loop will break and all windows will b e closed by the line cv2.destroAllWindows
        break

cv2.destroyAllWindows
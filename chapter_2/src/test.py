import numpy as np
x = np.pi/3
test1 = np.array([[np.cos(x),-np.sin(x),0],
                  [np.sin(x),np.cos(x),0],
                  [0,0,1]])
test2 = np.array([[1,0,0],
                  [0,np.cos(x),-np.sin(x)],
                  [0,np.sin(x),np.cos(x)]])
test3 = np.array([[np.cos(x),0,np.sin(x)],
                  [0,1,0],
                  [-np.sin(x),0,np.cos(x)]])


def is_rotation_valid(r):
    r12 = np.dot(r[:, 0], r[:, 1])
    r23 = np.dot(r[:, 1], r[:, 2])
    r31 = np.dot(r[:, 2], r[:, 0])
    det = np.linalg.det(r)

    if r12 == 0 and r23 == 0 and r31 == 0 and det == 1:
        return 1

n = is_rotation_valid(test1)
s = is_rotation_valid(test2)
a = is_rotation_valid(test3)
print(a)
print(s)
print(n)
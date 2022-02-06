import numpy as np
t=[4,2,1,3,5]
B = np.array([1,2,3])

m_list = [[t[2],t[1],t[0]],[t[3],t[2],t[1]],[t[4],t[3],t[2]]]
A = np.array(m_list)
print(A)
X = np.linalg.inv(A).dot(B)
print(X)
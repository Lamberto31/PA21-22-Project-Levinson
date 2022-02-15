import numpy as np
#t = [4,2,1,3,5]
t = [6,4,2,1,3,5,7]
B = np.array([1,2,3,4])
n = 4
m_list = []
for i in range(n):
    rows = []
    for j in range(n):
        rows.append(t[i-j+n-1])
    m_list.append(rows)

A = np.array(m_list)
print(A)
X = np.linalg.inv(A).dot(B)
print(X)
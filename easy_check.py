import numpy as np
#t = [4,2,1,3,5]
#t = [6,4,2,1,3,5,7]
#y = [1,23,,4]
t = [7,2,5,9,7,4,4,4,7,4,7,5,4,8,8,9,9,7,5]
y = [4,2,10,5,9,7,10,4,1,9]

B = np.array(y)
n = len(y)
m_list = []
for i in range(n):
    rows = []
    for j in range(n):
        rows.append(t[i-j+n-1])
    m_list.append(rows)

A = np.array(m_list)
print(A)
print(t)
print(B)
X = np.linalg.inv(A).dot(B)
print(X)
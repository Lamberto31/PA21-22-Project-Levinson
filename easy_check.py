import numpy as np
m_list = [[1,2,4,6],[3,1,2,4],[5,3,1,2],[7,5,3,1]]
B = np.array([1,2,3,4])
A = np.array(m_list)
X = np.linalg.inv(A).dot(B)
print(X)
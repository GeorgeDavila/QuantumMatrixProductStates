import numpy as np
import numpy.linalg as la

import mpnum as mp




CZ = np.array([[ 1.,  0.,  0.,  0.],
               [ 0.,  1.,  0.,  0.],
               [ 0.,  0.,  1.,  0.],
               [ 0.,  0.,  0., -1.]])
               
               





CZ_arr = CZ.reshape((2, 2, 2, 2))




print CZ_arr
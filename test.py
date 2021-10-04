import numpy as np
from fast_marching import fast_marching


if __name__ == '__main__':

    geo_mat = fast_marching.fast_marching(data_dir='./cat1.off', verbose=True)
    geo_mat = np.asarray(geo_mat)
    np.save('res', geo_mat)

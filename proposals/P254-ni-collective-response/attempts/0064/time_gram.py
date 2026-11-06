"""Ordinary all-mode bilinear time Gram; independent right leg, no adjoint."""
import numpy as np


def time_gram(left,right):
    """L[t,i,h,a] R[h,a,s,j] -> output[t,s,i,j], with no rank reduction.

    The owning banks use these contiguous layouts, so both matrix reshapes are
    views. The only dense allocation is the required time/output product.
    """
    nt,m,nh,modes=left.shape
    if right.shape!=(nh,modes,nt,m):
        raise ValueError('matching independent time/output/inner-mode axes required')
    matrix=left.reshape(nt*m,nh*modes)@right.reshape(nh*modes,nt*m)
    return matrix.reshape(nt,m,nt,m).transpose(0,2,1,3)

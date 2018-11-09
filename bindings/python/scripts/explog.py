#
# Copyright (c) 2015 CNRS
# Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
#

import math

import numpy as np

from . import libpinocchio_pywrap as se3


def exp(x):
    if isinstance(x, se3.Motion):
        return se3.exp6FromMotion(x)
    if np.isscalar(x):
        return math.exp(x)
    if isinstance(x, np.ndarray):
        if x.shape == (6, 1):
            return se3.exp6FromVector(x)
        if x.shape == (3, 1):
            return se3.exp3(x)
        raise ValueError('Error only 3 and 6 vectors are allowed.')
    raise ValueError('Error exp is only defined for real, vector3, vector6 and se3.Motion objects.')


def log(x):
    if isinstance(x, se3.SE3):
        return se3.log6FromSE3(x)
    if np.isscalar(x):
        return math.log(x)
    if isinstance(x, np.ndarray):
        if x.shape == (4, 4):
            return se3.log6FromMatrix(x)
        if x.shape == (3, 3):
            return se3.log3(x)
        raise ValueError('Error only 3 and 4 matrices are allowed.')
    raise ValueError('Error log is only defined for real, matrix3, matrix4 and se3.SE3 objects.')

__all__ = ['exp', 'log']

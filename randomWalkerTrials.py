#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:39:13 2019

@author: eg
"""

from skimage.segmentation import random_walker
import skimage.external.tifffile as tiffy
from skimage.filters import threshold_otsu
import numpy as np
import matplotlib.pyplot as plt

fglm = tiffy.imread('FilteredBox4B.tiff')
bi = tiffy.imread('Bin.tiff')
mrk = tiffy.imread('Markers.tiff')
seg = tiffy.imread('watershedSegmentation.tiff')

seg[np.where(seg==0)]=-1
seg[np.where(seg>115)]=-1
seg[np.where(seg<109)]=-1
seg[np.where(seg==109)]=0
seg[np.where(seg==110)]=-1
seg[np.where(seg==111)]=-1
seg[np.where(seg==112)]=-1
seg[np.where(seg==113)]=-1
seg[np.where(seg==114)]=-1
seg[np.where(seg==115)]=0
tiffy.imsave('trialSeg.tiff', seg)
mrk1=np.zeros_like(mrk)
mrk2=np.zeros_like(mrk)
mrk1[np.where(mrk==109)]=1
mrk2[np.where(mrk==115)]=2
seg2=seg+mrk1+mrk2
tiffy.imsave('trialSeg2.tiff', seg)


RW_Prob = random_walker(bi,seg2,return_full_prob=True)
RW_Prob[np.where(RW_Prob==-1)]=0 
tiffy.imsave('Prob0.tiff', RW_Prob[0])
tiffy.imsave('Prob1.tiff', RW_Prob[1])

RW_Seg = random_walker(bi,seg2,return_full_prob=False)
RW_Seg[np.where(RW_Seg==-1)]=0
tiffy.imsave('Segmentation2Particles.tiff', RW_Seg)

'''
RW
    Parameters
    ----------
    data : array_like
        Image to be segmented in phases. Gray-level `data` can be two- or
        three-dimensional; multichannel data can be three- or four-
        dimensional (multichannel=True) with the highest dimension denoting
        channels. Data spacing is assumed isotropic unless the `spacing`
        keyword argument is used.

    labels : array of ints, of same shape as `data` without channels dimension
            Array of seed markers labeled with different positive integers
            for different phases. 
            
            Zero-labeled pixels are unlabeled pixels.
            
            Negative labels correspond to inactive pixels that are not taken
            into account (they are removed from the graph). 
            
            If labels are not consecutive integers, the labels array will be 
            transformed so that labels are consecutive. In the multichannel case,
            `labels` should have the same shape as a single channel of `data`, 
            i.e. without the final dimension denoting channels.
            
    beta : float, optional
        Penalization coefficient for the random walker motion
        (the greater `beta`, the more difficult the diffusion).
    mode : string, available options {'cg_mg', 'cg', 'bf'}
        Mode for solving the linear system in the random walker algorithm.
        If no preference given, automatically attempt to use the fastest
        option available ('cg_mg' from pyamg >> 'cg' with UMFPACK > 'bf').
        - 'bf' (brute force): an LU factorization of the Laplacian is
          computed. This is fast for small images (<1024x1024), but very slow
          and memory-intensive for large images (e.g., 3-D volumes).
        - 'cg' (conjugate gradient): the linear system is solved iteratively
          using the Conjugate Gradient method from scipy.sparse.linalg. This is
          less memory-consuming than the brute force method for large images,
          but it is quite slow.
        - 'cg_mg' (conjugate gradient with multigrid preconditioner): a
          preconditioner is computed using a multigrid solver, then the
          solution is computed with the Conjugate Gradient method.  This mode
          requires that the pyamg module (http://pyamg.github.io/) is
          installed. For images of size > 512x512, this is the recommended
          (fastest) mode.
    tol : float, optional
        tolerance to achieve when solving the linear system, in
        cg' and 'cg_mg' modes.
    copy : bool, optional
        If copy is False, the `labels` array will be overwritten with
        the result of the segmentation. Use copy=False if you want to
        save on memory.
    multichannel : bool, optional
        If True, input data is parsed as multichannel data (see 'data' above
        for proper input format in this case).
    return_full_prob : bool, optional
        If True, the probability that a pixel belongs to each of the labels
        will be returned, instead of only the most likely label.
    spacing : iterable of floats, optional
        Spacing between voxels in each spatial dimension. If `None`, then
        the spacing between pixels/voxels in each dimension is assumed 1.
        
    Returns
    -------
    output : ndarray
        * If `return_full_prob` is False, array of ints of same shape as
          `data`, in which each pixel has been labeled according to the marker
          that reached the pixel first by anisotropic diffusion.
        * If `return_full_prob` is True, array of floats of shape
          `(nlabels, data.shape)`. `output[label_nb, i, j]` is the probability
          that label `label_nb` reaches the pixel `(i, j)` first.
'''
"""

:Functions:
  - grid_eval_1d: take single patch of 1d data and evaluate on new grid
  - grid_output_1d: take 1d (AMR) solution and evaluate on new grid
 
:Todo:
  - extend to 1d and 3d.

"""

from __future__ import print_function
import numpy as np

def grid_eval_1d(X, Q, xout, return_ma=True):
    """
    Utility function that takes a single patch of data in 1d and
    returns values on 1d grid specified by xout.

    Input:
        array X defining a grid patch and data Q on this patch,
        xout defining the points for output (1d array)
        return_ma (bool) determines if output is a masked array
    Returns:
        qout
        
    ndim(Q) is either 1 or 2.  If 2, then Q[m,i] is  m'th variable at i
    
    if return_ma==True then the result is masked at points outside
        the limits of X.   Otherwise result is NaN at these points.
    Uses zero-order interpolation, i.e.
        Sets value qout[i] to value in the finite volume grid cell
        of X that contains xout[i].
    Future: allow linear interpolation instead but this requires
        ghost cell values around Q grid.  (These are present in binary
        output fort.b files but normally thrown away.)
    """

    from scipy.interpolate import interp1d
    from numpy import ma  # for masked arrays
    
    Qdim = Q.ndim
    if Qdim == 1:
        # change to 2d array of shape (1, len(Q)]):
        Q = np.array([Q])
        
    nvars = Q.shape[0]  # number of arrays to interpolate

    dx = X[1] - X[0]

    # augment Q with border of values:
    x1 = np.hstack((X[0]-0.501*dx, X, X[-1]+0.501*dx))
    Q1 = np.empty((nvars,len(x1)))
    Q1[:,1:-1] = Q   # center portion
    Q1[:,0] = Q[:,0]
    Q1[:,-1] = Q[:,-1]

    qout = np.empty((nvars, len(xout)))
    for k in range(nvars):
        evalfunc = interp1d(x1, Q1[k,:], kind='nearest',
                bounds_error=False, fill_value=np.nan)
        #import pdb; pdb.set_trace()

        qout[k,:] = evalfunc(xout)
        if nvars == 1:
            qout = qout[0,:]

    if return_ma:
        # convert from an array with nan's to a masked array:
        qout = ma.masked_where(qout != qout, qout)

    #print('type is %s' % type(qout))

    return qout
    

def grid_output_1d(framesoln, out_var, xout, levels='all', 
                   return_ma=True, return_level=False):

    """
    :Input:
        framesoln:  One frame of Clawpack solution (perhaps with AMR),
                 An object of type pyclaw.Solution.solution.
        out_var: function that maps q to desired quantities Q[m,i,j] or
                 Q[i,j] if only one.  
                 If type(out_var) == int, then Q[i,j] = q[out_var,i,j]
        xout, yout: arrays of output points (1d or 1d arrays)
        levels: list of levels to use, or 'all'
        return_ma: True to return as masked_array, False to return with
                NaN in locations that framesoln doesn't cover.
        return_level: If False, just returns qout
                      If True, returns qout,level where level[i] is the AMR 
                         level of the patch used to fill qout[i].
    :Output:
        qout: Solution obtained on xout grid (if return_level==False)
        qout, level: if return_level==True

    Loop over all patches in framesoln and apply grid_eval function.
    Use non-NaN values that this returns to update qout array over the
        region covered by this patch

    :Example:

    Note that one frame of a Clawpack simulation can be loaded via, e.g.:

        from clawpack.pyclaw.solution import Solution
        framesoln = Solution(frameno=1, path='_output', 
                             file_format='ascii')

    Then define `xout, yout, out_var` and call this function.

    """
        
    from numpy import ma  # for masked arrays
    if levels == 'all':
        levels = range(1,100)  # more levels than will ever use

    qout = np.empty(xout.shape)
    qout[:] = np.nan
    if return_level:
        level = np.empty(xout.shape)
        level[:] = np.nan
    xmin = xout.min()
    xmax = xout.max()
    for stateno,state in enumerate(framesoln.states):
        state = framesoln.states[stateno]
        patch = state.patch
        #print('level = ',patch.level)

        if patch.level not in levels:
            # skip this patch
            continue

        if (xmin > state.grid.x.upper) or (xmax < state.grid.x.lower):
            # no overlap
            continue
            
        #print('overlap at level %i' % patch.level)
        xc = state.grid.c_centers[0]

        if type(out_var) == int:
            Q = state.q[out_var, :]
        else:
            Q = out_var(state.q)

        qout1 = grid_eval_1d(xc, Q, xout, return_ma=False)
        qout = np.where(np.isnan(qout1), qout, qout1)
        if return_level:
            level = np.where(np.isnan(qout1), level, patch.level)
            
    if return_ma:
        # convert from an array with nan's to a masked array:
        qout = ma.masked_where(qout != qout, qout)
            
    if return_level:
        return qout, level
    else:
        return qout


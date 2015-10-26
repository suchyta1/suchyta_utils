import _jkfunctions as _jk


def JKSphere(jarrs, jras, jdecs, jfunc, jargs=[], jkwargs={}, jtype='generate', jfile=None, njack=24, generateonly=False, gindex=0, varonly=False, save=None, rtype='n-1', itsave=False):
    """
    The first argument (jarrs) is an array of the tables (recarrays) you want to JK.
    The second argument (jras) is a list of the ra columns in these tables.
    The third argument (jdecs) is a list of the dec columns in these tables.
    The fourth agrument (jfunc) is the target function you want to JK.

    Optional Args:
        jargs: extra arguments to the target function, other than the required first first agrument, see comments above.
        jkwargs: keyword extra arguments to the target function

        jtype: 'generate' will use kmeans to generate the JK regions. 'read' reads them from a file
        jfile: used with jtype='read', the file to read from
        njack: number of JK regions, only relevant when generating, otherwise the file sets the number
        generateonly: returns just the name of the JK region file after generating it, without actually doing any JKing
        gindex: which element in jarrs/jras/jdecs to use for generating the JK regions.

        save: If given, a directory to save output
            numbering scheme is indexed, by the index it is in the return array
            cov/ has covariance matrix (i.e. Element 1 below)
            vec/ has results of each thing you return (i.e. Element 0 below)
            other/ has results of other things you return (i.e. Element 2 below)
        itsave: If True, save information about each JK iteration
            it/ will appear under vec/ and other/, with as many files as JK iterations


    Returns a 4-element array:
        Element 0: Array of the answers over the full area of the things you JKed
        Element 1: Covariances of the things you JKed
        Element 2: Array of whatever extra stuff you returned from your target function during the call for the full area
        Element 3: Like Element 2, but with arrays in each component for the extra returns of the target function in each JK realization.


    Look in _jkfunctions for some dummy examples. I recommend playing around with this function
    to get a better feel for what the arguments/returns mean. Once you understand it, you can
    get a lot of mileage out of this function.
    """
    jarrs, jras, jdec = _jk.ArrFunc(_jk.EnforceArray2D, jarrs, jras, jdecs)
    km, jfile = _jk.GenerateRegions(jarrs, jras, jdecs, jfile, njack, gindex, jtype)
    if generateonly:
        return jfile

    ind = _jk.AssignIndex(jarrs, jras, jdecs, km)
    save, itsave = _jk.SetupSave(save, itsave)
    full_j, full_other = _jk.FullResults(jfunc, jargs, jkwargs, jarrs) 
    it_j, it_other, frac = _jk.Iterate(jfunc, jargs, jkwargs, jarrs, ind, njack, full_j, full_other, rtype, itsave)
    cov_j = _jk.GetCov(full_j, it_j, njack, varonly, rtype)

    if itsave is not None:
        _jk.ItSave(itsave, [it_j, it_other])

    if save is not None:
        _jk.FullSave(save, [full_j, cov_j, full_other] )

    return [full_j, cov_j, full_other, it_other]

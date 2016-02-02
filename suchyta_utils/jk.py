"""
Since we often don't know what else to do with errors, astronomers like to get `jackknife <https://en.wikipedia.org/wiki/Jackknife_resampling>`_ (JK) covariances.
Many of us (including me) tend to work in angular coordinates. We get a bunch of RA and DECs -- coordinates on a unit sphere.
:mod:`suchyta_utils.jk` provides a way to generate JK realizations for an arbitrary function, which takes data on the sphere as an arugment.
The implementation uses `k-means clustering <https://en.wikipedia.org/wiki/K-means_clustering>`_ for region generation,
using Erin Sheldon's `kmeans_radec module <https://github.com/esheldon/kmeans_radec>`_.
:mod:`suchyta_utils.jk` provides functionality to use these regions for JK operations of datasets.

I know that I refer to examples that don't exist anywhere in the documenation yet. I'll fix this when I find the time.

"""

import _jkfunctions as _jk
import numpy as _np
import os
import kmeans_radec


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


class SphericalJK(object):
    """
    Class for executing JK realizations of an arbitrary target function.
    For N regions, calling the :func:`~suchyta_utils.jk.SphericalJK.DoJK` method will execute the target function  N times, leaving out one region in each case.
    The returns can then be used to compute the JK covariance. 

    The regions are in terms of coordinates on the sphere: RA and DEC.
    One gives a file generated from :func:`~suchyta_utils.jk.GenerateJKRegions` 
    (or equivalently from `kmeans_radec <https://github.com/esheldon/kmeans_radec>`_ itself),
    which specifies the N centers of k-means clusters. The data is then classified into these clusters,
    leaving out one cluster in each JK realization.

    Once :func:`~suchyta_utils.jk.SphericalJK.DoJK` executes, :func:`~suchyta_utils.jk.SphericalJK.GetResults` can be called to get the results.
    The results are whatever the return values are from the target function. There will be an N-dimensional array (N being the number of regions)
    of these return values from the target function. The returns from running the target function over the full area, without subsampling, are also available.

    The constructor arguments here are a bit tedious to explain in words, so I suggest you look at the examples. 
    The names are meant to be transparent, so once you understand what's going on you'll probably be able to just use those without much else.
    Nevertheless, the explanations below can be helpful too, so they're there.

    Parameters
    ----------
    target (func)
        The target function, which you want to perform the JKing with. 
    jkargs (strucutred or 2D array/array of structured or 2D arrays)
        Arguments to `target` which will change in each JK realization. These are either structured arrays or 2D arrays with RA/DEC columns.
        The arrays will be subsampled in each JK realization, based on the RA/DEC columns. 
        If you only need to subsample a single array in the JKing, `jkargs` can just be a single array.
        To subsample multiple quantites, put them in a list.
    jkargsby (array/array of arrays)
        The names/integers of the RA/DEC columns in each of the `jkargs`. Each is given as a 2D list, e.g. ['alphawin_j2000_i', 'deltawin_j2000_i'].
        For a single array in `jkargs` this single array can be given, otherwise put a set of them in a list.
        String column names are appropriate is that component of `jkargs` is a structured array.
        If it is a 2D array use integers, e.g. [0, 1].
    jkargspos (int/array of ints)
        Positions of `jkargs` in `target`. This can be a single number of a single `jkargs`, otherwise it's an array.
        For example, if you have a function ``Test(arr1, arr2, other1, other2, other3='bar')``,
        where you need to subsample ``arr1`` and ``arr2`` in JK realizations, you could use ``jkargs=[arr1,arr2], jkargspos=[0,1]`` or ``jkargs=[arr2,arr1], jkargspos=[1,0]``.
        The default assumes that ``jkargs`` are the first N arguments of your function, in the order you gave them.
        So in these examples, the default in either case would be: ``jkargspos=[0,1]``.
    nojkargs (list)
        This is a list, much like usual python *args to a function, for additional arguments which `target` takes, but should not be subsampled in the JKing.
        The default is an empty list, meaning no additional arguments.
    nojkargspos (int/array of ints)
        Analogous to `jkargspos`, but for the positions of `nojkargs` instead of `jkargs`.
        The intergers are for the positions in the full `target` argument list, not just amound `nojkargs`.
        For example for ``Test(arr1, other1, arr2, other2, other3='bar')``, where again you want to subsample ``arr1`` and ``arr2``
        could look like ``jkargs=[arr1,arr2], jkargspos=[0,2], nojkargs=[other1,other2], nojkargspos=[1,3]``.
        This would be the default for `nojkargspos`; the default "fills in" `nojkargs` 
        into the function arguments (after accounting for `jkargs`/`jkargspos`), in the order you gave them.
    nojkkwargs (dict)
        Keyword arguments to give to `target`. Behavior is like usual python **kwargs. These will not be subsampled.
    regions (str)
        Filename of the file with the k-means centers. Can be a file generated by :func:`~suchyta_utils.jk.GenerateJKRegions`.
    mpi (bool)
        Whether or not to use MPI multiprocessing. You must have mpi4py installed for this to work. See the examples for more about using MPI.

    
    Returns
    -------
    SphericalJK object

    """
    
    def __init__(self, target=None, jkargs=None, jkargsby=None, jkargspos=None, nojkargs=[], nojkargspos=[], nojkkwargs={}, regions=None, mpi=False):
        self.SetTarget(target=target, jkargs=jkargs, jkargsby=jkargsby, jkargspos=jkargspos, nojkargs=nojkargs, nojkkwargs=nojkkwargs)
        self.SetMPI(mpi=mpi)
        if regions is not None:
            self.SetRegions(regions)
    
    def SetMPI(self, mpi=True):
        """
        Set whether or not MPI will be used by default. :class:`Constructor <suchyta_utils.jk.SphericalJK>` calls this internally.

        Parameters
        ----------
        mpi (bool)
            
        Returns
        -------
        None

        """

        self.usempi = mpi

    def SetRegions(self, file):
        """
        Set the k-means cluster centers. :class:`Constructor <suchyta_utils.jk.SphericalJK>` calls this internally.

        Parameters
        ----------
        file (str)
            Filename of the file with the centers

        Returns
        -------
        None

        """

        self.regions = file

    def GetResults(self, jk=True, full=True):
        """
        Get the results of the JKing. To be called after :func:`~DoJK`

        Parameters
        ----------
        jk (bool)
            Whether or not to return the JK results
        full (bool)
            Whether or not to return the results over the full area

        Returns
        -------
        results (dict):
            Dictonary of the return values. ``if jk`` there is a ``"jk"`` key with an `njack` length array of the JK returns,
            as well as an ``"it"`` key which has an array with the integer iteration indexes.
            ``if full`` the dictionary has a ``"full"`` key with the return values over the full area.

        """

        ret = {}
        if jk:
            ret['jk'] = self.jkresults
            ret['it'] = self.jkindex
        if full:
            ret['full'] = self.fullresults
        return ret


    def DoJK(self, regions=None, mpi=None):
        """
        Do the JK iterations. An iteration with the full data will also be done. 
        After calling this, use the :func:`~GetResults` to return results.

        Parameters
        ----------
        regions (str)
            Filename of the file with the k-means centers. Can be a file generated by :func:`~suchyta_utils.jk.GenerateJKRegions`. 
            If given, will override `regions` given in :class:`constructor <suchyta_utils.jk.SphericalJK>` or :func:`~SetRegions`.
        mpi (bool)
            Whether or not to use MPI multiprocessing. You must have mpi4py installed for this to work. See the examples for more about using MPI.
            If given, Will override `mpi` setting givin in :class:`constructor <suchyta_utils.jk.SphericalJK>` or :func:`~SetMPI`.
           
        Returns
        -------
        None

        """

        if mpi is None:
            mpi = self.usempi

        if mpi:
            self.usempi = True
            from mpi4py import MPI
            import suchyta_utils.mpi as _mpi
            rank = MPI.COMM_WORLD.Get_rank()

        self._FindIndex(regions=regions)
        if mpi:
            index = None
            if rank==0:
                index = self._GenerateJKs()
            index = _mpi.Scatter(index)
        else:
            index = self._GenerateJKs()

        
        results = []
        for i in range(len(index)):
            args = self._GetArgs(index[i])
            res = self.target(*args, **self.nojkkwargs)
            results.append(res)
        results = _np.array(results)

        if mpi:
            results, index = _mpi.Gather(results, index)
            if rank==0:
                self._GetResults(results, index)
        else:
            self._GetResults(results, index)


    def _GetResults(self, results, index):
        full = (index==-1)
        self.jkresults = results[-full]
        self.jkindex = index[-full]
        self.fullresults = results[full][0]


    def _GenerateJKs(self):
        index = _np.arange(self.njack)
        index = _np.append(index, -1)
        return index


    def _GetArgs(self, jack):
       
        args = [None] * (len(self.jkargspos) + len(self.nojkargs))
        for i in range(len(self.jkargspos)):

            if jack > -1:
                cond = (self.index[i]!=jack)
            else:
                cond = _np.array( [True]*len(self.index[i]) )

            args[self.jkargspos[i]] = self.jkargs[i][cond]

        for i in range(len(self.nojkargspos)):
            args[self.nojkargspos[i]] = self.nojkargs[i]
        return args


    def _FindIndex(self, regions=None):
        if regions is None:
            try:
                regions = self.regions
            except:
                raise Exception('You must specify a regions file to use')

        self.index = []
        if type(regions)==str:
            regions = [regions] * len(self.jkargs)
        if len(regions)!=len(self.jkargs):
            raise Exception('Number or regions files (%i) does not match the number of jkargs (%i)' %(len(regions),len(self.jkargs)))

        for i in range(len(self.jkargs)):

            centers = _np.loadtxt(regions[i])
            self.njack = centers.shape[0]
            km = kmeans_radec.KMeans(centers)

            ra, dec = self._GetRaDec(i)
            rdi = _np.zeros( (len(ra),2) )
            rdi[:,0] = ra
            rdi[:,1] = dec
            index = km.find_nearest(rdi)
            self.index.append(index)


    def _GetRaDec(self, pos):
        if self.jkarrtype[pos]=='str':
            ra = self.jkargs[pos][self.jkargsby[pos][0]]
            dec = self.jkargs[pos][self.jkargsby[pos][1]]
        else:
            ra = self.jkargs[pos][:, self.jkargsby[pos][:,0]]
            dec = self.jkargs[pos][:, self.jkargsby[pos][:,1]]
        return ra, dec


    def SetTarget(self, target=None, jkargs=None, jkargsby=None, jkargspos=None, nojkargs=[], nojkargspos=[], nojkkwargs={}):
        """
        Set all the stuff that has to do with the target function. This is called internally by the :class:`constructor <suchyta_utils.jk.SphericalJK>`. 
        Refer there for help.

        """

        if target is None:
            raise Exception('Must set a target function')
        if jkargs is None:
            raise Exception('Must give some data to jackknife')
        if jkargsby is None:
            raise Exception('Must specify which part of the data are the coordinates to jackknife with')
        
        self.target = target
        self._FormatJKArgs(jkargs, jkargsby)
        if jkargspos is None:
            self.jkargspos = _np.arange(len(self.jkargs))
        if len(self.jkargs)!=len(self.jkargspos):
            raise Exception('Mismatch in length between jkargs and jkargspos')

        self.nojkkwargs = nojkkwargs
        self.nojkargs = self._EnsureArray(nojkargs)
        self.nojkargspos = self._EnsureArray(nojkargspos)

        if self.nojkargspos==[]:
            full = len(self.jkargs) + len(self.nojkargs)
            fullpos = _np.arange(full)
            found = _np.in1d(fullpos, self.jkargspos)
            self.nojkargspos = fullpos[-found]

        if len(self.nojkargs)!=len(self.nojkargspos):
            raise Exception('Mismatch in length between nojkargs and nojkargspos')


    def _EnsureArray(self, arr):
        if type(arr)==str:
            arr = [arr]

        try: 
            size = len(arr)
        except:
            arr = [arr]
        return arr


    def _FormatJKArgs(self, jkargs, jkargsby):
        try:
            _jkargsby = _np.array(jkargsby)
        except:
            raise Exception('jkargsby must be some kind of array: either shape (2,) or (N,2)')

        if _jkargsby.ndim==1:
            _jkargsby = _np.array( [_jkargsby] )
       
        self.jkarrtype = []
        for i in range(len(_jkargsby)):
            if len(_jkargsby[i])!=2:
                raise Exception('Each element of jkargsby must have shape (2,)')

            a = str(_jkargsby[i][0])
            try:
                a = int(a)
            except:
                pass
            b = str(_jkargsby[i][1])
            try:
                b = int(b)
            except:
                pass
            
            if type(a)!=type(b):
                raise Exception('Mismatching jargsby types, position: %i'%(i))

            if type(a)==int:
                self.jkarrtype.append('int')
            else:
                self.jkarrtype.append('str')

        self.jkargsby = _jkargsby
        if len(self.jkargsby)==1:
            if len(jkargs)!=1:
                jkargs = [jkargs]
        for i in range(len(self.jkargsby)):
            if self.jkarrtype[i]=='str':
                try:
                    jkargs[i][self.jkargsby[i][0]]
                except:
                    raise Exception('Column %s in jkargs position %i does not exist' %(self.jkargsby[i][0],i))
                try:
                    jkargs[i][self.jkargsby[i][1]]
                except:
                    raise Exception('Column %s in jkargs position %i does not exist' %(self.jkargsby[i][1],i))
            else:
                try:
                    jkargs[i][:, self.jkargsby[i][0]]
                except:
                    raise Exception('Position %i in jkargs position %i does not exist' %(self.jkargsby[i][0],i))
                try:
                    jkargs[i][:, self.jkargsby[i][1]]
                except:
                    raise Exception('Position %i in jkargs position %i does not exist' %(self.jkargsby[i][1],i))
        self.jkargs = jkargs



def GenerateJKRegions(ra, dec, njack, jfile, maxiter=200, tol=1.0e-5):
    """
    Generate k-means clusters from a set of data, using `kmeans_radec <https://github.com/esheldon/kmeans_radec>`_.
    If you're unfamilar with the k-means algorithm, the `wikipedia page <https://en.wikipedia.org/wiki/K-means_clustering>`_ is helpful.
    For roughly uniform data, it generates N-clusters of roughly equal cardinality.
    Here, distances are computed on the surface of the unit sphere, and coordinates are given as RA/DEC.

    Parameters
    ----------
    ra (float array)
        Right ascension values for each data point.
    dec (float array)
        Declination values for each data point.
    njack (int)
        Number of k-means clusters to generate, i.e. the number of JK regions.
    jfile (str)
        Output file name to save the regions.
    maxiter (int)
        Maximum number of iterations for the k-means generation, see `kmeans_radec documentation <https://github.com/esheldon/kmeans_radec>`_.
    tol (float)
        Tolerance level needed to be considered converged, see `kmeans_radec documentation <https://github.com/esheldon/kmeans_radec>`_.

    Returns
    -------
    None

    """

    rd = _np.zeros( (len(ra),2) )
    rd[:,0] = ra
    rd[:,1] = dec

    km = kmeans_radec.kmeans_sample(rd, njack, maxiter=maxiter, tol=tol)
    
    if not km.converged:
        raise RuntimeError("k means did not converge")

    dir = os.path.dirname(jfile)
    if not os.path.exists(dir):
        os.makedirs(dir)

    _np.savetxt(jfile, km.centers)

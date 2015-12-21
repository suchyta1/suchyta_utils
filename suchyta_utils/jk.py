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
    
    def __init__(self, target=None, jkargs=None, jkargsby=None, jkargspos=None, nojkargs=[], nojkargspos=[], nojkkwargs={}, mpi=False, regions=None):
        self.SetTarget(target=target, jkargs=jkargs, jkargsby=jkargsby, jkargspos=jkargspos, nojkargs=nojkargs, nojkkwargs=nojkkwargs)
        self.SetMPI(mpi=mpi)
        if regions is not None:
            self.SetRegions(regions)
    
    def SetMPI(self, mpi=True):
        self.usempi = mpi

    def SetRegions(self, file):
        self.regions = file

    def GetResults(self, jk=True, full=True):
        ret = {}
        if jk:
            ret['jk'] = self.jkresults
            ret['it'] = self.jkindex
        if full:
            ret['full'] = self.fullresults
        return ret


    def DoJK(self, mpi=None, regions=None):
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
       
        args = [None] * (len(self.jkargspos)+ len(self.nojkargs))
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

        centers = _np.loadtxt(regions)
        self.njack = centers.shape[0]
        km = kmeans_radec.KMeans(centers)

        self.index = []
        for i in range(len(self.jkargs)):
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

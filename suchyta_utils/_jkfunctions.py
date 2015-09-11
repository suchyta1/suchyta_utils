#!/usr/bin/env python

import kmeans_radec
import numpy as np
import os



def EnforceArray2D(arr):
    arr = EnforceArray(arr, arr)
    arr = EnforceArray(arr, arr[0])
    return arr

def EnforceArray1D(arr):
    arr = EnforceArray(arr, arr)
    return arr

def EnforceArray(arr, check):
    try:
        len(check)
    except:
        arr = [arr]
    return arr



"""
For the target function you give to JacknifeOnSphere, the first argument MUST be an array of the JKed tables. 
Even if you're only JKing a single data table, you'll receive this as an array like [JKed table].
I've done it like this, so it's easy to JK multiple tables simultaneously if you want.
You can define extra arguments if you want too.

You MUST return a 2D array from your target function. 
The first return element is an array of the quantities you want to find the covariances of.
The second return element is whatever else you want, maybe coordinates, or something else.
If you don't want anything, just set it to None. You should return it in an array though, this makes it easier to keep track of things.
"""

def Dummy(arr):
    a = arr[0]['mag_i']
    bins = np.arange(18, 25, 0.2)
    hist, bins = np.histogram(a, bins=bins)
    return [ [hist], [len(hist)] ]

def Dummy2(arr, thing, binsize=0.1):
    sim = arr[0][thing]
    des = arr[1][thing]
    bins = np.arange(18, 25, binsize)
    cent = (bins[1:] + bins[0:-1]) / 2

    shist, bins = np.histogram(sim, bins=bins)
    dhist, bins = np.histogram(des, bins=bins)

    ss = float(len(sim))
    dd = float(len(des))

    s = shist / ss
    d = dhist / dd

    return [ [s,d], [cent, ss, dd] ]


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

Returns a 4-element array:
    Element 0: Array of the answers over the full area of the things you JKed
    Element 1: Covariances of the things you JKed
    Element 2: Array of whatever extra stuff you returned from your target function during the call for the full area
    Element 3: Like Element 2, but with arrays in each component for the extra returns of the target function in each JK realization.


"""

def JackknifeOnSphere(jarrs, jras, jdecs, jfunc, jargs=[], jkwargs={}, jtype='generate', jfile=None, njack=24, generateonly=False, gindex=0, varonly=False, save=None):
    jarrs = EnforceArray2D(jarrs)
    jras = EnforceArray2D(jras)
    jdec = EnforceArray2D(jdecs)

    if jtype=='generate':
        rdi = np.zeros( (len(jarrs[gindex]),2) )
        rdi[:,0] = jarrs[gindex][jras[gindex]]
        rdi[:,1] = jarrs[gindex][jdecs[gindex]]

        if jfile is None:
            jfile = 'JK-{0}.txt'.format(njack)
        km = kmeans_radec.kmeans_sample(rdi, njack, maxiter=200, tol=1.0e-5)
        if not km.converged:
            raise RuntimeError("k means did not converge")
        np.savetxt(jfile, km.centers)
        if generateonly:
            return jfile

    elif jtype=='read':
        centers = np.loadtxt(jfile)
        km = kmeans_radec.KMeans(centers)
        njack = len(centers)
    
    ind = []
    for i in range(len(jarrs)):
        rdi = np.zeros( (len(jarrs[i]),2) )
        rdi[:,0] = jarrs[i][jras[i]]
        rdi[:,1] = jarrs[i][jdecs[i]]
        index = km.find_nearest(rdi)
        ind.append(index)
    
    full_j, full_other = jfunc(jarrs, *jargs, **jkwargs)
    full_j = EnforceArray2D(full_j)
    full_other = EnforceArray1D(full_other)

    it_j = []
    it_other = []
    frac = []
    for j in range(njack):
        print 'JK %i' %(j)

        ja = []
        f = []

        for i in range(len(full_other)):
            if j==0:
                it_other.append( [] )
        
        for i in range(len(full_j)):
            if j==0:
                it_j.append( [] )

        for i in range(len(jarrs)):
            if j==0:
                #it_j.append( [] )
                frac.append( [] )

            cut = (ind[i]==j)
            ja.append(jarrs[i][-cut])

            ff = np.sum(-cut)/float(len(cut))
            f.append(ff)

        i_j, i_other = jfunc(ja, *jargs, **jkwargs)
        i_j = EnforceArray2D(i_j)
        i_other = EnforceArray1D(i_other)

        for i in range(len(i_j)):
            it_j[i].append( np.copy(i_j[i]) )
        for i in range(len(jarrs)):
            frac[i].append( f[i] )
        for i in range(len(i_other)):
            it_other[i].append( np.copy(i_other[i]) )

    for i in range(len(it_j)):
        it_j[i] = np.array(it_j[i])

    for i in range(len(frac)):
        frac[i] = np.array(frac[i])


    cov_j = []
    for k in range(len(full_j)):

        if varonly:
            cov = np.power( np.std(it_j[k], axis=0), 2.0 ) * njack * float(njack-1)/njack
            cov_j.append(cov)

        else:
            csize = len(full_j[k])
            cov = np.zeros( (csize,csize) )
            
            for i in range(csize):
                for j in range(i, csize):
                    cov[i,j] =  np.sum( (it_j[k][:,i] - full_j[k][i]) * (it_j[k][:,j] - full_j[k][j]) ) * float(njack-1)/njack
                    #cov[i,j] =  np.sum( (it_j[k][:,i] - full_j[k][i]) * (it_j[k][:,j] - full_j[k][j])  * frac[k] )

                    if i!=j:
                        cov[j,i] = cov[i,j]
            cov_j.append(cov)

    if save is not None:
        
        vec = os.path.join(save, 'vec')
        cov = os.path.join(save, 'cov')
        other = os.path.join(save, 'other')
        
        Write2Dir(vec, full_j)
        Write2Dir(cov, cov_j)
        Write2Dir(other, full_other)

    return [full_j, cov_j, full_other, it_other]



def ArrFunc(func, *arr):
    if len(arr)==1:
        return func(arr[0])
    else:
        return [ func(arg) for arg in arr ]


def GenerateRegions(jarrs, jras, jdecs, jfile, njack, gindex, jtype):

    if jtype=='generate':
        rdi = np.zeros( (len(jarrs[gindex]),2) )
        rdi[:,0] = jarrs[gindex][jras[gindex]]
        rdi[:,1] = jarrs[gindex][jdecs[gindex]]

        if jfile is None:
            jfile = 'JK-{0}.txt'.format(njack)
        km = kmeans_radec.kmeans_sample(rdi, njack, maxiter=200, tol=1.0e-5)

        if not km.converged:
            raise RuntimeError("k means did not converge")
        np.savetxt(jfile, km.centers)

    elif jtype=='read':
        centers = np.loadtxt(jfile)
        km = kmeans_radec.KMeans(centers)
        njack = len(centers)

    return [km, jfile]


def AssignIndex(jarrs, jras, jdecs, km):
    ind = []
    for i in range(len(jarrs)):
        rdi = np.zeros( (len(jarrs[i]),2) )
        rdi[:,0] = jarrs[i][jras[i]]
        rdi[:,1] = jarrs[i][jdecs[i]]
        index = km.find_nearest(rdi)
        ind.append(index)
    return ind

def FullResults(jfunc, jargs, jkwargs, jarrs):
    full_j, full_other = jfunc(jarrs, *jargs, **jkwargs)
    full_j = EnforceArray2D(full_j)
    full_other = EnforceArray1D(full_other)
    return [full_j, full_other]


def Iterate(jfunc, jargs, jkwargs, jarrs, ind, njack, full_j, full_other, rtype, itsave):
    it_j = []
    it_other = []
    frac = []
    for j in range(njack):
        print 'JK %i' %(j)

        ja = []
        f = []

        for i in range(len(full_other)):
            if j==0:
                it_other.append( [] )
        
        for i in range(len(full_j)):
            if j==0:
                it_j.append( [] )

        for i in range(len(jarrs)):
            if j==0:
                #it_j.append( [] )
                frac.append( [] )
            
            if rtype=='n-1':
                cut = (ind[i]==j)
            elif rtype=='/n':
                cut = (ind[i]!=j)

            ja.append(jarrs[i][-cut])

            ff = np.sum(-cut)/float(len(cut))
            f.append(ff)

        i_j, i_other = jfunc(ja, *jargs, **jkwargs)
        i_j = EnforceArray2D(i_j)
        i_other = EnforceArray1D(i_other)

        for i in range(len(i_j)):
            it_j[i].append( np.copy(i_j[i]) )
        for i in range(len(jarrs)):
            frac[i].append( f[i] )
        for i in range(len(i_other)):
            it_other[i].append( np.copy(i_other[i]) )

    for i in range(len(it_j)):
        it_j[i] = np.array(it_j[i])

    for i in range(len(frac)):
        frac[i] = np.array(frac[i])

    return it_j, it_other, frac


def GetCov(full_j, it_j, njack, varonly, rtype):
    cov_j = []

    cfac = 1.0/njack
    if rtype=='n-1':
        norm = cfac * (njack-1)
    elif rtype=='/n':   
        norm = cfac / (njack-1)

    for k in range(len(full_j)):

        if varonly:
            cov = np.power( np.std(it_j[k], axis=0), 2.0 ) * njack * norm
            cov_j.append(cov)

        else:
            csize = len(full_j[k])
            cov = np.zeros( (csize,csize) )
            
            for i in range(csize):
                for j in range(i, csize):
                    cov[i,j] =  np.sum( (it_j[k][:,i] - full_j[k][i]) * (it_j[k][:,j] - full_j[k][j]) ) * norm
                    #cov[i,j] =  np.sum( (it_j[k][:,i] - full_j[k][i]) * (it_j[k][:,j] - full_j[k][j])  * frac[k] )

                    if i!=j:
                        cov[j,i] = cov[i,j]
            cov_j.append(cov)
    return cov_j

def SetupSave(save, itsave):
    r1 = None
    r2 = None

    if save is not None:
        vec = os.path.join(save, 'vec')
        cov = os.path.join(save, 'cov')
        other = os.path.join(save, 'other')
        r1 = [vec, cov, other]

        if itsave:
            vec_it = os.path.join(vec, 'it')
            other_it = os.path.join(other, 'it')
            r2 = [vec_it, other_it]

    return [r1, r2] 


def FullSave(save, full):
    for i in range(len(save)):
        Write2Dir(save[i], full[i])

def ItSave(itsave, it):
    for i in range(len(it)):
        for j in range(len(it[i])):
            dir = os.path.join(itsave[i], '%i'%j)
            Write2Dir(dir, it[i][j])



def Write2Dir(dir, arr):
    if not os.path.exists(dir):
        os.makedirs(dir)
    for i in range(len(arr)):
        file = os.path.join(dir, '%i.txt'%(i))
        np.savetxt(file, arr[i])

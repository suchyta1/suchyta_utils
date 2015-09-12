<h1> Contents </h1>

`suchyta_utils` is an assortment of some Python tools written by Eric Suchyta in the course of working in DES.
There are "submodules" for different kinds of utilities:

* `mpi` -- Functions for wrapping and simplifying syntax for common MPI functionality to distribute/collect data during MPI jobs
* `jk` -- Automatically generate (and save) jackknife realizations of an arbitrary function given data distributed on a sphere.
* `db` -- Talking to your area of the DES DB (or that of other users)


<h1> Examples </h1>

Refer to the [examples](https://github.com/suchyta1/suchyta_utils/tree/master/examples) directory for usage help.


<h1> Installation </h1>

Installation follows the usual GitHub --> Python pattern.
```
git clone https://github.com/suchyta1/suchyta_utils.git
cd suchyta_utils/
python setup.py install
```

If you don't have root privileges you can use `--user` flag:
```
python setup.py install --user
```


<h1> Requirements </h1>

`suchyta_utils.mpi`
* [mpi4py](https://pypi.python.org/pypi/mpi4py)
* [numpy](http://www.numpy.org/)

`suchyta_utils.jk`
* [kmeans_radec](https://github.com/esheldon/kmeans_radec)
* [numpy](http://www.numpy.org/)


`suchtya_utils.db`
*[desdb](https://github.com/esheldon/desdb)
*[cx_Oracle](http://cx-oracle.sourceforge.net/)

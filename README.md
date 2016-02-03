<h1> Documentation </h1>

Refer to the [Sphinx documenation](http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/doc/html/) for the python module.
Apologies if not everything is documented, I'm trying to keep it up to date.
The shell commands are things that could easily just be aliases in rc file,
but the ones included here are shell agnostic. They work in bash or tcsh.


<h1> Examples </h1>

Refer to the [examples](https://github.com/suchyta1/suchyta_utils/tree/master/examples) directory for some example scripts.


<h1> Installation </h1>

Installation follows the usual GitHub --> Python pattern.
```
git clone https://github.com/suchyta1/suchyta_utils.git
cd suchyta_utils/
python setup.py install
```

If you do not have root privileges you can use `--user` flag:
```
python setup.py install --user
```


<h1> Requirements </h1>

* [numpy](http://www.numpy.org/)
* [cx_Oracle](http://cx-oracle.sourceforge.net/)
* [esutil](https://github.com/esheldon/esutil)
* [desdb](https://github.com/esheldon/desdb)
* [kmeans_radec](https://github.com/esheldon/kmeans_radec)
* [mpi4py](https://pypi.python.org/pypi/mpi4py) (Package compiles without it, but will not be able to use the mpi submodule without it)




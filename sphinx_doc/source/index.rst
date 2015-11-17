.. Test documentation master file, created by
   sphinx-quickstart on Mon Nov  9 08:43:05 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``suchyta_utils`` Documentation
================================

This package is a set of python utilies written by Eric Suchyta while working with DES data and Balrog.
You shouldn't consider this an "offical" product of anything, but I do my best keep it bug-free and you're welcome to use it if you'd like
(be it for DES, astronomy, or whatever).
The code is hosted on `github <https://github.com/suchyta1/suchyta_utils>`_.

If you find a bug or have a pull request, feel free to post an issue on github.
You can also email me: eric.d.suchyta@gmail.com.

Throughout the documentation, I assume you're reasonably fluent in python, numpy, matplotlib, FITS files, etc.
I'm also lazy about some argument types: lists, tuples, numpy arrays, etc. I just call arrays.
Stuctured arrays can be actual numpy objects, or more or less equivalent things from pyfits or similar.

The package is organized into submodules, whose contents are limited to a farily specific scope.

.. toctree::
   :maxdepth: 2

   hp -- HEALPix <hp.rst>
   jk -- Jackknife <jk.rst>
   db -- Database <db.rst>
   slr -- Stellar locus regression <slr.rst>
   plot -- Matplotlib utils <plot.rst> 
   system <system.rst>
   balrog <balrog.rst>


Installation
============

Installation consists of cloning the code from `github <https://github.com/suchyta1/suchyta_utils>`_, then doing the usual ``setup.py`` incantations::
       
        git clone https://github.com/suchyta1/suchyta_utils.git
        cd suchyta_utils/
        python setup.py install

If you don't have root privileges you can use ``--user`` flag (or ``--prefix=/whatever/path/you/want``)::
        
        python setup.py install --user
        


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


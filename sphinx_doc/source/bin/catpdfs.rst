catpdfs
=======

``catpdfs`` concatenates multiple PDF files into a single file. This is implemented with ``gs``.

Usage
-----
::

$ catpdfs <inputs> <output>

* inputs -- (Space seperated) input files to concatenate. You can use globs too.
* output -- Output file name.

Example
-------
::

$ catpdfs ./input1.pdf subdir/*.pdf output.pdf

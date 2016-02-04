screen-retach
=============

``screen-retach`` reattaches to a named GNU screen session.
Used with sessions started by `screen-split <http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/doc/html/bin/screen-split.html>`_,
and detached with `screen-detach <http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/doc/html/bin/screen-detach.html>`_,
this will maintain any splitting to the appearace, which screen doesn't usually do by default.
Screen sessions persist even if you log out of your ssh session.

Usage
-----
::

$ screen-retach <name>

* name -- The name of the session.

Example
-------
::

$ screen-retach test

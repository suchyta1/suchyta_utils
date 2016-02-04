screen-detach
=============

``screen-detach`` detaches from a named GNU screen session. If started by `screen-split <http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/doc/html/bin/screen-split.html>`_,
and reattached using `screen-retach <http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/doc/html/bin/screen-retach.html>`_,
this maintains any splitting to the appearace, which screen doesn't usually do by default.
Screen sessions persist even if you log out of your ssh session.

Usage
-----
::

$ screen-detach <name>

* name -- The name of the session.

Example
-------
::

$ screen-detach test

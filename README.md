KMRNext (tentative name)
========================

BACKEND RUNTIME
---------------

This software has two backends to manage data and tasks.

* serial
* kmr

"serial" backend uses only one node and "kmr" backend uses multiple nodes to
run in parallel.

INSTALL
-------

Building this software with the "serial" backend requires a C++ compiler that
supports C++98 standard.  Doing with the "kmr" backend requires an MPI
library that supports MPI 2.2. and KMR (http://mt.aics.riken.jp/kmr).

This software can be build by a usual building flow, configure, make and
make install. To switch the backend, specify options to the configure script.
The default backend is "serial".

* configure with the "serial" backend

      $ ./configure --prefix=PATH_TO_INSTALL_DIR
      $ make
      $ make install

* configure with the "kmr" backend

      $ ./configure --prefix=PATH_TO_INSTALL_DIR \
        --with-backend=kmr --with-kmr=PATH_TO_KMR_DIR
      $ make
      $ make install

To enable DEBUG mode, specify "--enable-debug" option to the configure script.

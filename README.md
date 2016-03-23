scanalyzer
==========

An interactive visualizer for interferometric visibilities.

![Screenshot](http://newton.cx/~peter/files/scanalyzer-screenshot.png)


Installation
------------

The scanalyzer can be installed with a standard invocation of

```
python setup.py install
```

It depends on a few packages that you may not already have installed:

* The Gtk+3 graphical toolkit, accessible through the
  [PyGObject](https://wiki.gnome.org/Projects/PyGObject) layer. If you use
  [Anaconda Python](https://docs.continuum.io/anaconda/index), you can install
  the needed libraries with the command:

  ```
  conda install -c pkgw gtk3 pygobject3
  ```

* The `pwkit` support library, installable using PIP or Conda:

  ```
  conda install -c pkgw pwkit
  ```

* If you want to read MIRIAD data sets, you need the `miriad-python` package.
  **TODO**: this will be the pain point in the installation until I make some
  tidier packages for Anaconda. You can learn about downloading and installing
  `miriad-python`
  [here](https://www.cfa.harvard.edu/~pwilliam/miriad-python/). It requires an
  autotools-based CARMA MIRIAD install, as provided by
  [this tarball](https://www.cfa.harvard.edu/~pwilliam/miriad-macport/miriad-latest.tar.gz).

* If you want to read CASA datasets, you need the CASA libraries and their
  Python bindings. If you use Anaconda Python, you can install the needed
  libraries with the command:

  ```
  conda install -c pkgw casa-data casa-python
  ```

**TODO**: no reason we can't use AIPY to read in MIRIAD data.
[This function](https://github.com/HERA-Team/scanalyzer/blob/master/scanalyzer/transpose.py#L125)
is the one that needs changing.


Usage
-----

Launch with

```
scanalyzer go <visibility-data-set> [output-flag-file-name]
```

**TODO**: flags can be saved to an output file on the fly, but the file cannot
currently be reloaded at a later invocation!

Either MIRIAD or CASA visibility data sets are supported if the necessary
supporting modules are available. There are a few other modes that can be
learned about by simply running `scanalyzer` with no arguments.

Super-quick summary of the interactive interface:

* `q` and `w` to move between baselines
* `f` and `l` to go to the **f**irst and **l**ast basepols in the current
  selection
* `j` and `k` or left-button and right-button single clicks to move between
  looking at the real, imaginary, amplitude, and phase components.
* `a` to toggle showing of **a**ll or only unflagged data
* Hold down the `Ctrl` key to get a live readout of data properties
* `s` to flag (**s**uppress) the current basepol at all frequencies and times
* `Ctrl-w` or `Ctrl-q` to quit.

Draw boxes by clicking and dragging, then:

* `z` to **z**oom in to that box; `Z` (capital zee) to zoom out
* `c` to flag the box only on the **c**urrent baseline
* `b` to flag it on all **b**aselines
* `h` to extend the box infinitely in the **h**orizontal (frequency) direction
* `v` to extend the box infinitely in the **v**ertical (time) direction

Use `e` to toggle the box-**e**diting mode where you can see existing boxes,
apply or unapply them to the current baseline, and edit their shapes.

Many other undocumented features exist.


Authors
-------

Peter K. G. Williams and collaborators.


License
-------

The [MIT license](http://opensource.org/licenses/MIT). See the file LICENSE.

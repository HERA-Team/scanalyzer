scanalyzer
==========

An interactive visualizer for interferometric visibilities.


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

* The `miriad-python` package for reading MIRIAD datasets. **TODO**: this will
  be the pain point in the installation until I make some tidier packages for
  Anaconda. You can learn about downloading and installing `miriad-python`
  [here](https://www.cfa.harvard.edu/~pwilliam/miriad-python/). It requires an
  autotools-based CARMA MIRIAD install, as provided by
  [this tarball](https://www.cfa.harvard.edu/~pwilliam/miriad-macport/miriad-latest.tar.gz).

**TODO**: no reason we can't use AIPY to read in the data.
[This function](https://github.com/HERA-Team/scanalyzer/blob/master/scanalyzer/transpose.py#L125)
is the one that needs changing.


Usage
-----

Super-quick summary:

* `q` and `w` to move between baselines
* `j` and `k` or left-button and right-button single clicks to move between looking at the
  real, imaginary, amplitude, and phase components.
* `a` to toggle showing of flagged data

To flag interactively, click and drag to draw boxes, then:

* `c` to flag the box only on the current baseline
* `b` to flag it on all baselines
* `e` to toggle the box-editing mode where you can see existing boxes, apply
  or unapply them to the current baseline, and edit their shapes.


Authors
-------

Peter K. G. Williams and collaborators.


License
-------

The [MIT license](http://opensource.org/licenses/MIT). See the file LICENSE.

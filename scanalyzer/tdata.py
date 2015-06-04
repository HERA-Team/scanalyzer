# -*- mode: python; coding: utf-8 -*-
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""pwkit.scanalyzer.tdata - interfaces for accessing transposed UV data

In some ways, it might be more correct to call this gridded UV data, but we're
gridding in the freq/time plane, and the terminology might be confusing with
the sort of U/V gridding that happens during imaging.

This module requires Gtk+3 for the interactive transform controller widgets.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = (b'''GridAxis VisGrid TransposeData DataTransform
           MeanSubtractTransform TopocentrizeTransform
           AverageTransform SubsetTransform HacktasticDDRTransform
           CustomTransform''').split ()

import numpy as np
from gi.repository import Gtk


class GridAxis (object):
    __slots__ = ['start', 'step', 'size']

    def __init__ (self, start, step, size):
        self.start = start
        self.step = step
        self.size = size


    def __getstate__ (self):
        return self.start, self.step, self.size


    def __setstate__ (self, s):
        self.start, self.step, self.size = s


    def clone (self):
        return self.__class__ (self.start, self.step, self.size)


    def centers (self):
        return np.arange (self.size) * self.step + self.start


    def valsToScalars (self, vals):
        '''A "scalar" here is a mapping onto the grid such that 0 corresponds
        to the left-hand-side of the first grid box and 1 corresponds to the
        right-hand-side of the last one. That is, a scalar value of zero
        does *not* correspond to the center of the first grid box.'''
        vals = np.asarray (vals)
        return (vals - self.start + 0.5 * self.step) / (self.step * self.size)


    def valsToIndices (self, vals):
        '''An "index" here is a mapping onto the grid such that 0 corresponds
        to the center of the first grid box, 1 corresponds to the center of
        the second box, etc. These indices are NOT necessarily integers!'''
        vals = np.asarray (vals)
        return (vals - self.start) / self.step


    def scalarsToVals (self, s):
        s = np.asarray (s)
        return (s * self.size - 0.5) * self.step + self.start


    def indicesToVals (self, idxs):
        '''This maps index values into the grid (as if it were an array) onto
        coordinate values. Equivalent to centers()[idxs] for integer idxs.'''

        idxs = np.asarray (idxs)
        return idxs * self.step + self.start


    def getBounds (self):
        return (self.start - 0.5 * self.step,
                self.start + self.step * (self.size - 0.5))


    def overlaps (self, vstart, vend):
        if vend < vstart:
            vstart, vend = vend, vstart
        if vend < (self.start - 0.5 * self.step):
            return False
        if vstart > (self.start + self.step * (self.size - 0.5)):
            return False
        return True


    def overlapSlice (self, vstart, vend):
        """Returns a slice object representing which indices along the
        axis would fall within the specified (inclusive) value range,
        considering the centers of each bin. Returns None if no
        overlap."""

        if vend < vstart:
            vstart, vend = vend, vstart
        if vend < (self.start - 0.5 * self.step):
            return None
        if vstart > (self.start + self.step * (self.size - 0.5)):
            return None

        sidx = int (np.ceil ((vstart - self.start) / self.step))
        eidx = 1 + int (np.floor ((vend - self.start) / self.step))
        return slice (sidx, eidx)


class VisGrid (object):
    taxis = faxis = None
    data = flags = uvw = weights = None
    _where = None

    def alloc (self, taxis, faxis):
        nt, nf = taxis.size, faxis.size

        self.taxis = taxis
        self.faxis = faxis
        self.data = np.empty ((nt, nf), dtype=np.complex64)
        self.flags = np.empty ((nt, nf), dtype=np.bool_)
        self.uvw = np.empty ((nt, 3), dtype=np.double)
        self.weights = np.empty (nt, dtype=np.double)

        return self


    def valid (self):
        return True


    def invalidate (self):
        pass


    def where (self):
        if self._where is None:
            self._where = np.nonzero (self.flags)
        return self._where


    def allFlagged (self):
        return self.where ()[0].size == 0


    def ffdata (self):
        # flattened, flagged data
        return self.data[self.where ()]


    def clone (self):
        other = VisGrid ()

        if self.taxis is not None:
            other.taxis = self.taxis.clone ()
            other.faxis = self.faxis.clone ()
            other.data = self.data.copy ()
            other.flags = self.flags.copy ()
            other.uvw = self.uvw.copy ()
            other.weights = self.weights.copy ()

            if self._where is not None:
                other._where = tuple (x.copy () for x in self._where)

        return other


    def shallowcopy (self):
        other = VisGrid ()
        other.taxis = self.taxis
        other.faxis = self.faxis
        other.data = self.data
        other.flags = self.flags
        other.uvw = self.uvw
        other.weights = self.weights
        return other


    # Coordinate transforms -- the "x" coordinate is frequency
    # and the "y" coordinate is time, because that's how I want
    # to render things.

    def physical_to_scalar (self, time, ghz, clamp=False):
        sx = self.faxis.valsToScalars (ghz)
        sy = self.taxis.valsToScalars (time)

        if clamp:
            np.clip (sx, 0, 1, sx)
            np.clip (sy, 0, 1, sy)

        return sx, sy


    def scalar_to_physical (self, sx, sy):
        return self.taxis.scalarsToVals (sy), \
            self.faxis.scalarsToVals (sx)


    def effMaxCoord (self):
        w = self.where ()
        if w[0].size == 0:
            return None

        a = np.abs (self.data[w])
        amax = a.argmax ()
        mt = self.taxis.indicesToVals (w[0][amax])
        mf = self.faxis.indicesToVals (w[1][amax])
        return mt, mf


"""
TransposeData objects are just collections of VisGrids. More
sophisticated functionality can depend on various resources that may
or may not be available, so we have a feature-testing system for these
things to be determined on-the-fly. The features are named and
specified below.

delays
  The object has a getDelay(bp) method that takes a 2-tuple of
  antpols as an argument and returns the instrumental delay associated
  with that particular baseline in nanoseconds.

meanuvw
  The object has a getMeanUVW(bp) method that takes a 2-tuple of
  antpols as an argument and returns a 3-element double ndarray containing
  the mean UVW of the specified baseline within the considered data.
  The UVW coordinates are measured in nanoseconds.

  This function should be fast, because it can be used to sort the list
  of basepols in the scanalyzer. In these situations it will be called
  sequentially for every basepol.

controllerwidget
  The object has a getControlWidget(invalidate) method that returns
  a GTK+ widget with a user interface for controlling parameters that
  related to the data produced by the data source. *invalidate* is a
  callable, taking no arguments and returning None, that should be
  called when those parameters change in such a way as to affect the
  data.

checkallflagged
  The object has a method knownAllFlagged(bp) that takes a 2-tuple of
  antpols as an argument and returns True if all of the data for that
  basepol are currently flagged. As with the meanuvw feature, this function
  should be fast since it's used in filtering the list of basepols.
"""

class TransposeData (object):
    def getGrid (self, bp):
        raise NotImplementedError ()

    def getBPs (self):
        raise NotImplementedError ()

    def hasFeature (self, name):
        return False


class DataTransform (TransposeData):
    parent = None

    def setParent (self, parent):
        self.parent = parent
        return self

    def getGrid (self, bp):
        # Likely to be overridden ...
        return self.parent.getGrid (bp)

    def getBPs (self):
        return self.parent.getBPs ()

    # We explicitly list the features that we support, since they may
    # imply semantics that implementations must be aware of.

    _myfeatures = ()
    _passthroughfeatures = frozenset (('delays', 'meanuvw', 'checkallflagged'))

    def hasFeature (self, name):
        if name in self._myfeatures:
            return True
        if name not in self._passthroughfeatures:
            return False
        return self.parent.hasFeature (name)

    def getDelay (self, bp):
        return self.parent.getDelay (bp)

    def getMeanUVW (self, bp):
        return self.parent.getMeanUVW (bp)

    def knownAllFlagged (self, bp):
        return self.parent.knownAllFlagged (bp)


class MeanSubtractTransform (DataTransform):
    _myfeatures = frozenset (('controllerwidget', ))
    _removephase = True

    def getGrid (self, bp):
        vgrid = self.parent.getGrid (bp).shallowcopy ()
        vgrid.data = vgrid.data.copy ()
        m = vgrid.ffdata ().mean ()

        if self._removephase:
            vgrid.data *= m.conj () / np.abs (m)
            vgrid.data -= np.abs (m)
        else:
            vgrid.data -= m

        return vgrid

    def getControlWidget (self, invalidate):
        vb = Gtk.VBox ()

        l = Gtk.Label ()
        l.set_markup ('<b>Subtract mean</b>')
        l.set_alignment (0.0, 0.5)
        vb.pack_start (l, False, True, 2)

        cb = Gtk.CheckButton ('Remove mean phase')
        cb.set_active (self._removephase)
        vb.pack_start (cb, False, True, 2)

        def ontoggle (button):
            self._removephase = button.get_active ()
            invalidate ()

        cb.connect ('toggled', ontoggle)
        return vb


class TopocentrizeTransform (DataTransform):
    _myfeatures = frozenset (('controllerwidget', ))

    def setParent (self, parent):
        if not parent.hasFeature ('delays'):
            raise ValueError ('need delays to topocentrize')
        self.parent = parent

    def getGrid (self, bp):
        vgrid = self.parent.getGrid (bp).shallowcopy ()
        vgrid.data = vgrid.data.copy ()
        vgrid.uvw = vgrid.uvw.copy ()

        # Correction phases. The sign of 'w' is such that multiplying
        # the phased-up data by 2*pi*i*w takes them into the
        # topocentric frame. (As opposed to -2*pi*i*w.) There is also
        # the fixed delay back to the correlator. Defining the
        # relative antenna delay analogously to MIRIAD's UV coordinate
        # definition (coord = b2 - b1; tau = tau2 - tau1), the delay
        # correction can be un-done by multiplying by -2*pi*i*nu*tau.
        # I think.

        fixed = self.parent.getDelay (bp)
        weff = vgrid.uvw[:,2] - fixed
        ws = np.outer (weff, vgrid.faxis.centers ())
        vgrid.data *= np.exp (2 * np.pi * (0.+1j) * ws)
        # Not sure that this is right, but seems to make sense:
        vgrid.uvw[:,2] = 0
        return vgrid

    def getDelay (self, bp):
        return 0 # Also unsure about this one.

    def getControlWidget (self, invalidate):
        l = Gtk.Label ('Phase to topocenter')
        l.set_alignment (0.0, 0.5)
        return l


class AverageTransform (DataTransform):
    _myfeatures = frozenset (('controllerwidget', ))

    fbins = 2
    tbins = 2
    slop = 0.5

    def getGrid (self, bp):
        fbins = self.fbins
        tbins = self.tbins

        vgrid = self.parent.getGrid (bp)

        fstep = vgrid.faxis.step * fbins
        tstep = vgrid.taxis.step * tbins
        nfreq = int (vgrid.faxis.size / fbins)
        ntime = int (vgrid.taxis.size / tbins)
        sfreq = vgrid.faxis.start + vgrid.faxis.step * 0.5 * (fbins - 1)
        stime = vgrid.taxis.start + vgrid.taxis.step * 0.5 * (tbins - 1)

        faxis = GridAxis (sfreq, fstep, nfreq)
        taxis = GridAxis (stime, tstep, ntime)

        newgrid = VisGrid ().alloc (taxis, faxis)

        # I feel like there must be a better way to do this
        # box-averaging, somewhere out there. But I'm not even quite
        # sure what to call this operation. It's quite slow as
        # implemented.

        slarea = fbins * tbins * self.slop

        oldd = vgrid.data
        newd = newgrid.data
        oldf = vgrid.flags
        newf = newgrid.flags
        olduvw = vgrid.uvw
        newuvw = newgrid.uvw
        oldweights = vgrid.weights
        newweights = newgrid.weights

        uvwwork = np.empty ((3, ))

        for ntidx in xrange (ntime):
            stidx = ntidx * tbins
            etidx = stidx + tbins

            for nfidx in xrange (nfreq):
                sfidx = nfidx * fbins
                efidx = sfidx + fbins

                f = oldf[stidx:etidx,sfidx:efidx]
                d = oldd[stidx:etidx,sfidx:efidx]
                w = np.where (f)
                n = w[0].size

                if n < slarea:
                    # not enough subcells filled in to make the
                    # averaged cell valid.
                    newd[ntidx,nfidx] = 0
                    newf[ntidx,nfidx] = 0
                else:
                    newd[ntidx,nfidx] = d[w].mean ()
                    newf[ntidx,nfidx] = 1

            avgweight = 0
            uvwwork.fill (0)
            wtwork = 0

            for otidx in xrange (stidx, stidx + tbins):
                n = np.where (oldf[otidx])[0].size
                uvwwork += olduvw[otidx] * n
                avgweight += n
                wtwork += oldweights[otidx] * n

            if avgweight == 0:
                newuvw[ntidx].fill (0)
                newweights[ntidx].fill (0)
            else:
                newuvw[ntidx] = uvwwork / avgweight
                newweights[ntidx] = wtwork / avgweight

        return newgrid


    def getControlWidget (self, invalidate):
        t = Gtk.Table (4, 2)

        l = Gtk.Label ()
        l.set_markup (r'<b>Average</b>')
        l.set_alignment (0.0, 0.5)
        t.attach (l, 0, 2, 0, 1, Gtk.FILL, Gtk.FILL, 2, 2)

        l = Gtk.Label ('F binsize:')
        l.set_alignment (1.0, 0.5)
        t.attach (l, 0, 1, 1, 2, Gtk.FILL, Gtk.FILL, 2, 2)

        l = Gtk.Label ('T binsize:')
        l.set_alignment (1.0, 0.5)
        t.attach (l, 0, 1, 2, 3, Gtk.FILL, Gtk.FILL, 2, 2)

        l = Gtk.Label ('Slop:')
        l.set_alignment (1.0, 0.5)
        t.attach (l, 0, 1, 3, 4, Gtk.FILL, Gtk.FILL, 2, 2)

        s = Gtk.SpinButton ()
        # FIXME derive limits from actual grid (here and tbins)
        s.get_adjustment ().set_all (self.fbins, 1, 10000, 1, 8, 0)
        def changed (adj):
            self.fbins = int (adj.get_value ())
            invalidate ()
        s.get_adjustment ().connect ('value-changed', changed)
        t.attach (s, 1, 2, 1, 2, Gtk.EXPAND|Gtk.FILL, Gtk.FILL, 2, 2)

        s = Gtk.SpinButton ()
        s.get_adjustment ().set_all (self.tbins, 1, 10000, 1, 8, 0)
        def changed (adj):
            self.tbins = int (adj.get_value ())
            invalidate ()
        s.get_adjustment ().connect ('value-changed', changed)
        t.attach (s, 1, 2, 2, 3, Gtk.EXPAND|Gtk.FILL, Gtk.FILL, 2, 2)

        s = Gtk.SpinButton (digits=2)
        s.get_adjustment ().set_all (self.slop, 0, 1, 0.1, 0.3, 0)
        def changed (adj):
            self.slop = adj.get_value ()
            invalidate ()
        s.get_adjustment ().connect ('value-changed', changed)
        t.attach (s, 1, 2, 3, 4, Gtk.EXPAND|Gtk.FILL, Gtk.FILL, 2, 2)

        return t


class SubsetTransform (DataTransform):
    _bounds = None

    def setBounds (self, stime, etime, sfreq, efreq):
        self._bounds = (stime, etime, sfreq, efreq)

    # FIXME-ish: we just passthrough getMeanUVW. To be strictly
    # correct, we should recompute for the subset, but that's
    # slow, and getMeanUVW needs to be fast for sorting.

    def getGrid (self, bp):
        vgrid = self.parent.getGrid (bp)
        stime, etime, sfreq, efreq = self._bounds

        # Map the physical units of the bounds onto the grid,
        # rounding to integral grid elements. Have to be a bit careful
        # since the ideal coordinate system for doing the math ("scalars")
        # is not the same as how we index the data array.

        ntime, nfreq = vgrid.data.shape

        sts, ets = vgrid.taxis.valsToScalars ([stime, etime])
        sfs, efs = vgrid.faxis.valsToScalars ([sfreq, efreq])

        stidx = int (round (sts * ntime))
        etidx = int (round (ets * ntime))
        sfidx = int (round (sfs * nfreq))
        efidx = int (round (efs * nfreq))

        # Create and use new grid.

        nfreq = efidx - sfidx
        ntime = etidx - stidx

        if nfreq < 1 or ntime < 1:
            return vgrid

        fstep = vgrid.faxis.step
        tstep = vgrid.taxis.step

        newgrid = VisGrid ()
        # "stime" has different semantics for our bounds and for GridAxis.
        newgrid.taxis = GridAxis (stime + tstep * 0.5, tstep, ntime)
        newgrid.faxis = GridAxis (sfreq + fstep * 0.5, fstep, nfreq)
        newgrid.data = vgrid.data[stidx:etidx,sfidx:efidx]
        newgrid.flags = vgrid.flags[stidx:etidx,sfidx:efidx]
        newgrid.uvw = vgrid.uvw[stidx:etidx]
        newgrid.weights = vgrid.weights[stidx:etidx]

        return newgrid


class HacktasticDDRTransform (DataTransform):
    _myfeatures = frozenset (('controllerwidget', ))

    def getGrid (self, bp):
        vgrid = self.parent.getGrid (bp)
        omax = vgrid.ffdata ().max ()

        vgrid = vgrid.shallowcopy ()
        vgrid.data = vgrid.data.copy ()
        vgrid.flags = vgrid.flags.copy ()

        import numpy.fft
        d = vgrid.data
        d[np.where (np.logical_not (vgrid.flags))] = 0.0
        d = vgrid.data = numpy.fft.fft2 (d)
        vgrid.flags.fill (1)
        # tons of zero-frequency power makes viz hard
        vgrid.flags[0,0] = 0
        # todo: deconvolve
        # more natural ordering
        nt, nf = d.shape
        tmp = d[:,:(nf+1)//2].copy ()
        d[:,:(nf+1)//2] = d[:,nf//2:]
        d[:,nf//2:] = tmp
        tmp = d[:(nt+1)//2,:].copy ()
        d[:(nt+1)//2,:] = d[nt//2:,:]
        d[nt//2:,:] = tmp
        # compress scale hand-tunedly
        nmax = np.abs (vgrid.ffdata ()).max ()
        mag = np.abs (vgrid.data)
        scale = np.maximum (1, 70 * mag / nmax)
        vgrid.data /= scale
        vgrid.data *= 70 * omax / nmax

        return vgrid

    def getControlWidget (self, invalidate):
        l = Gtk.Label ('Hacktastic DDR')
        l.set_alignment (0.0, 0.5)
        return l


class CustomTransform (DataTransform):
    _myfeatures = frozenset (('controllerwidget', ))
    _gridfunc = None

    def _passthrough (self, bp):
        return self.parent.getGrid (bp)


    def _getfunc (self):
        import sys
        sys.path = ['.'] + sys.path

        try:
            import transform
            transform = reload (transform)
            self._gridfunc = transform.getGrid
        except ImportError:
            print ('Custom transform fail: no "transform.py" in current directory')
            self._gridfunc = self._passthrough
        finally:
            del sys.path[0]


    def getGrid (self, bp):
        if self._gridfunc is None:
            self._getfunc ()
        return self._gridfunc (self.parent, bp)


    def getControlWidget (self, invalidate):
        b = Gtk.Button ('Reload Custom')

        def clicked (*args):
            self._gridfunc = None
            invalidate ()
            return True

        b.connect ('clicked', clicked)
        return b


# TODO:
# Stokes-processing transform
# subtract mean along freq axis (or time axis)
# apply external gain/bandpass cal dataset
# apply additional (un-editable?) flags

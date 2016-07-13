# -*- mode: python; coding: utf-8 -*-
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""scanalyzer.ui - implementation of the UI.

This uses the Gtk+3 graphical toolkit.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = b'Scanalyzer'.split ()

from gi.repository import Gdk, Gtk
import cairo
import numpy as np
from . import mtutil as util
from .mtutil import POL_XX, POL_YY, POL_XY, POL_YX, POL_RR, POL_LL

def jdToFull (t):
    return '%.5f' % t


_phuncert_c1 = np.pi / np.sqrt (3)
_phuncert_c2 = np.sqrt (9. / 2 / np.pi**3)
_phuncert_poly = [2.73387475e-04,
                  -3.78537034e-03,
                  3.69150376e-03,
                  1.70698123e-01,
                  -9.86931168e-01,
                  1.92600603e+00]

def phaseUncert (snr):
    """\
Compute the uncertainty in a phase measurement, given a signal to
noise ratio for the visibility amplitude. Based on Equations 9.38,
9.48, and 9.53 in Thompson, Moran, & Swenson. I evaluated the full
probability distribution function numerically and fit a 5th-order
polynomial to it to get the expression to use for the intermediate-SNR
case, and chose the domain cutoffs by eye.

There are discontinuities at the domain cutoffs. For my current
purposes, they do not matter at all.
"""
    work = np.asarray (snr)

    # Not thrilled with evaluating all three cases for all values
    # but I don't know of a natural way to do the switch statement
    # in an array mode, and these should be fast operations

    low = _phuncert_c1 * (1 - _phuncert_c2 * work)
    high = 1. / np.maximum (work, 1e-2) # avoid div by zero
    mid = np.polyval (_phuncert_poly, work)

    ishigh = (work > 5)
    islow = (work < 0.7)

    work = np.where (ishigh, high, mid)
    return np.where (islow, low, work)



# Simple classes for "rendering" arrays of arbitary data (possibly
# complex-valued) into RGB values. Could use better terminology for
# this process. The destination is an array of int32s of the same
# shape as the input that should be treated as an RGB buffer.

class SimpleRenderer (object):
    def __call__ (self, data, destimg):
        raise NotImplementedError ()


class SingleRampRenderer (SimpleRenderer):
    '''Renders a single color ramp: minval is black, maxval is
    full-brightness basecolor (specifically, 255 * basecolor).'''

    def __init__ (self, minval, maxval, basecolor):
        self.minval = minval
        self.maxval = maxval
        self.basecolor = basecolor


    def __call__ (self, data, destimg):
        v = ((data - self.minval) * 255 / (self.maxval - self.minval))

        # If we multiply a floating-point value by basecolor, the
        # fractional parts pollute different color bits, making the
        # effect weird. So we have to round off to integers before
        # the multiplication for color. Same applies in the functions
        # below.

        np.clip (v, 0, 255, destimg)
        np.multiply (destimg, self.basecolor, destimg)


class DualRampRenderer (SimpleRenderer):
    '''Renders a dual-color ramp: minval is green, maxval is red,
    in between is blue.'''

    def __init__ (self, minval, maxval):
        self.minval = minval
        self.maxval = maxval


    def __call__ (self, d, destimg):
        d = ((d - self.minval) / (self.maxval - self.minval)).clip (0, 1)

        destimg[:] = (2*d - 1).clip (0, 1) * 255
        np.left_shift (destimg, 8, destimg)
        np.add ((1 - (2*d).clip (0, 1)) * 255, destimg, destimg, casting='unsafe')
        np.left_shift (destimg, 8, destimg)
        np.add ((0.5 - np.abs (d - 0.5)) * 255, destimg, destimg, casting='unsafe')


class WrappingRenderer (SimpleRenderer):
    '''Renders a wrapping value. Input is assumed to be between -period/2 and
    +period/2.'''

    def __init__ (self, period):
        self.period = period


    def __call__ (self, d, destimg):
        p2 = 0.5 * self.period
        tp2 = p2 / 3
        ttp2 = 2 * p2 / 3

        destimg[:] = (ttp2 - np.abs (d.clip (-p2, tp2) + tp2)) * 255 / ttp2
        np.left_shift (destimg, 8, destimg)
        np.add ((ttp2 - np.abs (d.clip (-tp2, p2) - tp2)) * 255 / ttp2,
               destimg, destimg, casting='unsafe')
        np.left_shift (destimg, 8, destimg)
        np.add ((ttp2 - np.abs (((d % (2 * p2)) - p2).clip (-ttp2, ttp2))) * 255 / ttp2,
               destimg, destimg, casting='unsafe')


# VGrid renderers take a vgrid as an argument and have the option of
# using more structured data

def simpleSampleData (vgrid):
    if vgrid.allFlagged ():
        return np.asarray (0.)
    return vgrid.ffdata ()


def wrapSimpleRenderer (preproc, simplerender):
    def rendervgrid (tdata, bp, vgrid, destimg, showall):
        if showall:
            simplerender (preproc (vgrid.data), destimg)
        else:
            w = vgrid.where ()
            temp = np.zeros (w[0].size, dtype=np.int32)
            simplerender (preproc (vgrid.data[w]), temp)
            destimg[w] = temp

    return rendervgrid


def makeAmpRenderer (samplevgrid):
    a = np.abs (simpleSampleData (samplevgrid))
    srr = SingleRampRenderer (a.min (), a.max (), 0x000100)
    return wrapSimpleRenderer (np.abs, srr)


def makeRealRenderer (samplevgrid):
    from operator import attrgetter
    r = simpleSampleData (samplevgrid).real
    drr = DualRampRenderer (r.min (), r.max ())
    return wrapSimpleRenderer (attrgetter ('real'), drr)


def makeImagRenderer (samplevgrid):
    from operator import attrgetter
    i = simpleSampleData (samplevgrid).imag
    drr = DualRampRenderer (i.min (), i.max ())
    return wrapSimpleRenderer (attrgetter ('imag'), drr)


def makePhaseRenderer (samplevgrid=None):
    p = lambda x: np.arctan2 (x.imag, x.real)
    return wrapSimpleRenderer (p, WrappingRenderer (2 * np.pi))


def makePhaseUncertRenderer (samplevgrid=None):
    if samplevgrid.allFlagged ():
        umax = 0
    else:
        signal = np.abs (samplevgrid.ffdata ())
        weights2d = np.outer (samplevgrid.weights, np.ones (samplevgrid.faxis.size))
        w = samplevgrid.where ()
        invnoise = np.sqrt (weights2d[w])
        temp = np.zeros (w[0].size, dtype=np.int32)
        umax = phaseUncert (signal * invnoise).max ()

    srr = SingleRampRenderer (0, umax, 0x10000)

    def rendervgrid (tdata, bp, vgrid, destimg, showall):
        if showall:
            # We need some contortions to use array broadcasting
            snr = (np.abs (vgrid.data).T * np.sqrt (vgrid.weights)).T
            srr (phaseUncert (snr), destimg)
        else:
            signal = np.abs (vgrid.ffdata ())
            weights2d = np.outer (vgrid.weights, np.ones (vgrid.faxis.size))
            w = vgrid.where ()
            invnoise = np.sqrt (weights2d[w])
            temp = np.zeros (w[0].size, dtype=np.int32)
            srr (phaseUncert (signal * invnoise), temp)
            destimg[w] = temp

    return rendervgrid


def makeSystemPhaseRenderer (samplevgrid=None):
    wr = WrappingRenderer (1)

    def rendervgrid (tdata, bp, vgrid, destimg, showall):
        # This code was derived from tdata.TopocentrizeTransform
        fixed = tdata.getDelay (bp)
        weff = vgrid.uvw[:,2] - fixed
        phases = np.outer (weff, vgrid.faxis.centers ())
        phases = ((phases + 0.5) % 1.) - 0.5
        wr (phases, destimg)

    return rendervgrid


# The main scanalyzer window

from .tdata import *

SEL_FILT_UNSUP = 0
SEL_FILT_SUP = 1
SEL_FILT_ALL = 2

SEL_BL_BOTH = 3 # bitfield interpretation
SEL_BL_CROSS = 1
SEL_BL_AUTO = 2

SEL_POL_ALL = 10
SEL_POL_INTENS = 11
SEL_POL_NONINTENS = 12
SEL_POL_XX = POL_XX
SEL_POL_YY = POL_YY
SEL_POL_RR = POL_RR
SEL_POL_LL = POL_LL

SEL_ANTPOL_ALL = 0
SEL_ANTPOL_ANT = 1
SEL_ANTPOL_ANTPOL = 2

SORT_NUMERICAL_BL_POL = 0
SORT_NUMERICAL_ANTPOL = 1
SORT_RANDOM = 2
SORT_UVDIST = 3

QUANT_AMP = 0
QUANT_PHA = 1
QUANT_REAL = 2
QUANT_IMAG = 3
QUANT_PHUNCERT = 4
QUANT_SYSPHASE = 5

RESIZE_LEFT = 0
RESIZE_RIGHT = 1
RESIZE_TOP = 2
RESIZE_BOTTOM = 3
RESIZE_NONE = 4

INVALIDATE_FLAGS = 0
INVALIDATE_VGRID = 1
INVALIDATE_RENDER = 2
INVALIDATE_SCALE = 3

SELFLAG_HSPAN = 1 << 0
SELFLAG_VSPAN = 1 << 1

_resizeCursors = { RESIZE_TOP: 'top_side', RESIZE_RIGHT: 'right_side',
                   RESIZE_BOTTOM: 'bottom_side', RESIZE_LEFT: 'left_side' }

_rendererMakers = {QUANT_AMP: makeAmpRenderer,
                   QUANT_PHA: makePhaseRenderer,
                   QUANT_REAL: makeRealRenderer,
                   QUANT_IMAG: makeImagRenderer,
                   QUANT_PHUNCERT: makePhaseUncertRenderer,
                   QUANT_SYSPHASE: makeSystemPhaseRenderer,
                   }

from os.path import dirname, join

class Scanalyzer (object):
    def __init__ (self, data, flagapi):
        if flagapi is not None:
            from .flag import FlagTransform
            data = FlagTransform (flagapi).setParent (data)

        self.flagapi = flagapi

        self.tdata_base = data
        self.vgrid_base = None
        self.tdata = None
        self.vgrid = None

        builder = self.builder = Gtk.Builder ()
        from pkg_resources import resource_filename
        builder.add_from_file (resource_filename ('scanalyzer', 'scanalyzer.ui'))

        events = {'on_sa_darea_draw': self._draw,
                  'on_sa_darea_button_press_event': self._button_press,
                  'on_sa_darea_button_release_event': self._button_release,
                  'on_sa_darea_motion_notify_event': self._motion_notify,
                  'on_darea_key_press_event': self._on_darea_key_press,
                  'on_darea_key_release_event': self._on_darea_key_release,
                  'on_darea_leave_notify_event': self._on_darea_leave_notify,
                  'on_help_btn_clicked': self._help_clicked,
                  'on_help_deleted': self._help_deleted,
                  'on_scanalyzer_key_press_event': self._on_keypress,
                  'on_scanalyzer_delete_event': self._on_deleted,
                }
        builder.connect_signals (events)

        self.win = builder.get_object ('scanalyzer')
        #self.win.set_title ('Scanalyzer - specific info')
        self.darea = builder.get_object ('sa_darea')

        self.helpwin = builder.get_object ('help_win')

        self._bpsel_init ()
        self._xform_init ()
        self._zoom_init ()
        self._shownquant_init ()
        self._darea_init ()
        self._selection_init ()
        self._tooltip_init ()


    def _help_clicked (self, widget):
        self.helpwin.show_all ()

    def _help_deleted (self, widget, event):
        widget.hide ()
        return True

    # ############################################################
    # Cleanup ...

    def close (self):
        self.win.destroy ()

        self.tdata_base = None
        self.vgrid_base = None
        self.vgrid = None
        self.win = None
        self.darea = None
        self.builder = None


    def _on_deleted (self, window, event):
        self.close ()


    # ############################################################
    # Invalidation of data structures for redraw

    def invalidate (self, itype):
        if itype <= INVALIDATE_FLAGS:
            if self.inboxsel:
                self._compute_boxsel_boxes ()
        if itype <= INVALIDATE_VGRID:
            self.vgrid = None
        if itype <= INVALIDATE_RENDER:
            # Effective flags still valid, but rendering mode has changed.
            # (Some changes to rendering mode do not require a call to
            # this function since they're caught automatically via
            # img_rendered_for.)
            self.img_rendered_for = None
            self.img_data = None
        if itype <= INVALIDATE_SCALE:
            # Rendering mode still valid, but scaling to window has changed.
            # (Again, some changes to scaling do not require a call to this
            # function since they're caught via img_scaled_for.)
            self.img_scaled_for = None
            self.img_scaledsurf = None

        self.darea.queue_draw ()


    # ############################################################
    # Basepol selection management: filtering, sorting

    def _bpsel_init (self):
        self.bpsel_index = 0
        self.bpsel_list = None

        seenaps = set ()

        # FIXME: keep listings updated with changes imposed by
        # transformations (e.g. Stokes processing)
        for bp in self.tdata_base.getBPs ():
            seenaps.add (bp[0])
            seenaps.add (bp[1])

        get = self.builder.get_object

        # Combobox of flaggedness-filtering choices. We get underlying
        # info from the transpose file, so we can check for unflagged
        # basepols even when flagapi is None.
        ls = Gtk.ListStore (int, str)
        ls.append ((SEL_FILT_UNSUP, 'Unflagged basepols'))
        ls.append ((SEL_FILT_SUP, 'Flagged basepols'))
        ls.append ((SEL_FILT_ALL, 'Both'))
        cb = get ('sel_sup_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.set_model (ls)
        cb.connect ('changed', self._bpsel_update)
        self._bpsel_sup_cb = cb

        # Auto/cross-correlation type filtering
        ls = Gtk.ListStore (int, str)
        ls.append ((SEL_BL_BOTH, 'Auto/cross'))
        ls.append ((SEL_BL_CROSS, 'Cross-corrs'))
        ls.append ((SEL_BL_AUTO, 'Auto-corrs'))
        cb = get ('sel_bltype_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.set_model (ls)
        cb.connect ('changed', self._bpsel_update)
        self._bpsel_bltype_cb = cb

        # Combobox of polarization-filtering choices
        ls = Gtk.ListStore (int, str)
        ls.append ((SEL_POL_ALL, 'All polarizations'))
        ls.append ((SEL_POL_INTENS, 'Intensity-type'))
        ls.append ((SEL_POL_NONINTENS, 'Non-intensity-type'))
        ls.append ((SEL_POL_XX, 'XX'))
        ls.append ((SEL_POL_YY, 'YY'))
        ls.append ((SEL_POL_RR, 'RR'))
        ls.append ((SEL_POL_LL, 'LL'))
        cb = get ('sel_pol_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.set_model (ls)
        cb.connect ('changed', self._bpsel_update)
        self._bpsel_pol_cb = cb

        # Ant component of antpol filtering
        ls = Gtk.ListStore (int, str)
        ls.append ((-1, '(All ants)'))
        for antnum in sorted (set (util.apAnt (ap) for ap in seenaps)):
            ls.append ((antnum, self.format_antnum (antnum)))
        cb = get ('sel_apfilter_ant_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.set_model (ls)
        cb.connect ('changed', self._bpsel_apfilter_ant_changed)
        self._bpsel_apfilter_ant_cb = cb

        # Pol component
        ls = Gtk.ListStore (int, str)
        ls.append ((-1, '(pols)'))
        for fpol in sorted (set (util.apFPol (ap) for ap in seenaps)):
            ls.append ((fpol, util.fPolNames[fpol]))
        cb = get ('sel_apfilter_pol_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.set_model (ls)
        cb.connect ('changed', self._bpsel_update)
        self._bpsel_apfilter_pol_cb = cb

        # Sorting
        ls = Gtk.ListStore (int, str)
        ls.append ((SORT_NUMERICAL_BL_POL, 'Numerical (BL/Pol)'))
        ls.append ((SORT_NUMERICAL_ANTPOL, 'Numerical (Antpol)'))
        ls.append ((SORT_RANDOM, 'Random'))
        ls.append ((SORT_UVDIST, 'UV Distance'))
        cb = get ('sort_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.set_model (ls)
        cb.connect ('changed', self._bpsel_update)
        self._bpsel_sort_cb = cb

        # Labels showing prev/current/next basepol
        self._bpsel_labels = [get ('info_bp1_l'), get ('info_bp2_l'),
                              get ('info_bp3_l')]

        # Set up everything ...
        self._bpsel_sup_cb.set_active (0)
        self._bpsel_bltype_cb.set_active (1)
        self._bpsel_pol_cb.set_active (1)
        self._bpsel_apfilter_ant_cb.set_active (0)
        self._bpsel_apfilter_pol_cb.set_active (0)
        self._bpsel_sort_cb.set_active (0)


    def _bpsel_apfilter_ant_changed (self, widget):
        iter = self._bpsel_apfilter_ant_cb.get_active_iter ()
        if iter is None:
            antsel = -1
        else:
            antsel = self._bpsel_apfilter_ant_cb.get_model ().get_value (iter, 0)

        sensitive = antsel != -1
        self._bpsel_apfilter_pol_cb.set_sensitive (sensitive)
        self._bpsel_update ()


    def _bpsel_get_selected (self):
        # Load up filtering state

        iter = self._bpsel_sup_cb.get_active_iter ()
        if iter is None:
            supsel = SEL_FILT_ALL
        else:
            supsel = self._bpsel_sup_cb.get_model ().get_value (iter, 0)

        iter = self._bpsel_bltype_cb.get_active_iter ()
        if iter is None:
            bltsel = SEL_BL_BOTH
        else:
            bltsel = self._bpsel_bltype_cb.get_model ().get_value (iter, 0)

        iter = self._bpsel_pol_cb.get_active_iter ()
        if iter is None:
            polsel = SEL_POL_ALL
        else:
            polsel = self._bpsel_pol_cb.get_model ().get_value (iter, 0)

        iter = self._bpsel_apfilter_ant_cb.get_active_iter ()
        if iter is None:
            antsel = -1
        else:
            antsel = self._bpsel_apfilter_ant_cb.get_model ().get_value (iter, 0)

        if antsel == -1:
            antpolsel = SEL_ANTPOL_ALL
        else:
            iter = self._bpsel_apfilter_pol_cb.get_active_iter ()
            if iter is None:
                appolsel = -1
            else:
                appolsel = self._bpsel_apfilter_pol_cb.get_model ().get_value (iter, 0)

            if appolsel == -1:
                antpolsel = SEL_ANTPOL_ANT
            else:
                antpolsel = SEL_ANTPOL_ANTPOL

        # Set up filtering functions

        if supsel == SEL_FILT_UNSUP:
            supfunc = lambda bp: not self.tdata_base.knownAllFlagged (bp)
        elif supsel == SEL_FILT_SUP:
            supfunc = lambda bp: self.tdata_base.knownAllFlagged (bp)
        elif supsel == SEL_FILT_ALL:
            supfunc = lambda bp: True
        else:
            assert False

        if bltsel == SEL_BL_BOTH:
            bltfunc = lambda bp: True
        elif bltsel == SEL_BL_CROSS:
            def bltfunc (bp):
                ant1, ant2, pol = util.bp2aap (bp)
                return ant1 != ant2
        elif bltsel == SEL_BL_AUTO:
            def bltfunc (bp):
                ant1, ant2, pol = util.bp2aap (bp)
                return ant1 == ant2
        else:
            assert False

        if polsel == SEL_POL_ALL:
            polfunc = lambda bp: True
        elif polsel == SEL_POL_INTENS:
            polfunc = lambda bp: util.bpIsInten (bp)
        elif polsel == SEL_POL_NONINTENS:
            polfunc = lambda bp: not util.bpIsInten (bp)
        else:
            polfunc = lambda bp: (util.bp2aap (bp)[2] == polsel)

        if antpolsel == SEL_ANTPOL_ALL:
            apfunc = lambda bp: True
        elif antpolsel == SEL_ANTPOL_ANT:
            apfunc = lambda bp: (antsel in (util.bp2aap (bp)[0:2]))
        elif antpolsel == SEL_ANTPOL_ANTPOL:
            theantpol = util.antpol2ap (antsel, appolsel)
            apfunc = lambda bp: (theantpol in bp)

        # Run the filter!

        for bp in self.tdata_base.getBPs ():
            if supfunc (bp) and bltfunc (bp) and polfunc (bp) and apfunc (bp):
                yield bp


    def _bpsel_get_sorted (self):
        iter = self._bpsel_sort_cb.get_active_iter ()
        if iter is None:
            sorttype = SORT_NUMERICAL_BL_POL
        else:
            sorttype = self._bpsel_sort_cb.get_model ().get_value (iter, 0)

        if sorttype == SORT_NUMERICAL_BL_POL:
            def key (bp):
                bl, pol = util.bp2blpol (bp)
                return bl * 9 + abs (pol)
        elif sorttype == SORT_NUMERICAL_ANTPOL:
            key = lambda bp: bp
        elif sorttype == SORT_RANDOM:
            from random import random
            key = lambda bp: random ()
        elif sorttype == SORT_UVDIST:
            def key (bp):
                uvw = self.tdata.getMeanUVW (bp)
                return uvw[0]**2 + uvw[1]**2
        else:
            assert False

        return sorted (self._bpsel_get_selected (), key=key)


    def _bpsel_update (self, *args):
        # *args so we can be used as a widget change callback
        oldlist = self.bpsel_list
        oldindex = self.bpsel_index

        if oldlist is not None and len (oldlist) > 0:
            oldbp = oldlist[oldindex]
        else:
            oldbp = None

        newlist = self._bpsel_get_sorted ()

        if len (newlist) == 0:
            newbp = None
            newindex = 0
        else:
            if oldbp is None:
                newindex = 0
            else:
                # Idea: if the current basepol isn't in the new list,
                # we could check for other basepols with the same baseline
                # but different polarizations.
                checkindex = oldindex
                while checkindex < len (oldlist):
                    try:
                        newindex = newlist.index (oldlist[checkindex])
                        break
                    except ValueError:
                        pass
                    checkindex += 1
                else:
                    newindex = 0

            newbp = newlist[newindex]

        self.bpsel_list = newlist
        self.bpsel_index = newindex

        self._bpsel_labels_update ()

        if newbp != oldbp:
            self.invalidate (INVALIDATE_VGRID)


    def format_antnum (self, antnum):
        "Note that this operates on 1-based antenna names."
        return self.tdata_base.getAntName (antnum)


    def format_ap (self, ap):
        antidx = ap >> 3
        fp = ap & 0x7
        s = self.tdata_base.getAntName (antidx + 1)
        return s + util.fPolNames[fp]


    def format_bp (self, bp):
        return '-'.join (self.format_ap (ap) for ap in bp)


    def _bpsel_labels_stringify (self, index):
        if index < 0 or index >= len (self.bpsel_list):
            return '(none)'

        s = self.format_bp (self.bpsel_list[index])

        if index != self.bpsel_index:
            return s

        return '<b>' + s + '</b>'


    def _bpsel_labels_update (self):
        idx = self.bpsel_index
        self._bpsel_labels[0].set_markup (self._bpsel_labels_stringify (idx - 1))
        self._bpsel_labels[1].set_markup (self._bpsel_labels_stringify (idx))
        self._bpsel_labels[2].set_markup (self._bpsel_labels_stringify (idx + 1))


    def _bpsel_get_current (self):
        if len (self.bpsel_list) == 0:
            return None
        return self.bpsel_list[self.bpsel_index]


    def _bpsel_set_current (self, newindex):
        if self.bpsel_index == newindex:
            return

        self.bpsel_index = newindex
        self._bpsel_labels_update ()
        self.invalidate (INVALIDATE_VGRID)


    def _bpsel_select_next (self):
        if self.bpsel_index >= len (self.bpsel_list) - 1:
            return False
        self._bpsel_set_current (self.bpsel_index + 1)
        return True


    def _bpsel_select_prev (self):
        if self.bpsel_index == 0:
            return False
        self._bpsel_set_current (self.bpsel_index - 1)
        return True


    def _bpsel_select_next_or_different (self):
        if self.bpsel_index < len (self.bpsel_list) - 1:
            self._bpsel_set_current (self.bpsel_index + 1)
        elif self.bpsel_index > 0:
            self._bpsel_set_current (self.bpsel_index - 1)


    def _bpsel_select_first (self):
        self._bpsel_set_current (0)
        return True


    def _bpsel_select_last (self):
        if len (self.bpsel_list) == 0:
            return False
        self._bpsel_set_current (len (self.bpsel_list) - 1)
        return True


    # ############################################################
    # Transformation stack

    def _xform_init (self):
        get = self.builder.get_object

        ls = Gtk.ListStore (object, str)
        ls.append ((None, '(Add transform)'))
        ls.append ((MeanSubtractTransform, 'Subtract mean'))
        ls.append ((AverageTransform, 'Freq/Time Average'))
        ls.append ((TopocentrizeTransform, 'Topocentrize phases'))
        ls.append ((HacktasticDDRTransform, 'Hacktastic DDR'))
        ls.append ((CustomTransform, 'Custom'))

        cb = get ('xform_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.set_model (ls)
        cb.connect ('changed', self._xform_changed)

        self._xform_vb = get ('xform_vb')
        self._xform_list = []
        self._xform_tempskip = False

        cb.set_active (0)


    def _xform_changed (self, widget):
        iter = widget.get_active_iter ()
        if iter is None:
            return

        klass = widget.get_model ().get_value (iter, 0)
        if klass is None:
            return

        widget.set_active (0) # reset to "Add" item

        xform = klass ()
        assert xform.hasFeature ('controllerwidget')

        cbutton = Gtk.Button ('>')
        cbutton.connect ('popup-menu', lambda w: self._xform_popup_menu (xform, w, None))
        cbutton.connect ('clicked', lambda w: self._xform_popup_menu (xform, w, None))

        control = xform.getControlWidget (lambda: self.invalidate (INVALIDATE_VGRID))
        hb = Gtk.HBox (spacing=4)
        hb.pack_start (cbutton, False, False, 0)
        hb.pack_start (control, True, True, 0)
        hb.show_all ()

        self._xform_vb.pack_start (hb, False, True, 0)
        self._xform_list.append ((xform, hb, True))
        self.invalidate (INVALIDATE_VGRID)


    def _xform_popup_menu (self, xform, widget, event):
        menu = Gtk.Menu ()
        xl = self._xform_list

        for i, (ixform, hbox, enabled) in enumerate (xl):
            if ixform is xform:
                idx = i
                break
        else:
            print ('Cannot find clicked xform???')
            return False

        def moveup (item):
            if idx == 0:
                return False

            self._xform_vb.reorder_child (xl[idx][1], idx - 1)
            t = xl[idx]
            xl[idx] = xl[idx - 1]
            xl[idx - 1] = t
            self.invalidate (INVALIDATE_VGRID)
            return True

        def movedown (item):
            if idx == len (self._xform_list) - 1:
                return False

            self._xform_vb.reorder_child (xl[idx][1], idx + 1)
            t = xl[idx]
            xl[idx] = xl[idx + 1]
            xl[idx + 1] = t
            self.invalidate (INVALIDATE_VGRID)
            return True

        def delete (item):
            self._xform_vb.remove (xl[idx][1])
            del xl[idx]
            self.invalidate (INVALIDATE_VGRID)
            return True

        def togable (item):
            xform, hbox, enabled = xl[idx]
            enabled = not enabled
            xl[idx] = xform, hbox, enabled

            first = True
            for c in hbox.get_children ():
                if not first:
                    c.set_sensitive (enabled)
                first = False

            self.invalidate (INVALIDATE_VGRID)
            return True

        def additem (label, func, sensitive=True):
            item = Gtk.MenuItem (label)
            item.connect ('activate', func)
            item.set_sensitive (sensitive)
            menu.append (item)

        additem ('Move Up', moveup, i > 0)

        if enabled:
            additem ('Disable', togable)
        else:
            additem ('Enable', togable)

        additem ('Delete', delete)
        additem ('Move Down', movedown, i < len (xl) - 1)

        menu.show_all ()

        if event is None:
            button = 0
            time = Gtk.get_current_event_time ()
        else:
            button = event.button
            time = event.time

        menu.attach_to_widget (widget, None)
        menu.popup (None, None, None, button, time)
        return True


    def _xform_toggle_tempskip (self):
        self._xform_tempskip = not self._xform_tempskip
        sens = not self._xform_tempskip
        self._xform_vb.set_sensitive (sens)
        self.builder.get_object ('xform_cb').set_sensitive (sens)
        self.invalidate (INVALIDATE_VGRID)


    # ############################################################
    # Management of the visgrid as it applies to our UI

    def _vgrid_update (self, bp):
        if self.vgrid is not None:
            # To force an update, call self.invalidate (INVALIDATE_VGRID)
            return

        if bp is None:
            self.vgrid = self.vgrid_base = None
            return

        self.vgrid_base = self.tdata_base.getGrid (bp)
        self.effmaxcoord = self.vgrid_base.effMaxCoord ()

        # Set up the pipeline of transforms, with possible zoom at the end

        self.tdata = self.tdata_base

        if not self._xform_tempskip:
            for xform, hbox, enabled in self._xform_list:
                if enabled:
                    xform.setParent (self.tdata)
                    self.tdata = xform

        if self._zoom_bounds is not None:
            self._zoom_transform.setBounds (*self._zoom_bounds)
            self._zoom_transform.setParent (self.tdata)
            self.tdata = self._zoom_transform

        # Get transformed grid.

        self.vgrid = self.tdata.getGrid (bp)

        # Recompute boxsel boxes and selection since the parameters
        # may have changed if the transforms  changed the axes.
        # FIXME: boxsel should just stored physical coordinates,
        # not scalars.

        if self.boxsel_selected is not None:
            x1, x2, y1, y2, clipping, ident = self.boxsel_selected
            t1, t2, f1, f2 = self.flagapi.get (ident)[:4]

        if self.inboxsel:
            self._compute_boxsel_boxes ()

        if self.boxsel_selected is not None:
            x1, x2, y1, y2, clipping = self._boxsel_compute_boxinfo (t1, t2, f1, f2)

            if x1 > 1 or x2 < 0 or y1 > 1 or y2 < 0:
                # out of bounds
                self.boxsel_selected = None
            else:
                self.boxsel_selected = (x1, x2, y1, y2, clipping, ident)

        # Must assume that we have to rerender.
        self.invalidate (INVALIDATE_RENDER)


    # ############################################################
    # State of viewed data item

    def _shownquant_init (self):
        self._shownquant_enable_unusual = False

        get = self.builder.get_object

        ls = Gtk.ListStore (int, str, bool, bool)
        ls.append ((QUANT_REAL, 'Real', True, True))
        ls.append ((QUANT_IMAG, 'Imaginary', True, True))
        ls.append ((QUANT_AMP, 'Amplitude', True, True))
        ls.append ((QUANT_PHA, 'Phase', True, True))
        ls.append ((0, '', True, True))
        ls.append ((QUANT_PHUNCERT, 'Phase Uncert.', False, False))
        ls.append ((QUANT_SYSPHASE, 'System Phase', False, False))

        cb = get ('show_cb')
        cr = Gtk.CellRendererText ()
        cb.pack_start (cr, True)
        cb.add_attribute (cr, 'text', 1)
        cb.add_attribute (cr, 'sensitive', 2)
        def _issep (model, iter, unused=None):
            return model.get_value (iter, 1) == ''
        cb.set_row_separator_func (_issep)
        cb.set_model (ls)
        cb.connect ('changed', self._shownquant_changed)
        self._shownquant_cb = cb

        cb.set_active (2)


    def _shownquant_toggle_unusual (self):
        enabled = not self._shownquant_enable_unusual
        self._shownquant_enable_unusual = enabled

        model = self._shownquant_cb.get_model ()
        iter = model.get_iter_first ()

        while iter is not None:
            always = model.get_value (iter, 3)
            model.set_value (iter, 2, always or enabled)
            iter = model.iter_next (iter)

        return True


    def _shownquant_get (self):
        iter = self._shownquant_cb.get_active_iter ()
        if iter is None:
            return None

        return self._shownquant_cb.get_model ().get_value (iter, 0)


    def _shownquant_changed (self, widget):
        self.darea.queue_draw ()


    def _shownquant_show_next (self):
        iter = self._shownquant_cb.get_active_iter ()
        if iter is None:
            return False

        model = self._shownquant_cb.get_model ()
        iter = model.iter_next (iter)
        while iter is not None:
            if not model.get_value (iter, 2):
                pass
            elif model.get_value (iter, 1) == '':
                pass
            else:
                self._shownquant_cb.set_active_iter (iter)
                return True
            iter = model.iter_next (iter)

        return False


    def _shownquant_show_prev (self):
        active = self._shownquant_cb.get_active_iter ()
        if active is None:
            return False

        model = self._shownquant_cb.get_model ()
        actpath = model.get_path (active)
        iter = model.get_iter_first ()
        lastgood = None

        while iter is not None:
            next = model.iter_next (iter)
            if next is None:
                break

            if not model.get_value (iter, 2):
                pass
            elif model.get_value (iter, 1) == '':
                pass
            else:
                lastgood = iter

            nextpath = model.get_path (next)

            if nextpath == actpath:
                if lastgood is not None:
                    self._shownquant_cb.set_active_iter (lastgood)
                    return True
                return False

            iter = next

        return False


    # ############################################################
    # Coordinates

    def _widget_to_scalar (self, widget, x, y, clamp=False, lattice=False):
        w = widget.get_allocated_width ()
        h = widget.get_allocated_height ()

        sx = 1.0 * x / w
        sy = 1.0 * y / h

        if clamp:
            sx = min (max (sx, 0), 1)
            sy = min (max (sy, 0), 1)

        if lattice:
            # Map onto lattice points of our current grid system
            self._vgrid_update (self._bpsel_get_current ())
            assert self.vgrid is not None

            nt, nf = self.vgrid.data.shape

            sx = round (sx * nf) / nf
            sy = round (sy * nt) / nt

        return sx, sy


    def _event_to_scalar (self, widget, event, **kwargs):
        return self._widget_to_scalar (widget, event.x, event.y, **kwargs)


    def _pointer_to_scalar (self, widget, **kwargs):
        x, y = widget.get_pointer ()
        return self._widget_to_scalar (widget, x, y, **kwargs)


    def _scalar_to_widget (self, widget, sx, sy):
        w = widget.get_allocated_width ()
        h = widget.get_allocated_height ()
        return sx * w, sy * h


    # ############################################################
    # Display area rendering

    def _darea_init (self):
        self.img_rendered_for = None
        self.img_scaledsurf = None
        self.img_scaled_for = None
        self.img_data = None

        self.showall = False
        self.effmaxcoord = None

        self.darea.grab_focus ()


    def _draw (self, widget, ctxt):
        w = widget.get_allocated_width ()
        h = widget.get_allocated_height ()

        bp = self._bpsel_get_current ()
        self._vgrid_update (bp)

        if bp is None:
            # No data, or something
            text = 'No basepols selected'
            extents = ctxt.text_extents (text)
            ctxt.move_to ((w - extents[2]) / 2 - extents[0],
                          (h - extents[3]) / 2 - extents[1])
            ctxt.set_source_rgb (0, 0, 0)
            ctxt.show_text (text)
            return

        # Create the buffer for the image data, if needed
        if self.img_data is None:
            self.img_scaledsurf = self.img_scaled_for = None

            nrec, nchan = self.vgrid.data.shape

            stride = cairo.ImageSurface.format_stride_for_width (cairo.FORMAT_RGB24, nchan)
            imgdata = np.empty ((nrec, (stride + 3) // 4), dtype=np.int32)
            self.img_surf = cairo.ImageSurface.create_for_data (imgdata, cairo.FORMAT_RGB24,
                                                                nchan, nrec, stride)
            self.img_data = imgdata[:,:nchan]

        # Render the data into the image buffer, if needed
        self._img_update_render ()

        # Rescale the rendered image, if needed
        if self.img_scaled_for != (w, h):
            ssurf = cairo.ImageSurface (cairo.FORMAT_RGB24, w, h)
            sctxt = cairo.Context (ssurf)
            sctxt.scale (1.0 * w / self.img_data.shape[1], 1.0 * h / self.img_data.shape[0])
            sctxt.set_source_surface (self.img_surf, 0, 0)
            pat = sctxt.get_source ()
            pat.set_extend (cairo.EXTEND_NONE)
            pat.set_filter (cairo.FILTER_NEAREST)
            sctxt.paint ()

            self.img_scaledsurf = ssurf
            self.img_scaled_for = (w, h)

        # Now actually draw.
        ctxt.save ()
        ctxt.set_source_surface (self.img_scaledsurf, 0, 0)
        ctxt.paint ()
        ctxt.restore ()

        if bp is None:
            return True

        ctxt.set_line_width (1)

        if self.effmaxcoord is not None:
            sx, sy = self.vgrid.physical_to_scalar (*self.effmaxcoord)

            if sx >= 0 and sx <= 1 and sy >= 0 and sy <= 1:
                ctxt.arc (sx * w, sy * h, 8, 0, 2 * np.pi)
                ctxt.set_source_rgba (1., 1., 1., 1.)
                ctxt.stroke ()

        if self.inboxsel:
            for x1, x2, y1, y2, clipping, ident in self.boxsel_boxes:
                x1 *= w
                x2 *= w
                y1 *= h
                y2 *= h

                ctxt.rectangle (x1, y1, x2 - x1, y2 - y1)
                if self.boxsel_shade:
                    ctxt.set_source_rgba (0.8, 0.2, 0.2, 0.3)
                    ctxt.fill_preserve ()
                ctxt.set_source_rgba (0.8, 0.2, 0.2, 0.7)
                ctxt.stroke ()

                if self.flagapi.get (ident)[4] is None:
                    # An all-bp record
                    ctxt.move_to (x1, y1)
                    ctxt.line_to (x2, y2)
                    ctxt.stroke ()
                    ctxt.move_to (x2, y1)
                    ctxt.line_to (x1, y2)
                    ctxt.stroke ()

        if self.boxsel_selected is not None:
            x1, x2, y1, y2, clipping, ident = self.boxsel_selected

            x1 *= w
            x2 *= w
            y1 *= h
            y2 *= h

            ctxt.rectangle (x1, y1, x2 - x1, y2 - y1)
            if self.boxsel_shade:
                ctxt.set_source_rgba (0.9, 0.1, 0.1, 0.5)
                ctxt.fill_preserve ()
            ctxt.set_source_rgba (0.9, 0.1, 0.1, 0.8)
            ctxt.set_line_width (1)
            ctxt.stroke ()
            if not self.boxsel_shade:
                mx = 0.5 * (x1 + x2)
                my = 0.5 * (y1 + y2)
                ctxt.move_to (x1, my)
                ctxt.line_to (mx, y1)
                ctxt.line_to (x2, my)
                ctxt.line_to (mx, y2)
                ctxt.line_to (x1, my)
                ctxt.stroke ()

        if self.selend is not None:
            sx, sy = self.vgrid.physical_to_scalar (*self.selstart)
            ex, ey = self.vgrid.physical_to_scalar (*self.selend)

            sx *= w
            sy *= h
            ex *= w
            ey *= h

            ctxt.rectangle (sx, sy, ex - sx, ey - sy)
            ctxt.set_source_rgba (0.2, 0.2, 0.8, 0.3)
            ctxt.fill_preserve ()
            ctxt.set_source_rgba (0.2, 0.2, 0.8, 0.7)
            ctxt.stroke ()

        return True


    def _img_update_render (self):
        bp = self._bpsel_get_current ()
        shownquant = self._shownquant_get ()

        if bp is None:
            rkey = None
        else:
            rkey = (bp, shownquant, self.showall)

        if self.img_rendered_for == rkey:
            return

        img = self.img_data
        img.fill (0)
        self.img_rendered_for = rkey

        if shownquant is None or bp is None:
            return

        renderer = _rendererMakers[shownquant] (self.vgrid_base)
        renderer (self.tdata, bp, self.vgrid, img, self.showall)
        self.invalidate (INVALIDATE_SCALE)


    def _toggle_showall_mode (self):
        self.showall = not self.showall
        self.darea.queue_draw ()


    # ############################################################
    # Keyboard commands

    def _on_darea_key_press (self, widget, event):
        kn = Gdk.keyval_name (event.keyval)

        if kn == 'Control_L' or kn == 'Control_R':
            self._tooltip_update (widget, event)
        return False


    def _on_darea_key_release (self, widget, event):
        kn = Gdk.keyval_name (event.keyval)
        tooltip = False
        handled = False

        if kn == 'Control_L' or kn == 'Control_R':
            tooltip = False
            handled = True
        elif event.state & Gdk.ModifierType.CONTROL_MASK:
            tooltip = True
        else:
            tooltip = False

        if not tooltip:
            self._tooltip_destroy ()

        return False


    def _on_keypress (self, widget, event):
        kn = Gdk.keyval_name (event.keyval)
        modmask = Gtk.accelerator_get_default_mod_mask ()
        isctrl = (event.state & modmask) == Gdk.ModifierType.CONTROL_MASK

        if kn == 'q':
            if isctrl:
                # Close the window
                self.close ()
                return True
            rv = self._bpsel_select_prev ()
        elif kn == 'w':
            if isctrl:
                # Close the window
                self.close ()
                return True
            rv = self._bpsel_select_next ()
        elif kn == 'f':
            rv = self._bpsel_select_first ()
        elif kn == 'l':
            rv = self._bpsel_select_last ()
        elif kn == 'j':
            rv = self._shownquant_show_prev ()
        elif kn == 'k':
            rv = self._shownquant_show_next ()
        elif kn == 'U':
            rv = self._shownquant_toggle_unusual ()
        elif kn == 'v':
            rv = self._sel_span_vert ()
        elif kn == 'h':
            rv = self._sel_span_horz ()
        elif kn == 'z':
            rv = self._zoom_in ()
        elif kn == 'Z':
            rv = self._zoom_out ()
        elif kn == 'a':
            rv = self._toggle_showall_mode ()
        elif kn == 't':
            rv = self._xform_toggle_tempskip ()
        elif kn == 'p':
            rv = self._sel_print ()
        elif kn == '1':
            rv = self._resize_win (1, 1)
        elif kn == '2':
            rv = self._resize_win (2, 2)
        elif kn == '3':
            rv = self._resize_win (1, 2)
        else:
            #print 'unhandled key:', kn
            rv = False

        if self.flagapi is not None:
            if kn == 'b':
                rv = self._sel_flag_allbps ()
            elif kn == 'c':
                rv = self._sel_flag_curbp ()
            elif kn == 'm':
                rv = self._toggle_effmax ()
            elif kn == 's':
                rv = self._flag_curbp ()
            elif kn == 'u':
                rv = self._boxsel_update ()
            elif kn == 'e':
                rv = self._toggle_boxsel_mode ()
            elif kn == 'E':
                rv = self._toggle_boxsel_shading ()
            elif kn == 'space':
                rv = self._boxsel_toggle ()
            elif kn == 'Delete':
                rv = self._boxsel_delete ()
            elif kn == 'M':
                rv = self._magic_rfi_flag ()

        if rv:
            self.darea.grab_focus ()
        return rv


    def _toggle_boxsel_mode (self):
        self.inboxsel = not self.inboxsel
        if self.inboxsel:
            self._compute_boxsel_boxes ()
        else:
            self.boxsel_boxes = None
            self._boxsel_set_selection (None)
        self.darea.queue_draw ()
        return True


    def _toggle_boxsel_shading (self):
        self.boxsel_shade = not self.boxsel_shade
        self.darea.queue_draw ()
        return True


    def _boxsel_compute_boxinfo (self, tst, tend, fst, fend):
        clipping = [False] * 4

        # Clamp the scalar coordinates ourselves, recording which box
        # edges go beyond our viewport so we know which ones to make
        # un-editable and to preserve if the box is edited. Coordinate
        # transforms can add some floating-point jitter so allow a bit
        # of slop.

        if fst is None:
            x1, x2 = 0., 1.
            clipping[RESIZE_LEFT] = clipping[RESIZE_RIGHT] = True
        else:
            x1, x2 = self.vgrid.faxis.valsToScalars([fst, fend])

            if x2 < x1:
                x1, x2 = x2, x1

            if x1 < 0:
                if x1 < -1e-5:
                    clipping[RESIZE_LEFT] = True
                x1 = 0.

            if x2 > 1:
                if x2 > 1.00001:
                    clipping[RESIZE_RIGHT] = True
                x2 = 1.

        if tst is None:
            y1, y2 = 0., 1.
            clipping[RESIZE_TOP] = clipping[RESIZE_BOTTOM] = True
        else:
            y1, y2 = self.vgrid.taxis.valsToScalars([tst, tend])

            if y2 < y1:
                y1, y2 = y2, y1

            if y1 < 0:
                if y1 < -1e-5:
                    clipping[RESIZE_TOP] = True
                y1 = 0.

            if y2 > 1:
                if y2 > 1.00001:
                    clipping[RESIZE_BOTTOM] = True
                y2 = 1.

        return x1, x2, y1, y2, clipping


    def _compute_boxsel_boxes (self):
        boxes = []
        self._vgrid_update (self._bpsel_get_current ())

        def addbox (ident, tst, tend, fst, fend):
            if tst is None and fst is None:
                # universal basepol flagging, ignore for boxsel
                return

            x1, x2, y1, y2, clipping = self._boxsel_compute_boxinfo (tst, tend,
                                                                     fst, fend)
            boxes.append ((x1, x2, y1, y2, clipping, ident))

        for k in self.flagapi.getRecords (self.vgrid.taxis,
                                          self.vgrid.faxis):
            addbox (*k)

        self.boxsel_boxes = boxes


    def _boxsel_get_hits (self, boxlist, sx, sy, allok):
        def contains (tup):
            x1, x2, y1, y2, clipping, ident = tup
            if not allok:
                # Check whether the bpmatch is None = match-all
                if self.flagapi.get (ident)[4] is None:
                    return False
            if sx < x1 or sx >= x2:
                return False
            if sy < y1 or sy >= y2:
                return False
            return True

        def area (tup):
            return (tup[1] - tup[0]) * (tup[3] - tup[2])

        return sorted ((t for t in boxlist if contains (t)),
                       key=area)


    def _boxsel_sel_coords (self, sx, sy, allok):
        hits = self._boxsel_get_hits (self.boxsel_boxes, sx, sy, allok)

        if len (hits) == 0:
            return None

        if self.boxsel_selected is None or self.boxsel_selected not in hits:
            selarea = 0
        else:
            selarea = (self.boxsel_selected[1] - self.boxsel_selected[0]) * \
                (self.boxsel_selected[3] - self.boxsel_selected[2])

        for i in xrange (len (hits)):
            hitarea = (hits[i][1] - hits[i][0]) * (hits[i][3] - hits[i][2])

            if hitarea > selarea:
                return hits[i]
            if hitarea < selarea:
                continue
            if hits[i] != self.boxsel_selected:
                return hits[i]
            # FIXME: this will fail if there are two boxes with
            # identical areas.

        # The biggest box had been selected, cycle around to no selection.
        return None


    def _boxsel_set_selection (self, newsel):
        if self.boxsel_selected == newsel:
            return

        self.boxsel_selected = newsel

        if newsel is None:
            if self.boxsel_resizemode != RESIZE_NONE:
                self.boxsel_resizemode = RESIZE_NONE

                if self.grabbed:
                    self.darea.grab_remove ()
                    self.grabbed = False

                self.darea.get_window ().set_cursor (None)

        self.darea.queue_draw ()


    def _boxsel_toggle (self):
        bp = self._bpsel_get_current ()

        if not self.inboxsel:
            # There was some code here to toggle the flagging of boxsel_selected
            # if it was defined and a non-all box. I'm not sure how that was
            # supposed to happen and couldn't figure out how to trigger that
            # condition in the UI; so, I deleted it.
            return False

        sx, sy = self._pointer_to_scalar (self.darea)

        if sx < 0 or sy < 0 or sx > 1 or sy > 1:
            return False

        if self.boxsel_selected is not None:
            hits = self._boxsel_get_hits ([self.boxsel_selected], sx, sy, True)
            if len (hits) > 0:
                # If still in previously-selected box, probably doing
                # a repeated toggle, and we should unsel it.
                self.flagapi.toggleBasepol (self.boxsel_selected[5], bp)
                self.flagapi.save ()
                self.invalidate (INVALIDATE_FLAGS)


        self._boxsel_set_selection (self._boxsel_sel_coords (sx, sy, False))

        if self.boxsel_selected is not None:
            self.flagapi.toggleBasepol (self.boxsel_selected[5], bp)
            self.flagapi.save ()
            self.invalidate (INVALIDATE_FLAGS)

        return True


    def _boxsel_resize_motion (self, mapx, mapy):
        rm = self.boxsel_resizemode
        x1, x2, y1, y2, clipping, ident = self.boxsel_selected

        if rm == RESIZE_TOP:
            if mapy <= y2:
                y1 = mapy
            else:
                # We have moved to below the bottom edge of the original box;
                # switch to RESIZE_BOTTOM mode
                self.boxsel_resizemode = RESIZE_BOTTOM
                y1 = y2
                y2 = mapy
        elif rm == RESIZE_BOTTOM:
            if mapy >= y1:
                y2 = mapy
            else:
                # Moved above top edge.
                self.boxsel_resizemode = RESIZE_TOP
                y2 = y1
                y1 = mapy
        elif rm == RESIZE_LEFT:
            if mapx <= x2:
                x1 = mapx
            else:
                # moved beyond right edge
                self.boxsel_resizemode = RESIZE_RIGHT
                x1 = x2
                x2 = mapx
        elif rm == RESIZE_RIGHT:
            if mapx >= x1:
                x2 = mapx
            else:
                # moved to left of left edge
                self.boxsel_resizemode = RESIZE_LEFT
                x2 = x1
                x1 = mapx

        self.boxsel_selected = (x1, x2, y1, y2, clipping, ident)

        darea = self.darea

        if self.boxsel_resizemode != rm:
            # Updated direction.
            curs = Gdk.Cursor.new_from_name (darea.get_display (),
                                                 _resizeCursors[self.boxsel_resizemode])
            darea.get_window ().set_cursor (curs)

        self.darea.queue_draw ()


    def _boxsel_resize_end (self, widget, event):
        sx, sy = self._event_to_scalar (widget, event, clamp=True, lattice=True)
        self._boxsel_resize_motion (sx, sy)
        self.boxsel_resizemode = RESIZE_NONE
        self.darea.get_window ().set_cursor (None)


    def _boxsel_update (self):
        if self.boxsel_selected is None:
            return False

        x1, x2, y1, y2, clipping, ident = self.boxsel_selected

        # Compute new physical coordinates of the box, taking care
        # not to modify any edges whose locations are clipped out
        # of our viewport.

        t1, f1 = self.vgrid.scalar_to_physical (x1, y1)
        t2, f2 = self.vgrid.scalar_to_physical (x2, y2)
        newbounds = [t1, t2, f1, f2]
        oldbounds = self.flagapi.get (ident)[:4]

        if clipping[RESIZE_LEFT]: # left <-> f1
            newbounds[2] = oldbounds[2]
        if clipping[RESIZE_RIGHT]: # right <-> f2
            newbounds[3] = oldbounds[3]
        if clipping[RESIZE_TOP]: # top <-> t1
            newbounds[0] = oldbounds[0]
        if clipping[RESIZE_BOTTOM]: # bottom <-> t2
            newbounds[1] = oldbounds[1]

        # Apply the change.

        self.flagapi.resize (ident, *newbounds)
        self.flagapi.save ()
        self.invalidate (INVALIDATE_FLAGS)
        return True


    def _boxsel_delete (self):
        if self.boxsel_selected is None:
            return False

        x1, x2, y1, y2, clipping, ident = self.boxsel_selected
        self.flagapi.delete (ident)
        self.flagapi.save ()
        self.invalidate (INVALIDATE_FLAGS)
        self._boxsel_set_selection (None)


    def _flag_curbp (self):
        bp = self._bpsel_get_current ()
        self.flagapi.add (None, None, None, None, set ((bp, )))
        self._bpsel_select_next_or_different ()
        self._bpsel_update () # get old bps out of the list
        self.flagapi.save ()
        self.invalidate (INVALIDATE_FLAGS)
        return True


    def _toggle_effmax (self):
        if self.effmaxcoord is None:
            return False

        if not self.inboxsel:
            self._compute_boxsel_boxes ()

        x, y = self.vgrid.physical_to_scalar (*self.effmaxcoord)
        hits = self._boxsel_get_hits (self.boxsel_boxes, x, y, False)

        if not self.inboxsel:
            self.boxsel_boxes = None

        if len (hits) < 1:
            return False

        bp = self._bpsel_get_current ()
        self.flagapi.toggleBasepol (hits[0][5], bp)
        self.flagapi.save ()
        self.invalidate (INVALIDATE_FLAGS)
        return True


    def _magic_rfi_flag (self):
        if not self.inboxsel:
            return False
        if self.boxsel_selected is None:
            return False

        import sys
        sys.path = ['.'] + sys.path

        try:
            import rfitest
            rfitest = reload (rfitest)
        except ImportError:
            print ('Cannot magic flag: no "rfitest.py" in current directory')
            return False

        x1, x2, y1, y2, clipping, ident = self.boxsel_selected
        t1, t2, f1, f2, bpmatch = self.flagapi.get (ident)

        if bpmatch is None:
            print ('Cannot magic flag: selected box applies to all basepols')
            return False

        print ('Magic flag!!!')
        newbpmatch = set ()
        sst = SubsetTransform ()
        sst.setBounds (t1, t2, f1, f2)
        sst.setParent (self.tdata_base)

        try:
            newbpmatch = rfitest.decide (self.bpsel_list, sst.getGrid)
        except Exception, e:
            print ('*** Exception raise in magic flagging routine!')
            import traceback
            traceback.print_exc ()
            print ('*** Ignoring magicflag results')
            return

        if newbpmatch is None:
            print ('Done. No flagging information returned.')
            return

        if not len (newbpmatch):
            print ('   Not updating flags since nothing was selected')
            print ('   The logical equivalent would be to delete the region')
        else:
            self.flagapi.setBPMatch (ident, newbpmatch)
        self.flagapi.save ()
        self.invalidate (INVALIDATE_FLAGS)


    # ############################################################
    # Data tooltip

    def _tooltip_init (self):
        self.tooltip_shown = False
        self.tooltip_window = None
        self.tooltip_label = None


    def _tooltip_update (self, widget, event):
        if hasattr (event, 'x'):
            sx, sy = self._event_to_scalar (widget, event)
            lx, ly = self._event_to_scalar (widget, event, clamp=True, lattice=True)
        else:
            # When we're created by a key-press-event, there's
            # no position information by default.
            sx, sy = self._pointer_to_scalar (self.darea)
            lx, ly = self._pointer_to_scalar (self.darea, clamp=True, lattice=True)

        if sx < 0 or sy < 0 or sx > 1 or sy > 1:
            # Out of bounds - can happen with sloppy focus
            self._tooltip_destroy ()
            return

        if self.tooltip_shown:
            label = self.tooltip_label
        else:
            win = Gtk.Window (Gtk.WindowType.POPUP)
            win.set_type_hint (Gdk.WindowTypeHint.TOOLTIP)
            #win.set_gravity (Gdk.GRAVITY_NORTH_WEST)
            win.set_transient_for (self.win)
            win.set_position (Gtk.WindowPosition.CENTER)

            label = self.tooltip_label = Gtk.Label ()
            win.add (label)
            win.show_all ()

            self.tooltip_window = win
            self.tooltip_shown = True

        nt, nf = self.vgrid.data.shape
        i, j = int (np.floor (sy * nt)), int (np.floor (sx * nf))
        flag = self.vgrid.flags[i,j]
        datum = self.vgrid.data[i,j]
        amp = np.abs (datum)
        pha = np.arctan2 (datum.imag, datum.real) * 180 / np.pi

        time, freq = self.vgrid.scalar_to_physical (sx, sy)
        tstr = jdToFull (time)
        f0 = self.vgrid.faxis.start
        uvw0 = self.vgrid.uvw[i]
        uvw = uvw0 * freq # ns to wavelengths
        uvd = np.sqrt (uvw[0]**2 + uvw[1]**2)

        # Given noise, amp is a biased estimator of the true
        # visibility amplitude, so our assessment of the phase
        # noise is a bit off. The maximum-likelihood estimate
        # of the *actual* visibility amplitude is sqrt(amp^2 - sigma^2),
        # cf. Thompson, Moran, & Swenson section 9.3. Whatever.

        if self.vgrid.weights[i] == 0:
            sigma = np.inf
            sigph = phaseUncert (0) * 180 / np.pi
        else:
            sigma = self.vgrid.weights[i]**-0.5
            sigph = phaseUncert (amp / sigma) * 180 / np.pi

        bnt, bnf = self.vgrid_base.data.shape
        bsx, bsy = self.vgrid_base.physical_to_scalar (time, freq)
        bi, bj = int (np.floor (bsy * bnt)), int (np.floor (bsx * bnf))

        text = """<tt><b>%8.4f</b> amp <b>%8.0f</b>°pha
<b>%8.4f</b> re  <b>%8.4f</b> im
<b>%8.5f</b> σ   <b>%8.2f</b>°σφ
""" % (amp, pha, datum.real, datum.imag, sigma, sigph)

        text += """

<b>%s</b> UT, <b>%.5f</b> GHz
uvw/dist = <b>%.0f %.0f %.0f</b> / <b>%.0f</b> λ
<b>%4d</b> row <b>%4d</b> col</tt>""" % (tstr, freq,
                                         uvw[0], uvw[1], uvw[2], uvd,
                                         bi, bj)

        # TODO: HA, LST, parang, fringerate, ant{1,2}{az,el}
        label.set_markup (text)


    def _tooltip_destroy (self):
        if not self.tooltip_shown:
            return

        self.tooltip_window.destroy ()
        self.tooltip_window = None
        self.tooltip_label = None
        self.tooltip_shown = False


    # ############################################################
    # Selection management

    def _selection_init (self):
        self.selstart = None
        self.selend = None
        self.selflags = 0
        self.grabbed = False
        self.inboxsel = False
        self.sel_had_prev = False
        self.boxsel_selected = None
        self.boxsel_resizemode = RESIZE_NONE
        self.boxsel_shade = True


    def _button_press (self, widget, event):
        if self._bpsel_get_current () is None:
            return True

        if event.type == Gdk.EventType.BUTTON_PRESS and event.button == 1:
            widget.grab_focus ()

            if self.boxsel_resizemode != RESIZE_NONE:
                self.grabbed = True
                widget.grab_add ()
            else:
                # Record whether we just started with a selection to distinguish
                # between selection-clearing single left-clicks and shownquant-changing
                # single left-clicks.
                self.sel_had_prev = self.selend is not None
                sxy = self._event_to_scalar (widget, event, clamp=True, lattice=True)
                self.selstart = self.vgrid.scalar_to_physical (*sxy)
                self.selend = None
                self.selflags = 0
                self.grabbed = True
                widget.grab_add ()
            return True

        return False


    def _button_release (self, widget, event):
        if self._bpsel_get_current () is None:
            return True

        if self.grabbed:
            widget.grab_remove ()
            self.grabbed = False

        if event.button == 3:
            # (Presumed single) right-click
            self._shownquant_show_next ()
            return True

        if event.button == 1:
            # We'll always do something with this click that requires a redraw.
            self.darea.queue_draw ()

            if self.boxsel_resizemode != RESIZE_NONE:
                # Highest-priority: finish resizing a box
                self._boxsel_resize_end (widget, event)
                return True

            exy = self._event_to_scalar (widget, event, clamp=True, lattice=True)
            self.selend = self.vgrid.scalar_to_physical (*exy)

            if self.selend == self.selstart:
                # We just a single left-click. Don't set up drawn selection.
                self.selstart = self.selend = None
                self.selflags = 0
            else:
                # Next-highest priority: was not a single click;
                # finish drawing the selection
                return True

            if self.sel_had_prev:
                # Next-highest: user was just clicking to clear previous drawn
                # selection
                return True

            if self.inboxsel:
                # Next-highest: select new box, or clear box selection
                sx, sy = self._event_to_scalar (widget, event)
                newsel = self._boxsel_sel_coords (sx, sy, True)

                if newsel is not None or self.boxsel_selected is not None:
                    self._boxsel_set_selection (newsel)
                    return True

            # Lowest priority: cycle shown quantity
            self._shownquant_show_prev ()
            return True

        return False


    def _motion_notify (self, widget, event):
        if self._bpsel_get_current () is None:
            return True
        if self.vgrid is None:
            # Can happen if the window pops open on top of the cursor
            return False

        sx, sy = self._event_to_scalar (widget, event)
        lx, ly = self._event_to_scalar (widget, event, clamp=True, lattice=True)

        if self.grabbed:
            if self.boxsel_resizemode != RESIZE_NONE:
                self._boxsel_resize_motion (lx, ly)
            else:
                self.selend = self.vgrid.scalar_to_physical (lx, ly)
                self.darea.queue_draw ()
            return True

        if self.boxsel_selected is not None:
            x1, x2, y1, y2, clipping, ident = self.boxsel_selected

            in_horz = sy >= y1 and sy <= y2
            in_vert = sx >= x1 and sx <= x2

            # proximity testing should be done in pixel-space
            wx1, wy1 = self._scalar_to_widget (widget, x1, y1)
            wx2, wy2 = self._scalar_to_widget (widget, x2, y2)

            rm = RESIZE_NONE

            if in_horz and abs (event.x - wx1) < 4:
                rm = RESIZE_LEFT
            elif in_horz and abs (event.x - wx2) < 4:
                rm = RESIZE_RIGHT
            elif in_vert and abs (event.y - wy1) < 4:
                rm = RESIZE_TOP
            elif in_vert and abs (event.y - wy2) < 4:
                rm = RESIZE_BOTTOM

            if self.boxsel_resizemode != rm:
                if rm == RESIZE_NONE:
                    curs = None
                else:
                    try:
                        curs = Gdk.Cursor.new_from_name (widget.get_display (), _resizeCursors[rm])
                    except TypeError:
                        curs = None # HACK: cursor lookup busted on Anaconda Gtk install?

                widget.get_window ().set_cursor (curs)
                self.boxsel_resizemode = rm

        if event.state & Gdk.ModifierType.CONTROL_MASK:
            self._tooltip_update (widget, event)

        return True


    def _on_darea_leave_notify (self, widget, event):
        # Noop if no tooltip.
        self._tooltip_destroy ()
        return True


    def _sel_span_vert (self):
        if self.selend is None:
            return False

        self.selstart = (self.vgrid.taxis.scalarsToVals (0), self.selstart[1])
        self.selend = (self.vgrid.taxis.scalarsToVals (1), self.selend[1])
        self.selflags |= SELFLAG_VSPAN
        self.darea.queue_draw ()
        return True


    def _sel_span_horz (self):
        if self.selend is None:
            return False

        self.selstart = (self.selstart[0], self.vgrid.faxis.scalarsToVals (0))
        self.selend = (self.selend[0], self.vgrid.faxis.scalarsToVals (1))
        self.selflags |= SELFLAG_HSPAN
        self.darea.queue_draw ()
        return True


    def _sel_canonical (self, spans=False):
        if self.selend is None:
            return None

        st, sf = self.selstart
        et, ef = self.selend

        if st == et or sf == ef:
            return None

        if et < st:
            st, et = et, st
        if ef < sf:
            sf, ef = ef, sf

        if spans:
            if self.selflags & SELFLAG_VSPAN:
                st, et = None, None
            if self.selflags & SELFLAG_HSPAN:
                sf, ef = None, None

        return st, et, sf, ef


    def _sel_print (self):
        if self.selend is None:
            return False

        bp = self._bpsel_get_current ()

        print ('Current selection:')
        print ('   ', 'basepol', self.format_bp (bp))

        t1, t2, f1, f2 = self._sel_canonical (spans=True)
        if t1 is None and f1 is None:
            print ('   ', 'match-nothing')
        elif t1 is None:
            print ('   ', 'freq(%.18e,%.18e)' % (f1, f2))
        elif f1 is None:
            print ('   ', 'time(%s,%s)' % (jdToFull (t1), jdToFull (t2)))
        else:
            print ('   ', 'time(%s,%s)' % (jdToFull (t1), jdToFull (t2)),
                   'freq(%.18e,%.18e)' % (f1, f2))

        subset = SubsetTransform ()
        subset.setBounds (t1, t2, f1, f2)
        subset.setParent (self.tdata)
        sgrid = subset.getGrid (bp)
        d = sgrid.ffdata ()
        m = d.mean ()
        ph = np.arctan2 (m.imag, m.real) * 180 / np.pi
        print ('   ', 'mean %.6e @ %.1f°' % (np.abs (m), ph))
        s = 0.5 * (d.real.std () + d.imag.std ())
        print ('   ', 'stddev %.6e' % s)
        return True


    def _sel_flag_allbps (self):
        if self.selend is not None:
            t1, t2, f1, f2 = self._sel_canonical (spans=True)
            self.flagapi.add (t1, t2, f1, f2, None)
            self.flagapi.save ()
            self.selstart = self.selend = None
            self.selflags = 0
            self.invalidate (INVALIDATE_FLAGS)
            return True

        if self.boxsel_selected is not None:
            # Upgrade selected box to apply to all bps?
            # This logic used to only fire if the box was an all-bp box.
            # Changed to fire even if it's already all-bp
            x1, x2, y1, y2, clipping, ident = self.boxsel_selected
            self.flagapi.setBPMatch (ident, None)
            self.flagapi.save ()
            self.invalidate (INVALIDATE_FLAGS)
            return True

        return False


    def _sel_flag_curbp (self):
        if self.selend is None:
            return False

        bp = self._bpsel_get_current ()
        t1, t2, f1, f2 = self._sel_canonical (spans=True)
        self.flagapi.add (t1, t2, f1, f2, set ((bp, )))
        self.flagapi.save ()
        self.selstart = self.selend = None
        self.selflags = 0
        self.invalidate (INVALIDATE_FLAGS)
        return True


    # ############################################################
    # Zooming

    def _zoom_init (self):
        self._zoom_bounds = None
        self._zoom_transform = SubsetTransform ()


    def _zoom_in (self):
        if self.selend is None:
            return False

        self._zoom_bounds = self._sel_canonical ()
        self.selstart = self.selend = None
        self.selflags = 0
        self.invalidate (INVALIDATE_VGRID)
        return True


    def _zoom_out (self):
        self._zoom_bounds = None
        self.invalidate (INVALIDATE_VGRID)
        return True


    def _resize_win (self, widthfactor, heightfactor):
        ww = self.win.get_allocated_width ()
        wh = self.win.get_allocated_height ()
        dw = self.darea.get_allocated_width ()
        dh = self.darea.get_allocated_height ()

        pad_w = ww - dw
        pad_h = wh - dh

        neww = pad_w + self.vgrid.data.shape[1] * widthfactor
        newh = pad_h + self.vgrid.data.shape[0] * heightfactor
        self.win.resize (neww, newh)
        return True

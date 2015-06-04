# -*- mode: python; coding: utf-8 -*-
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""pwkit.scanalyzer - interactive analysis of interferometric visbilities.

This uses the Gtk+3 graphical toolkit. This is inherited from some old
MIRIAD-targeted code so it's pretty messy.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = (b'get_transposed launch').split ()

from os.path import exists, join

from pwkit import PKError


def get_transposed (path, transpose_args):
    from . import transpose
    import os.path

    if exists (join (path, 'table.dat')):
        # It's a MeasurementSet!! Hacks.
        import hashlib
        h = hashlib.sha1 ()
        h.update (path)
        h.update ('args')
        for name in sorted (transpose_args.iterkeys ()):
            if name[0] != '_':
                h.update (name)
                h.update ('1')
                h.update (transpose_args[name])
                h.update ('2')
        vhash = h.digest ()
        tfunc = transpose.ms_transpose
    else:
        # TODO: put MIRIAD support back in.
        raise PKError ('unsupported dataset %r (MIRIAD support has been disabled)', path)

    tpath = 'transpose.' + vhash.encode ('hex') + '.dat'

    if exists (tpath):
        print ('Loading', tpath, '...')
    else:
        print ('Transposing', path, 'into', tpath, '...')
        from time import time
        tst = time ()
        nrec, nvis, nbyte = tfunc (path, tpath, transpose_args)
        elapsed = time () - tst
        print ('  Total time elapsed: %.2f s' % elapsed)
        print ('  Record rate: %.0f/s' % (nrec / elapsed))
        print ('  Visibility rate: %.0f/s' % (nvis / elapsed))
        print ('  Write rate: %.0f kiB/s' % (nbyte / elapsed / 1024))

    return transpose.TransposeFile (open (tpath, 'r+'))


def launch (path, flagpath=None, transpose_args={}):
    tfile = get_transposed (path, transpose_args)

    if flagpath is None:
        flags = None
    else:
        from .flag import FlagImplementation
        flags = FlagImplementation (flagpath)
        flags.try_load ()

    from gi.repository import Gtk
    from .ui import Scanalyzer

    sa = Scanalyzer (tfile, flags)
    sa.win.show_all ()
    sa.win.connect ('destroy', Gtk.main_quit)

    Gtk.main ()

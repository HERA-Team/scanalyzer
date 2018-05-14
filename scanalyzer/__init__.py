# -*- mode: python; coding: utf-8 -*-
# Copyright 2015-2016 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""scanalyzer - interactive analysis of interferometric visbilities.

This uses the Gtk+3 graphical toolkit. This is inherited from some old
MIRIAD-targeted code so it's pretty messy.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = str('get_transposed launch').split ()

from os.path import exists, join
from six import iterkeys

from pwkit import PKError


def get_transposed (path, transpose_args):
    from . import transpose
    import os.path

    if exists (join (path, 'table.dat')):
        # It's a MeasurementSet!! Hacks.
        import hashlib
        h = hashlib.sha1 ()
        h.update (path.encode('utf8'))
        h.update (b'args')
        for name in sorted (iterkeys(transpose_args)):
            if name[0] != '_':
                h.update (name.encode('utf8'))
                h.update (b'1')
                h.update (transpose_args[name].encode('utf8'))
                h.update (b'2')
        vhash = h.digest ()
        tfunc = transpose.ms_transpose
    else:
        # HAAACK due to disappearance of my quickHash() function.
        import hashlib
        h = hashlib.sha1 ()
        h.update (path)
        vhash = h.digest ()
        tfunc = transpose.mir_transpose

    tpath = 'transpose.' + vhash.hex() + '.dat'

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

    return transpose.TransposeFile (open (tpath, 'rb+'))


def launch (path, flagpath=None, transpose_args={}):
    tfile = get_transposed (path, transpose_args)

    # NOTE: if flagpath is None, FlagImplementation does stuff but will treat
    # save/load as noops.
    from .flag import FlagImplementation
    flags = FlagImplementation (flagpath)
    flags.try_load ()

    from gi.repository import Gtk
    from .ui import Scanalyzer

    sa = Scanalyzer (tfile, flags)
    sa.win.show_all ()
    sa.win.connect ('destroy', Gtk.main_quit)

    Gtk.main ()

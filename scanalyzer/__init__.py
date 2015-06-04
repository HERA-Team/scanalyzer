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

from .. import PKError


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


# The CLI interface

from ..cli import multitool

class Go (multitool.Command):
    name = 'go'
    argspec = '<vis> [flagfile] [kw1=val1] [kw2=val2...]'
    summary = 'Interactively view and edit visibility data.'
    more_help = """Keyword arguments are used to filter the data when transposing."""
    help_if_no_args = True

    def invoke (self, args, **kwargs):
        vis = args[0]
        flags = None
        kws = {}

        for rest in args[1:]:
            if '=' in rest:
                k, v = rest.split ('=', 1)
                kws[k] = v
            elif flags is not None:
                raise multitool.UsageError ('expected 1 or 2 non-keyword '
                                            'arguments; got at least 3')
            else:
                flags = rest

        launch (vis, flags, kws)


class Tinfo (multitool.Command):
    name = 'tinfo'
    argspec = '<transpose-file>'
    summary = 'Print out information about a transpose file.'

    def invoke (self, args, **kwargs):
        if len (args) != 1:
            raise multitool.UsageError ('expected exactly one argument')

        from . import transpose

        with open (args[0]) as f:
            tf = transpose.TransposeFile (f)

            if 'vispath' not in tf.vars:
                print ('This file does not record the name of the dataset '
                       'from which it was generated.')
            else:
                print ('Generated from:', tf.vars['vispath'].tostring ())

            if 'transargs' in tf.vars:
                print ('     with args:', tf.vars['transargs'].tostring ())

            print ('Number of basepols:', tf.nbp)
            print ('Number of time slots:', tf.nt)
            print ('Number of freq slots:', tf.nchan)
            print ('Visibility data size:', tf.slice_size * tf.nbp, 'bytes')


from ..cli.multitool import HelpCommand

class Scanalyzer (multitool.Multitool):
    cli_name = 'scanalyzer'
    summary = 'Interactively view and edit interferometric visibility data.'

def commandline ():
    multitool.invoke_tool (globals ())

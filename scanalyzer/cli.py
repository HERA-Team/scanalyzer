# -*- mode: python; coding: utf-8 -*-
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""scanalyzer.cli - the command-line interface to the scanalyzer

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = b'commandline'.split ()

from pwkit.cli import multitool

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

        from . import launch
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


from pwkit.cli.multitool import HelpCommand

class Scanalyzer (multitool.Multitool):
    cli_name = 'scanalyzer'
    summary = 'Interactively view and edit interferometric visibility data.'


def commandline ():
    multitool.invoke_tool (globals ())

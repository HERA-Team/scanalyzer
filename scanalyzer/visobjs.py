# -*- mode: python; coding: utf-8 -*-
# Copyright 2016 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""Utility objects for tracking state in MIRIAD datasets.

(Extracted from PKGW's "arf" package.)

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = str ('Window FreqConfig Pointing').split ()

import numpy as np


class Window (object):
    """A Window is a region of spectrum.  It has three
    characteristics: a center frequency, cfreq, measured in GHz; a
    frequency width, width, measured in GHz, and a number of channels,
    nchan, an integer greater than zero."""

    __slots__ = ['cfreq', 'width', 'nchan']

    def __init__ (self, cfreq=0., width=0., nchan=0):
        self.cfreq = cfreq
        self.width = width
        self.nchan = nchan


    def __eq__ (self, other):
        if not isinstance (other, Window):
            return False
        if self.cfreq != other.cfreq:
            return False
        if self.width != other.width:
            return False
        if self.nchan != other.nchan:
            return False
        return True


    def __ne__ (self, other):
        return not self.__eq__ (other)


    def __nonzero__ (self):
        return self.nchan != 0


    def __hash__ (self):
        return hash (self.cfreq) ^ hash (self.width) ^ hash (self.nchan)


    def __cmp__ (self, other):
        if not isinstance (other, Window):
            return -1
        if self.cfreq < other.cfreq:
            return -1
        if self.cfreq > other.cfreq:
            return 1
        if self.width < other.width:
            return -1
        if self.width > other.width:
            return 1
        if self.nchan < other.nchan:
            return -1
        if self.nchan > other.nchan:
            return 1
        return 0


    def clone (self):
        other = Window ()
        other.cfreq = self.cfreq
        other.width = self.width
        other.nchan = self.nchan
        return other


    def __getstate__ (self):
        return (self.cfreq, self.width, self.nchan)


    def __setstate__ (self, s):
        self.cfreq, self.width, self.nchan = s


    def __str__ (self):
        return '%.3f(%.3f)@%d' % (self.cfreq, self.width * 1e3, self.nchan)


    def formatLine (self):
        return '%.18e %.18e %d' % (self.cfreq, self.width, self.nchan)


    def parseLine (self, line):
        a = line.strip ().split ()
        self.cfreq = float (a[0])
        self.width = float (a[1])
        self.nchan = int (a[2])


    def getFreqs (self):
        sdf = self.width / self.nchan
        # sfreq - width / 2 gives us the left edge of the first channel,
        # whereas we want its middle; hence the extra "- sdf".
        sfreq = self.cfreq - 0.5 * (self.width - sdf)
        return np.arange (self.nchan) * sdf + sfreq


    def getFreqBounds (self):
        return self.cfreq - 0.5 * self.width, self.cfreq + 0.5 * self.width


    def mapFreq (self, ghz):
        # We'll map any frequency into a channel number; the channel
        # may just be negative, fractional, etc.

        sdf = self.width / self.nchan
        # See note in getFreqs re: sfreq
        sfreq = self.cfreq - 0.5 * (self.width - sdf)
        return (ghz - sfreq) / sdf


    def mapChan (self, chan):
        # As above, we'll map any "channel number" into a frequency.
        sdf = self.width / self.nchan
        sfreq = self.cfreq - 0.5 * (self.width - sdf)
        return sfreq + chan * sdf


class FreqConfig (object):
    __slots__ = ['nchan', 'sdata', 'wdata']

    WINTYPE_SPECTRAL = 0
    WINTYPE_MFS = 1
    WINTYPE_WIDE = 2

    def __init__ (self):
        self.nchan = 0
        self.sdata = []
        self.wdata = []


    def makeTracker (self, handle):
        t = handle.makeVarTracker ()
        t.track ('nchan', 'nspect', 'nwide', 'ischan',
                 'nschan', 'sfreq', 'sdf', 'wfreq', 'wwidth')
        return t


    def fill (self, handle):
        self.nchan = handle.getScalar ('nchan', 0)
        nspect = handle.getScalar ('nspect', 0)
        nwide = handle.getScalar ('nwide', 0)

        if nspect == 0:
            self.sdata = []
        else:
            ischan = np.atleast_1d (handle.getVarInt ('ischan', nspect))
            nschan = np.atleast_1d (handle.getVarInt ('nschan', nspect))
            sfreq = np.atleast_1d (handle.getVarDouble ('sfreq', nspect))
            sdf = np.atleast_1d (handle.getVarDouble ('sdf', nspect))

            # Fortran 1-based indices to 0-based.
            ischan -= 1

            self.sdata = zip (ischan, nschan, sfreq, sdf)

        if nwide == 0:
            self.wdata = []
        else:
            wfreq = np.atleast_1d (handle.getVarFloat ('wfreq', nwide))
            wwidth = np.atleast_1d (handle.getVarFloat ('wwidth', nwide))

            self.wdata = zip (wfreq, wwidth)

        return self


    def almosteq (self, other, tol=1e-4):
        if other.nchan != self.nchan:
            return False

        if len (other.sdata) != len (self.sdata):
            return False
        if len (other.wdata) != len (self.wdata):
            return False

        for e1, e2 in zip (self.sdata, other.sdata):
            i1, n1, sf1, sd1 = e1
            i2, n2, sf2, sd2 = e2

            if i1 != i2 or n1 != n2:
                return False
            if np.abs (sf1 / sf2 - 1) > tol:
                return False
            if np.abs (sd1 / sd2 - 1) > tol:
                return False

        for (f1, w2), (f2, w2) in zip (self.wdata, other.wdata):
            if np.abs (f1 / f2 - 1) > tol:
                return False
            if np.abs (w1 / w2 - 1) > tol:
                return False

        return True


    def __eq__ (self, other):
        if not isinstance (other, FreqConfig):
            return False
        return self.almosteq (other, tol=0)


    def __ne__ (self, other):
        return not self.__eq__ (other)


    def __nonzero__ (self):
        return self.nchan != 0 or len (self.wdata) > 0


    def __hash__ (self):
        h = hash (self.nchan)
        for i, n, sfreq, sdf in self.sdata:
            h ^= hash (i)
            h ^= hash (n)
            h ^= hash (sfreq)
            h ^= hash (sdf)
        for freq, width in self.wdata:
            h ^= hash (freq)
            h ^= hash (width)
        return h


    def clone (self):
        other = FreqConfig ()
        other.nchan = self.nchan
        other.sdata = list (self.sdata)
        other.wdata = list (self.wdata)
        return other


    def __getstate__ (self):
        return (self.nchan, self.sdata, self.wdata)


    def __setstate__ (self, s):
        self.nchan, self.sdata, self.wdata = s


    def formatLine (self):
        formats = lambda sdata: '%d/%d/%.18e/%.18e' % sdata
        formatw = lambda wdata: '%.18e/%.18e' % wdata

        sinfo = ';'.join (formats (x) for x in self.sdata)
        winfo = ';'.join (formatw (x) for x in self.wdata)
        return '%d %s %s' % (self.nchan, sinfo, winfo)


    def parseLine (self, line):
        a = line.strip ().split ()

        self.nchan = int (a[0])

        def parses (piece):
            b = piece.split ('/')
            ischan = int (b[0])
            nschan = int (b[1])
            sfreq = float (b[2])
            sdf = float (b[3])
            return ischan, nschan, sfreq, sdf

        def parsew (piece):
            b = piece.split ('/')
            freq = float (b[0])
            width = float (b[1])
            return freq, width

        self.sdata = [parses (x) for x in a[1].split (';')]
        self.wdata = [parsew (x) for x in a[1].split (';')]


    def __str__ (self):
        w = '/'.join ('%.3f(%.3f)' % (w[0], w[1] * 1e3) for w in self.wdata)
        s = '/'.join ('%.3f(%.3f)@%d(%d)' % (s[2], s[3] * 1e3, s[0], s[1]) for s in self.sdata)
        return '[%s;%s;%d]' % (s, w, self.nchan)


    def numSpectralWindows (self):
        return len (self.sdata)


    def numWideChannels (self):
        return len (self.wdata)


    def getFreqs (self):
        freqs = np.zeros (self.nchan, dtype=np.double)

        for (ischan, nschan, sfreq, sdf) in self.sdata:
            freqs[ischan:ischan+nschan] = np.arange (nschan) * sdf + sfreq

        return freqs


    def getFreqBounds (self):
        minmax = [None, None]

        def check (ghz):
            if minmax[0] is None:
                minmax[0] = minmax[1] = ghz
            else:
                minmax[0] = min (minmax[0], ghz)
                minmax[1] = max (minmax[1], ghz)

        for (ischan, nschan, sfreq, sdf) in self.sdata:
            # SDF may be < 0 so we can't assume that sfreq
            # is the smallest frequency present.
            check (sfreq - 0.5 * sdf)
            check (sfreq + (0.5 + nschan) * sdf)

        for (freq, width) in self.wdata:
            check (freq - 0.5 * width)
            check (freq + 0.5 * width)

        return tuple (minmax)


    def mapFreq (self, ghz):
        for (ischan, nschan, sfreq, sdf) in self.sdata:
            if ghz < sfreq or ghz > sfreq + sdf * nschan:
                continue

            return (ghz - sfreq) / sdf + ischan
        return None


    def allWinIdents (self):
        for i in xrange (len (self.sdata)):
            yield (self.WINTYPE_SPECTRAL, i)
            yield (self.WINTYPE_MFS, i)

        for i in xrange (len (self.wdata)):
            yield (self.WINTYPE_WIDE, i)


    def fundamentalWinIdents (self):
        for i in xrange (len (self.sdata)):
            yield (self.WINTYPE_SPECTRAL, i)
        for i in xrange (len (self.wdata)):
            yield (self.WINTYPE_WIDE, i)


    def windowFromIdent (self, ident):
        try:
            wtype, num = ident
        except StandardError:
            raise ValueError ()

        if wtype == self.WINTYPE_SPECTRAL:
            if num < 0 or num >= len (self.sdata):
                raise ValueError ()
            ischan, nschan, sfreq, sdf = self.sdata[num]
            return Window (sfreq + sdf * (nschan - 1) * 0.5, sdf * nschan, nschan)

        if wtype == self.WINTYPE_MFS:
            if num < 0 or num >= len (self.sdata):
                raise ValueError ()
            ischan, nschan, sfreq, sdf = self.sdata[num]
            return Window (sfreq + sdf * (nschan - 1) * 0.5, sdf * nschan, 1)

        if wtype == self.WINTYPE_WIDE:
            if num < 0 or num >= len (self.wdata):
                raise ValueError ()
            wfreq, wwidth = self.wdata[num]
            return Window (wfreq, wwidth, 1)

        raise ValueError ()


    def linespecFromIdent (self, ident):
        try:
            wtype, num = ident
        except StandardError:
            raise ValueError ()

        if wtype == self.WINTYPE_SPECTRAL:
            if num < 0 or num >= len (self.sdata):
                raise ValueError ()
            ischan, nschan, sfreq, sdf = self.sdata[num]
            return 'chan,%d,%d' % (nschan, ischan + 1)
        if wtype == self.WINTYPE_MFS:
            if num < 0 or num >= len (self.sdata):
                raise ValueError ()
            ischan, nschan, sfreq, sdf = self.sdata[num]
            return 'chan,1,%d,%d' % (ischan + 1, nschan)
        if wtype == self.WINTYPE_WIDE:
            if num < 0 or num >= len (self.wdata):
                raise ValueError ()
            return 'wide,%d' % (num)

        raise ValueError ()


class Pointing (object):
    __slots__ = ['rarad', 'decrad', 'name', 'aliases']

    def __init__ (self, rarad=None, decrad=None, name=None):
        self.rarad = rarad
        self.decrad = decrad
        self.name = name
        self.aliases = set ()


    def makeTracker (self, uvset):
        t = uvset.makeVarTracker ()
        t.track ('ra', 'dec', 'dra', 'ddec', 'source')
        return t


    def fill (self, uvset):
        ra = uvset.getVarDouble ('ra')
        dec = uvset.getVarDouble ('dec')

        dra = uvset.getScalar ('dra')
        if dra is not None:
            ddec = uvset.getScalar ('ddec', 0.)
            ra += dra / np.cos (dec)
            dec += ddec

        self.rarad = ra
        self.decrad = dec
        self.name = uvset.getVarString ('source')
        self.aliases.clear ()
        return self


    def __eq__ (self, other):
        """Not recommended for typical use due to coordinate
        jitter and potential name changes."""
        if not isinstance (other, Pointing):
            return False
        if other.rarad != self.rarad:
            return False
        if other.decrad != self.decrad:
            return False
        if other.name != self.name:
            return False
        return other.aliases == self.aliases


    def __ne__ (self, other):
        return not self.__eq__ (other)


    def __hash__ (self):
        h = hash (self.rarad) ^ hash (self.decrad) ^ hash (self.name)
        for a in self.aliases:
            h ^= hash (a)
        return h


    def __nonzero__ (self):
        return self.rarad is not None and self.decrad is not None and \
            self.name is not None


    def clone (self):
        other = Pointing ()
        other.rarad = self.rarad
        other.decrad = self.decrad
        other.name = self.name
        other.aliases = set (self.aliases)
        return other


    def __getstate__ (self):
        return (self.rarad, self.decrad, self.name, self.aliases)


    def __setstate__ (self, s):
        self.rarad, self.decrad, self.name, self.aliases = s


    def formatLine (self):
        a = ','.join (sorted (self.aliases))
        return '%.18e %.18e %s %s' % (self.rarad, self.decrad, self.name, a)


    def parseLine (self, line):
        a = line.strip ().split ()
        self.rarad = float (a[0])
        self.decrad = float (a[1])
        self.name = a[2]

        if len (a) < 4:
            self.aliases = set ()
        else:
            self.aliases = set (a for a in a[3].split (','))


    def sep2 (self, other, *args):
        """Return the squared separation, in squared radians, between
        this pointing and another pointing, using the small-angle
        approximation as centered on this pointing.

        If one argument is supplied, it is assumed to be another
        Pointing object. If two are supplied, they are assumed to be
        an RA and a dec, measured in radians.
        """

        if len (args) == 0:
            rarad, decrad = other.rarad, other.decrad
        elif len (args) == 1:
            rarad, decrad = other, args[0]
        else:
            raise ValueError ('expect exactly one or two arguments to sep2')

        return (self.decrad - decrad)**2 + \
            ((self.rarad - rarad) * np.cos (self.decrad))**2

# -*- mode: python; coding: utf-8 -*-
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""scanalyzer.mtutil - miriad-python utilities

This file is copied from miriad-python's mirtask.util, preserving only the
bits needed for the scanalyzer.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

# Polarizations. From subs/uvdat.h

POL_II = 0
POL_I = 1
POL_Q = 2
POL_U = 3
POL_V = 4
POL_RR = -1
POL_LL = -2
POL_RL = -3
POL_LR = -4
POL_XX = -5
POL_YY = -6
POL_XY = -7
POL_YX = -8
POL_QQ = 5
POL_UU = 6

_polNames = { POL_II: 'II', POL_I: 'I', POL_Q: 'Q',
              POL_U: 'U', POL_V: 'V', POL_RR: 'RR',
              POL_LL: 'LL', POL_RL: 'RL', POL_LR: 'LR',
              POL_XX: 'XX', POL_YY: 'YY', POL_XY: 'XY',
              POL_YX: 'YX', POL_QQ: 'QQ', POL_UU: 'UU' }

def polarizationName (polnum):
    """Return the textual description of a MIRIAD polarization type
    from its number."""

    return _polNames[polnum]

def polarizationNumber (polname):
    for (num, name) in _polNames.iteritems ():
        if name.lower () == polname.lower (): return num

    raise Exception ('Unknown polarization name \'%s\'' % polname)


# And, merging them together: antpol and basepol handling.
#
# In the following, "M" stands for the MIRIAD antenna number
# of an antenna. These numbers are 1-based. "P" stands for
# a MIRIAD polarization number, values given above.
#
# First, a "feed polarization" (f-pol) is a polarization that an
# individual feed can respond to. We include Stokes parameters here,
# even though such feeds can't operate physically, to allow sensible
# roundtripping with MIRIAD/FITS polarization values in the code
# below.

FPOL_X = 0
FPOL_Y = 1
FPOL_R = 2
FPOL_L = 3
FPOL_I = 4
FPOL_Q = 5
FPOL_U = 6
FPOL_V = 7

fPolNames = 'XYRLIQUV'

# This table helps split a MIRIAD/FITS pol code into a pair of f-pol
# values.  The pair is packed into 8 bits, the upper 3 being for the
# left pol and the lower 4 being for the right. An offset of 8 is
# required because the pol codes range from -8 to +6

_polToFPol = [0x10, 0x01, 0x11, 0x00, # YX XY YY XX
              0x32, 0x23, 0x33, 0x22, # LR RL LL RR
              0x44, # II
              0x44, 0x55, 0x66, 0x77, # I Q U V
              0x55, 0x66] # QQ UU

# This table performs the reverse mapping, with index being the two
# f-pol values packed into four bits each. A value of 0xFF indicates
# an illegal pairing. Correlations written in Stokes space are
# indicated with the single-letter FITS codes; the "II", "QQ", and
# "UU" codes are only used during pol conversion inside UVDAT.

_fpolToPol = np.ndarray (128, dtype=np.int8)
_fpolToPol.fill (0xFF)
_fpolToPol[0x00] = POL_XX
_fpolToPol[0x01] = POL_XY
_fpolToPol[0x10] = POL_YX
_fpolToPol[0x11] = POL_YY
_fpolToPol[0x22] = POL_RR
_fpolToPol[0x23] = POL_RL
_fpolToPol[0x32] = POL_LR
_fpolToPol[0x33] = POL_LL
_fpolToPol[0x44] = POL_I
_fpolToPol[0x55] = POL_Q
_fpolToPol[0x66] = POL_U
_fpolToPol[0x77] = POL_V

# A "antpol" (AP) is an integer identifying an antenna/f-pol pair. It
# can be decoded without any external information.  The translation
# between AP and M,FP is:
#
#   AP = (M - 1) << 3 + FP
#
# or
#
#   M = AP >> 3 + 1
#   P = AP & 0x7
#
# Note that arbitrarily-large antenna numbers can be encoded
# if sufficiently many bits are used to store the AP. Also note that
# if you think of the antpol "antenna number" as AP >> 3, antpol
# antenna numbers start at zero, while MIRIAD antenna numbers start at
# one.

def fmtAP (ap):
    m = (ap >> 3) + 1
    fp = ap & 0x7
    return '%d%c' % (m, fPolNames[fp])

def apAnt (ap):
    return (ap >> 3) + 1

def apFPol (ap):
    return ap & 0x7

def antpol2ap (m, fpol):
    return ((m - 1) << 3) + fpol

def parseAP (text):
    try:
        polcode = text[-1].upper ()
        fpol = fPolNames.find (polcode)
        assert fpol >= 0

        m = int (text[:-1])
        assert m > 0
    except:
        raise ValueError ('text does not specify an antpol: ' + text)

    return antpol2ap (m, fpol)


# A "basepol" is a baseline between two antpols. It is expressed as a
# 2-tuple of antpols.

def fmtBP (bp):
    ap1, ap2 = bp

    if ap1 < 0:
        raise ValueError ('first antpol %d is negative' % ap1)
    if ap2 < 0:
        raise ValueError ('second antpol %d is negative' % ap2)

    m1 = (ap1 >> 3) + 1
    fp1 = ap1 & 0x7
    m2 = (ap2 >> 3) + 1
    fp2 = ap2 & 0x7

    return '%d%c-%d%c' % (m1, fPolNames[fp1], m2, fPolNames[fp2])


def bp2aap (bp):
    """Converts a basepol into a tuple of (ant1, ant2, pol)."""

    ap1, ap2 = bp

    if ap1 < 0:
        raise ValueError ('first antpol %d is negative' % ap1)
    if ap2 < 0:
        raise ValueError ('second antpol %d is negative' % ap2)

    m1 = (ap1 >> 3) + 1
    m2 = (ap2 >> 3) + 1
    pol = _fpolToPol[((ap1 & 0x7) << 4) + (ap2 & 0x7)]

    if pol == 0xFF:
        raise ValueError ('no FITS polarization code for pairing '
                          '%c-%c' % (fPolNames[ap1 & 0x7],
                                     fPolNames[ap2 & 0x7]))

    return m1, m2, pol


def aap2bp (m1, m2, pol):
    """\
Create a basepol from antenna numbers and a FITS/MIRIAD polarization
code.

:arg int m1: the first antenna number; *one*-based as used
  internally by MIRIAD, not zero based
:arg int m2: the second antenna number; also one-based
:type pol: FITS/MIRIAD polarization code
:arg pol: the polarization
:returns: the corresponding basepol
:raises: :exc:`ValueError` if *m1* or *m2* is below one, or
  *pol* is not a known polarization code.

Note that the input antenna numbers should be one-based, not
zero-based as more conventional for C and Python. (This
is consistent with :func:`bp2aap`.) *m1* need not be
smaller than *m2*, although this is the typical convention.
"""

    if m1 < 1:
        raise ValueError ('first antenna is below 1: %s' % m1)
    if m2 < 0:
        raise ValueError ('second antenna is below 1: %s' % m2)
    if pol < POL_YX or pol > POL_UU:
        raise ValueError ('illegal polarization code %s' % pol)

    fps = _polToFPol[pol + 8]
    ap1 = ((m1 - 1) << 3) + ((fps >> 4) & 0x07)
    ap2 = ((m2 - 1) << 3) + (fps & 0x07)
    return ap1, ap2


def bp2blpol (bp):
    """Converts a basepol into a tuple of (bl, pol) where
'bl' is the MIRIAD-encoded baseline number."""

    m1, m2, pol = bp2aap (bp)
    ##return encodeBaseline (m1, m2), pol
    ## only used for sorting in ui.py
    return m1 * 100000 + m2, pol


def bpIsInten (bp):
    ap1, ap2 = bp

    if ap1 < 0:
        raise ValueError ('first antpol %d is negative' % ap1)
    if ap2 < 0:
        raise ValueError ('second antpol %d is negative' % ap2)

    fp1, fp2 = ap1 & 0x7, ap2 & 0x7
    return (fp1 >= 0 and fp1 < 5 and fp2 == fp1)


def parseBP (text):
    t1, t2 = text.split ('-', 1)

    try:
        fp1 = fPolNames.find (t1[-1].upper ())
        assert fp1 >= 0

        m1 = int (t1[:-1])
        assert m1 > 0

        fp2 = fPolNames.find (t2[-1].upper ())
        assert fp2 >= 0

        m2 = int (t2[:-1])
        assert m2 > 0
    except Exception:
        raise ValueError ('text does not specify a basepol: ' + text)

    return ((m1 - 1) << 3) + fp1, ((m2 - 1) << 3) + fp2


# A "packed 32-bit basepol" (PBP32) encodes a basepol in a single
# 32-bit integer. It can be decoded without any external
# information. The translation between PBP32 and M1,M2,FP1,FP2 is:
#
#  PBP32 = ((M1 - 1) << 19) + (FP1 << 16) + ((M2 - 1) << 3) + FP2
#
# or
#
#  M1 = (PBP32 >> 19) + 1
#  FP1 = (PBP32 >> 16) & 0x7
#  M2 = (PBP32 >> 3 & 0x1FFF) + 1
#  FP2 = PBP32 & 0x7
#
# This encoding allocates 13 bits for antenna number, which gets us up
# to 8192 antennas. This should be sufficient for most applications.

def fmtPBP32 (pbp32):
    if pbp32 < 0 or pbp32 > 0xFFFFFFFF:
        raise ValueError ('illegal PBP32 0x%x' % pbp32)

    m1 = ((pbp32 >> 19) & 0x1FFF) + 1
    fp1 = (pbp32 >> 16) & 0x7
    m2 = ((pbp32 >> 3) & 0x1FFF) + 1
    fp2 = pbp32 & 0x7

    return '%d%c-%d%c' % (m1, fPolNames[fp1], m2, fPolNames[fp2])


def pbp32ToBP (pbp32):
    if pbp32 < 0 or pbp32 > 0xFFFFFFFF:
        raise ValueError ('illegal PBP32 %x' % pbp32)
    return ((pbp32 >> 16) & 0xFFFF, pbp32 & 0xFFFF)


def bpToPBP32 (bp):
    ap1, ap2 = bp

    if ap1 < 0:
        raise ValueError ('first antpol %d is negative' % ap1)
    if ap2 < 0:
        raise ValueError ('second antpol %d is negative' % ap2)
    if ap1 > 0xFFFF:
        raise ValueError ('cannot store first antpol 0x%x in PBP32: '
                          'a1 > 0xFFFF' % ap1)
    if ap2 > 0xFFFF:
        raise ValueError ('cannot store second antpol 0x%x in PBP32: '
                          'a2 > 0xFFFF' % ap2)

    return (ap1 << 16) + (ap2 & 0xFFFF)


def pbp32IsInten (pbp32):
    if pbp32 < 0 or pbp32 > 0xFFFFFFFF:
        raise ValueError ('illegal PBP32 %x' % pbp32)

    return (pbp32 & 0x70007) in (0, 0x10001, 0x20002, 0x30003, 0x40004)


def parsePBP32 (text):
    t1, t2 = text.split ('-', 1)

    try:
        fp1 = fPolNames.find (t1[-1].upper ())
        assert fp1 >= 0

        m1 = int (t1[:-1])
        assert m1 > 0

        fp2 = fPolNames.find (t2[-1].upper ())
        assert fp2 >= 0

        m2 = int (t2[:-1])
        assert m2 > 0
    except Exception:
        raise ValueError ('text does not encode a basepol: ' + text)

    if m1 > 0x2000 or m2 > 0x2000:
        raise ValueError ('basepol cannot be encoded in a PBP32: ' + text)

    return ((m1 - 1) << 19) + (fp1 << 16) + ((m2 - 1) << 3) + fp2

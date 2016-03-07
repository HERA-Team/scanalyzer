# -*- mode: python; coding: utf-8 -*-
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""scanalyzer.transpose - interactive analysis of interferometric visbilities.

This uses the Gtk+3 graphical toolkit. This is inherited from some old
MIRIAD-targeted code so it's pretty messy.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = (b'mir_transpose ms_transpose TransposeFile').split ()

import sys, os.path, time
from struct import Struct
import numpy as np

from pwkit import binary_type
from pwkit.cli import die, warn
from . import mtutil
from .tdata import GridAxis, VisGrid, TransposeData


CACHE_SIZE = 256 * 1024 * 1024

def _sfmt (t):
    # Format timespan (t in seconds)
    if t < 90:
        return '%.0fs' % t
    if t < 4000:
        return '%.1fm' % (t / 60)
    if t < 90000:
        return '%.1fh' % (t / 3600)
    return '%.1fd' % (t / 86400)

"""
Header:
  byte order marker
  version
  nbp
  ntime
  nchan
  center frequency of first channel [GHz]
  channel width [GHz]
  midpoint time of first integration [Julian Date]
  dump cadence [Julian Date]
  correlation/flag/uvw/weight data offset
  variable data offset
"""

BYTE_ORDER_MARKER = 0x01020304
FORMAT_VERSION = 0x4
header = Struct ('=iiiiiddddQQ')

"""
Variables are just numpy arrays stored in the dataset. Below we
list the names, types, dimensions, physical units, and purposes
of standard variables.

basepols [miriad-python 32-bit packed basepol] [int: nbp]
  The basepols in the dataset for which we have data.

delays [ns] [float: (nants, 2)]
  Fixed delays for each ant/pol, used to topocentrize phases.

flaggedbps [mirpy 32-bit packed basepol] [int: nflaggedbp]
  Basepols for which all of their data are flagged.

lsts [radians] [double: nt]
  The Local Sidereal Time at each timedump.

meanuvws [ns] [double: (nbp, 3)]
  The mean UVW coordinate of each baseline, averaged over
  frequencies and times ignoring flagging. Convertible to
  wavelengths by multiplying by the frequency in GHz.

vispath [bytes] [byte: *]
  The path name of the visibility dataset used to create the
  transpose file, encoded in whatever system was used at the
  time of transposition.

transargs [bytes] [byte: *]
  A string representation of the "transpose arguments" used
  to select subsets of the input dataset. Encoded in
  whatever system was used at the time of transposition.
  Informational only and not meant to be parsable to
  recover the arguments in any reliable way.

antnames [null-padded ASCII] [byte: (nants, maxantnamelen)]
  Symbolic names of each of the antennas. This is currently
  only extracted from MeasurementSets although ATA datasets
  contain this information too. This information isn't
  actually used yet.
"""

# magic, variable name
VARIABLE_MAGIC = 0x11223344
variable = Struct ('=i32s')


# XXX MIRIAD code stale and likely busted!

def mir_transpose (vpath, tpath, transpose_args):
    try:
        return _mir_transpose (vpath, tpath, transpose_args)
    except:
        # The unlink can overwrite the exception info that
        # an argumentless "raise" would reraise.
        t, v, tb = sys.exc_info ()
        try:
            os.unlink (tpath)
        except:
            pass
        raise t, v, tb


def savevariable (stream, name, array):
    if len (name) > 32:
        raise ValueError (name)
    stream.write (variable.pack (VARIABLE_MAGIC, binary_type (name)))
    np.save (stream, array)


def _mir_transpose (vpath, tpath, unused_transpose_args):
    vis = UVData (vpath) # XXX approximate

    # Pass 1: build up list of basepols, times

    first = True
    nrecs = 0
    delays = None
    window = None
    fc = visobjs.FreqConfig ()
    times = set ()
    pbps = set ()
    visgen = vis.readLowlevel ('3', False)

    print ('pass 1 ...')

    for inp, pream, data, flags in visgen:
        t = pream[3]
        pbp = mir2pbp32 (inp, pream)
        nrecs += 1

        if first:
            ftrack = fc.makeTracker (inp)
            first = False

        if ftrack.updated ():
            fc.fill (inp)

            if fc.numSpectralWindows () != 1:
                die ('cannot transpose: need exactly one spectral window')

            idents = list (fc.fundamentalWinIdents ())
            newwindow = fc.windowFromIdent (idents[0])

            if window is not None and newwindow != window:
                die ('cannot transpose: frequency config changes inside dataset')

            window = newwindow

        if delays is None:
            nants = inp.getVarInt ('nants')
            dinfo = inp.probeVar ('delay0')

            if dinfo is None:
                delays = False
            elif dinfo[1] == 2 * nants:
                # An ATA extension: one fixed delay per antpol. Reshape
                # to be a bit more sensible
                delays = inp.getVarFloat ('delay0', 2 * nants)
                delays = delays.reshape ((2, nants)).T
            else:
                delays = inp.getVarFloat ('delay0', nants)
                delays = np.vstack ((delays, delays)).T

        times.add (t)
        pbps.add (pbp)

    # Get the timestamps onto a nice even grid, checking that our
    # gridding is decent.

    datatimes = np.asarray (sorted (times), dtype=np.double)
    nt = datatimes.size
    time0 = datatimes[0]
    cadence = np.median (datatimes[1:] - datatimes[:-1])
    tidxs = (datatimes - time0) / cadence
    timemap = np.empty (nt, dtype=np.int)
    nslot = int (round (tidxs[-1])) + 1
    scale = nslot * 1. / nt
    noff = 0

    for i in xrange (nt):
        timemap[i] = int (round (tidxs[i]))
        if (tidxs[i] - timemap[i]) > 0.01:
            noff += 1

    if noff > 0:
        warn ('had %d timestamps (out of %d) with poor '
              'mapping onto the grid', noff, nt)

    if scale > 1.05:
        warn ('data size increasing by factor of %.2f '
              'to get everything onto the time grid', scale)

    times = np.arange (nslot) * cadence + time0
    nt = nslot

    # Compute a few other things

    pbps = np.asarray (sorted (pbps), dtype=np.int32)
    nbp = pbps.size

    # Without the int(), nchan is a numpy.int32, the type of which
    # propagates along to various byte counts and offsets which end up
    # overflowing for sufficiently large datasets and causing
    # exceptions on negative values getting passed to the various
    # system calls used below.
    nchan = int (window.nchan)
    sdf = window.width / nchan
    sfreq = window.cfreq - 0.5 * (window.width - sdf)

    corr_bytes = 8 * nchan
    uvww_bytes = 4 * 8
    flag_bytes = nchan
    slice_bytes = (corr_bytes + flag_bytes + uvww_bytes) * nt
    dump_bytes = (corr_bytes + flag_bytes + uvww_bytes) * nbp

    nsimult = CACHE_SIZE // dump_bytes

    # Write out header info
    # Write-then-seek seems to break if buffering is used???

    data_offset = ((header.size + 7) // 8) * 8
    data_size = slice_bytes * nbp

    vars_offset = ((data_offset + data_size + 7) // 8) * 8

    f = open (tpath, 'w+', 0)
    f.truncate (vars_offset) # hint how big the file will be
    f.write (header.pack (BYTE_ORDER_MARKER,
                          FORMAT_VERSION,
                          nbp, nt, nchan,
                          sfreq, sdf,
                          time0, cadence,
                          data_offset,
                          vars_offset))

    # Pass 2: write data. Below we cast the tidx variables to ints for
    # the same reason as with nchan above.

    def corr_offset (bpidx, tidx):
        return data_offset + bpidx * slice_bytes + corr_bytes * int (tidx)
    def flag_offset (bpidx, tidx):
        return (data_offset + bpidx * slice_bytes + corr_bytes * nt +
                flag_bytes * int (tidx))
    def uvww_offset (bpidx, tidx):
        return (data_offset + bpidx * slice_bytes + (corr_bytes + flag_bytes) * nt +
                uvww_bytes * int (tidx))

    lsts = np.empty (nt, dtype=np.double)

    corrs = np.empty ((nsimult, nbp, nchan), dtype=np.complex64)
    flags = np.empty ((nsimult, nbp, nchan), dtype=np.int8)
    uvwws = np.empty ((nsimult, nbp, 4), dtype=np.double)
    seen = np.empty ((nsimult, nbp), dtype=np.bool)
    lstbuf = np.empty (nsimult, dtype=np.double)

    empty_corr = np.zeros (nchan, dtype=np.complex64)
    empty_flags = np.zeros (nchan, dtype=np.int8)
    empty_uvww = np.zeros (4, dtype=np.double)

    # Progress reporting:
    unbufout = os.fdopen (os.dup (1), 'w', 0)
    currec = 0
    tstart = time.time ()
    tlastprint = 0

    def dump (curtimes):
        nbatch = len (curtimes)

        tidxs = np.empty (nbatch, dtype=np.int)
        for time, sidx in curtimes.iteritems ():
            tidxs[sidx] = timemap[datatimes.searchsorted (time)]
            lsts[tidxs[sidx]] = lstbuf[sidx]

        info = np.empty ((nbatch, 3), dtype=np.int)
        info[:,0] = tidxs.argsort ()
        info[:,1] = tidxs[info[:,0]]
        info[0,2] = 1

        for i in xrange (1, nbatch):
            info[i,2] = (info[i,1] != info[i-1,1] + 1)

        for bpidx in xrange (nbp):
            for sidx, tidx, seek in info:
                if seek:
                    f.seek (corr_offset (bpidx, tidx))
                if seen[sidx,bpidx]:
                    f.write (corrs[sidx,bpidx])
                else:
                    f.write (empty_corr)

            for sidx, tidx, seek in info:
                if seek:
                    f.seek (flag_offset (bpidx, tidx))
                if seen[sidx,bpidx]:
                    f.write (flags[sidx,bpidx])
                else:
                    f.write (empty_flags)

            for sidx, tidx, seek in info:
                if seek:
                    f.seek (uvww_offset (bpidx, tidx))
                if seen[sidx,bpidx]:
                    f.write (uvwws[sidx,bpidx])
                else:
                    f.write (empty_uvww)

    newchunk = True
    curtimes = {}
    nrec = nvis = 0
    seenany = np.zeros (nbp, dtype=np.bool)
    meanuvw = np.zeros ((nbp, 3), dtype=np.double)
    muvwcounts = np.zeros (nbp, dtype=np.int)
    visgen = vis.readLowlevel ('3', False)

    print ('pass 2 ...')

    for inp, pream, data, recflags in visgen:
        uvw = pream[:3]
        t = pream[3]
        pbp = mir2pbp32 (inp, pream)
        weight = 1. / inp.getVariance ()

        if currec % 500 == 0 and currec:
            now = time.time ()

            if now - tlastprint > 1:
                pct = 100. * currec / nrecs
                elapsed = now - tstart
                total = 1. * elapsed * nrecs / currec
                eta = total - elapsed

                msg = '   %3.1f%% (%d/%d) elapsed %s ETA %s total %s' % \
                    (pct, currec, nrecs, _sfmt (elapsed), _sfmt (eta), _sfmt (total))
                print (msg.ljust (60) + '\r', end='', file=unbufout)
                tlastprint = now

        currec += 1

        if t not in curtimes and len (curtimes) == nsimult:
            msg = '   %3.1f%% (%d/%d) writing ...' % (pct, currec, nrecs)
            print (msg.ljust (60) + '\r', file=unbufout)
            dump (curtimes)
            newchunk = True

        if newchunk:
            curtimes = {}
            newchunk = False

        sidx = curtimes.get (t)

        if sidx is None:
            sidx = len (curtimes)
            curtimes[t] = sidx
            seen[sidx].fill (False)

        bpidx = pbps.searchsorted (pbp)

        seen[sidx,bpidx] = True
        uvwws[sidx,bpidx,:3] = uvw
        uvwws[sidx,bpidx,3] = weight
        corrs[sidx,bpidx] = data
        flags[sidx,bpidx] = recflags.astype (np.int8)
        lstbuf[sidx] = inp.getVarDouble ('lst')
        muvwcounts[bpidx] += 1
        meanuvw[bpidx] += uvw

        if recflags.any ():
            seenany[bpidx] = 1

        nrec += 1
        nvis += data.size

    if len (curtimes):
        msg = '   100%% (%d/%d) writing ...' % (currec, nrecs)
        print (msg.ljust (60) + '\r', end='', file=unbufout)
        dump (curtimes)

    tfinish = time.time ()
    elapsed = tfinish - tstart
    print ('   100%% (%d/%d) elapsed %s ETA 0s total %s   ' % \
           (currec, nrecs, _sfmt (elapsed), _sfmt (elapsed)))
    unbufout.close ()

    # Finally, write out variables

    f.seek (vars_offset)
    savevariable (f, 'vispath', np.fromstring (str (vis), dtype=np.byte))
    savevariable (f, 'basepols', pbps)
    if delays is not False:
        savevariable (f, 'delays', delays)
    flaggedbps = pbps[np.where (seenany == 0)]
    savevariable (f, 'flaggedbps', flaggedbps)
    savevariable (f, 'lsts', lsts)

    wbad = np.where (muvwcounts == 0)
    muvwcounts[wbad] = 1
    meanuvw[:,0] /= muvwcounts # apparently broadcasting doesn't
    meanuvw[:,1] /= muvwcounts # do what you'd want here. Not sure
    meanuvw[:,2] /= muvwcounts # why, but it's only two extra lines.
    meanuvw[wbad] = 0
    # Take the mean across the spectral window, as well as in time:
    meanuvw *= window.cfreq / sfreq
    savevariable (f, 'meanuvws', meanuvw)

    f.close ()
    return nrec, nvis, data_size


# MeasurementSet transposition

def ms_transpose (vpath, tpath, transpose_args):
    try:
        return _ms_transpose (vpath, tpath, transpose_args)
    except:
        # The unlink can overwrite the exception info that
        # an argumentless "raise" would reraise.
        t, v, tb = sys.exc_info ()
        try:
            os.unlink (tpath)
        except:
            pass
        raise t, v, tb


from pwkit.environments.casa import util as casautil
b = casautil.sanitize_unicode

def _ms_transpose (vpath, tpath, transpose_args):
    def vispath (*args):
        return b(os.path.join (vpath, *args))

    # TODO: I think that with ms.nrow() and ms.range() we can do this
    # while taking only one pass through the data.

    tb = casautil.tools.table ()
    ms = casautil.tools.ms ()
    print ('pass 1 ...')

    # Load polarization stuff we need

    tb.open (vispath ('DATA_DESCRIPTION'))
    ddid_to_pid = tb.getcol (b'POLARIZATION_ID')
    ddid_to_spwid = tb.getcol (b'SPECTRAL_WINDOW_ID')
    tb.close ()

    tb.open (vispath ('POLARIZATION'))
    numcorrs = tb.getcol (b'NUM_CORR')
    npids = numcorrs.size
    prodinfo = [None] * npids

    for i in xrange (npids):
        corrtypes = tb.getcell (b'CORR_TYPE', i)
        prodinfo[i] = [casautil.pol_to_miriad[c] for c in corrtypes]

    tb.close ()

    ddprods = [prodinfo[p] for p in ddid_to_pid]

    # Load spw configuration stuff we need. Don't grid the info yet
    # since many of the spws may be filtered out by the selection
    # setup.

    tb.open (vispath ('SPECTRAL_WINDOW'))
    nspws = tb.getcol (b'NUM_CHAN').size
    sfreqs = []

    for i in xrange (nspws):
        sfreqs.append (tb.getcell (b'CHAN_FREQ', i) * 1e-9) # Hz -> GHz

    tb.close ()

    # Antenna info

    tb.open (vispath ('ANTENNA'))
    nants = tb.getcol (b'DISH_DIAMETER').size
    names = tb.getcol (b'NAME')
    stations = tb.getcol (b'STATION')
    fullnames = []
    maxnamelen = 0

    for i in xrange (nants):
        f = '%s@%s' % (names[i], stations[i])
        fullnames.append (f)
        maxnamelen = max (maxnamelen, len (f))

    antnames = np.zeros ((nants, maxnamelen), dtype=np.byte)

    for i in xrange (nants):
        f = fullnames[i]
        n = len (f)
        antnames[i,:n] = np.fromstring (f, dtype=np.byte)

    # Open and set up filtering. msselect() says it supports
    # 'polarization' as a field, but it doesn't seem to do anything?

    ms.open (vispath ())
    ms_selectors = frozenset ('array baseline field observation polarization '
                              'scan scanintent spw taql time uvdist'.split ())
    mssel = dict (kv for kv in transpose_args.iteritems ()
                  if kv[0] in ms_selectors)
    # ms.selectinit () needed for selectpolarization() below
    ms.msselect (b(mssel))

    # Changes shape of 'data' column below. Disable for now since
    # I don't feel like debugging it.
    if 'polarization' in transpose_args:
        warn ('polarization selection not implemented for MS data')
        pass #ms.selectpolarization (transpose_args['polarization'].split (','))

    # Get table of times and basepols

    ms.iterinit (maxrows=65536) # XXX semi-arbitrary constant
    ms.iterorigin ()
    colnames = b('time antenna1 antenna2 data_desc_id'.split ())
    nrecs = 0
    times = set ()
    pbps = set ()
    seenspws = set ()

    while True:
        cols = ms.getdata (items=colnames)
        # time is (chunksize)

        for i in xrange (cols['time'].size):
            t = cols['time'][i] / 86400. + 2400000.5 # CASA to miriad timesystems

            ddid = cols['data_desc_id'][i]

            pi = ddprods[ddid]
            a1 = cols['antenna1'][i] + 1 # 0-based -> 1-based
            a2 = cols['antenna2'][i] + 1

            seenspws.add (ddid_to_spwid[ddid])

            for j in xrange (len (pi)):
                nrecs += 1
                pbp = mtutil.bpToPBP32 (mtutil.aap2bp (a1, a2, pi[j]))
                times.add (t)
                pbps.add (pbp)

        if not ms.iternext ():
            break

    # Get the timestamps onto a nice even grid, checking that our
    # gridding is decent.

    datatimes = np.asarray (sorted (times), dtype=np.double)
    nt = datatimes.size
    time0 = datatimes[0]
    cadence = np.median (datatimes[1:] - datatimes[:-1])
    tidxs = (datatimes - time0) / cadence
    timemap = np.empty (nt, dtype=np.int)
    ntslot = int (round (tidxs[-1])) + 1
    tscale = ntslot * 1. / nt
    ntoff = 0

    for i in xrange (nt):
        timemap[i] = int (round (tidxs[i]))
        if (tidxs[i] - timemap[i]) > 0.01:
            ntoff += 1

    if ntoff > 0:
        warn ('had %d timestamps (out of %d) with poor mapping onto the grid',
              ntoff, nt)

    if tscale > 1.05:
        warn ('data size increasing by factor of %.2f to get everything onto '
              'the time grid', tscale)

    times = np.arange (ntslot) * cadence + time0
    nt = ntslot

    # Now do the same thing for the spectral windows that are actually used,
    # computing lookup info for fast mapping of DDID to our frequency grid.

    freqs = set ()

    for spwid in seenspws:
        freqs.update (sfreqs[spwid])

    datafreqs = np.asarray (sorted (freqs), dtype=np.double)
    nf = datafreqs.size
    freq0 = datafreqs[0]
    sdf = np.median (datafreqs[1:] - datafreqs[:-1])
    nfslot = int (round ((datafreqs[-1] - freq0) / sdf)) + 1
    fscale = nfslot * 1. / nf
    ddfreqmap = []
    nfoff = 0
    maxnchan = 0

    for i in xrange (len (ddid_to_spwid)):
        spwid = ddid_to_spwid[i]
        if spwid not in seenspws:
            ddfreqmap.append (None)
            continue

        # If more than one DDID shares a SPWID, we're recomputing this stuff.
        # Oh well.

        ddfreqs = sfreqs[spwid]
        ddidx0 = None
        ddprevidx = None

        if ddfreqs.size > 1 and ddfreqs[1] < ddfreqs[0]:
            ddstep = -1
        else:
            ddstep = 1

        for j in xrange (ddfreqs.size):
            trueidx = (ddfreqs[j] - freq0) / sdf
            ddidx = int (round (trueidx))

            if (ddidx - trueidx) > 0.01:
                nfoff += 1

            if j == 0:
                ddidx0 = ddidx
            elif ddidx != ddprevidx + ddstep:
                die ('cannot transpose: spw must map directly onto freq grid '
                     '(spw #%d, chan %d->%d, %d->%d)', spwid, j - 1, j,
                     ddprevidx, ddidx)

            ddprevidx = ddidx

        ddfreqmap.append ((ddidx0, ddfreqs.size, ddstep))
        maxnchan = max (maxnchan, ddfreqs.size)

    if nfoff > 0:
        warn ('had %d frequencies (out of %d) with poor mapping onto the grid',
              nfoff, nf)

    if fscale > 1.05:
        warn ('data size increasing by factor of %.2f to get everything onto '
              'the frequency grid', fscale)

    freqs = np.arange (nfslot) * sdf + freq0
    nf = nfslot

    # Compute offsets and record sizes for our output file, and write
    # the header. Write-then-seek seems to break if buffering is used???

    pbps = np.asarray (sorted (pbps), dtype=np.int32)
    nbp = pbps.size

    corr_bytes = 8 * nf
    uvww_bytes = 4 * 8
    flag_bytes = nf
    slice_bytes = (corr_bytes + flag_bytes + uvww_bytes) * nt

    data_offset = ((header.size + 7) // 8) * 8
    data_size = slice_bytes * nbp

    vars_offset = ((data_offset + data_size + 7) // 8) * 8

    def corr_offset (bpidx, tidx, fidx):
        return data_offset + bpidx * slice_bytes + corr_bytes * tidx + 8 * fidx

    def flag_offset (bpidx, tidx, fidx):
        return (data_offset + bpidx * slice_bytes + corr_bytes * nt +
                flag_bytes * tidx + fidx)

    def uvww_offset (bpidx, tidx):
        return (data_offset + bpidx * slice_bytes + (corr_bytes + flag_bytes) * nt +
                uvww_bytes * tidx)

    f = open (tpath, 'w+', 0)
    f.truncate (vars_offset) # hint how big the file will be
    f.write (header.pack (BYTE_ORDER_MARKER,
                          FORMAT_VERSION,
                          nbp, nt, nf,
                          freq0, sdf,
                          time0, cadence,
                          data_offset,
                          vars_offset))

    # Our little system for buffering/writing data. Given how the CASA Python
    # interface works, I don't think we can preallocate a huge buffer that
    # everything gets stuffed in. Which is sad. TODO: smarter data structure
    # that sorts the keys as we insert them.

    buffer_size = [0] # hack so we can modify value in the funcs below
    buffer_info = {}
    buffer_data = np.empty (CACHE_SIZE, dtype=np.uint8)
    currec = 0

    def dump ():
        if not len (buffer_info):
            return

        pct = 100. * currec / nrecs
        msg = '   %3.1f%% (%d/%d) writing ...' % (pct, currec, nrecs)
        print (msg.ljust (60) + '\r', end='', file=unbufout)

        offsets = sorted (buffer_info.iterkeys ())
        curofs = None

        for offset in offsets:
            bofs, blen = buffer_info[offset]

            if curofs is None or offset != curofs:
                f.seek (offset)
            f.write (buffer_data[bofs:bofs+blen])

            curofs = offset + blen

        buffer_size[0] = 0
        buffer_info.clear ()

    def bufferview (offset, dtype, nelem):
        bofs = (buffer_size[0] + 7) & (~7) # align for safety
        blen = dtype ().nbytes * nelem

        if bofs + blen > CACHE_SIZE:
            dump ()
            bofs = 0

        # if paranoid, check that offset not already in buffer_data
        buffer_size[0] = bofs + blen
        buffer_info[offset] = (bofs, blen)
        return buffer_data[bofs:bofs+blen].view (dtype)

    # Pass 2: write data. Set up some stuff for progress reporting.
    # NOTE: we're going to keep on rewriting uvw once for each spw

    print ('pass 2 ...')

    unbufout = os.fdopen (os.dup (1), 'w', 0)
    tstart = time.time ()
    tlastprint = 0
    nvis = 0
    seenany = np.zeros (nbp, dtype=np.bool)
    meanuvw = np.zeros ((nbp, 3), dtype=np.double)
    muvwcounts = np.zeros (nbp, dtype=np.int)

    datacol = transpose_args.get ('datacol', 'data')
    colnames = b([datacol] +
                 'time antenna1 antenna2 data_desc_id flag uvw sigma'.split ())
    maxrows = CACHE_SIZE // (2 * maxnchan * 16) # 128 bits per viz.; factor of 2 safety margin
    ms.iterinit (maxrows=maxrows)
    ms.iterorigin ()

    while True:
        cols = ms.getdata (items=colnames)
        # flag and data are (npol, nchan, chunksize)
        # uvw is (3, chunksize)
        # sigma is (npol, chunksize)
        # rest are scalars, shape (chunksize)
        # data is complex128!!! converting is super slow and sad :-(

        data = cols[datacol]
        flags = cols['flag']

        for i in xrange (cols['time'].size):
            t = cols['time'][i] / 86400. + 2400000.5 # CASA to miriad timesystems
            tidx = timemap[datatimes.searchsorted (t)]
            ddid = cols['data_desc_id'][i]
            pi = ddprods[ddid]
            npol = len (pi)
            a1 = cols['antenna1'][i] + 1 # 0-based -> 1-based
            a2 = cols['antenna2'][i] + 1
            freqidx0, nchan, step = ddfreqmap[ddid]

            if currec % 100 == 0 and currec:
                now = time.time ()

                if now - tlastprint > 1:
                    pct = 100. * currec / nrecs
                    elapsed = now - tstart
                    total = 1. * elapsed * nrecs / currec
                    eta = total - elapsed

                    msg = '   %3.1f%% (%d/%d) elapsed %s ETA %s total %s' % \
                        (pct, currec, nrecs, _sfmt (elapsed), _sfmt (eta), _sfmt (total))
                    print (msg.ljust (60) + '\r', end='', file=unbufout)
                    tlastprint = now

            nvis += npol * nchan

            for j in xrange (npol):
                currec += 1
                pbp = mtutil.bpToPBP32 (mtutil.aap2bp (a1, a2, pi[j]))
                bpidx = pbps.searchsorted (pbp)

                uvww = bufferview (uvww_offset (bpidx, tidx), np.double, 4)
                uvww[:3] = cols['uvw'][:,i] * casautil.INVERSE_C_MNS
                uvww[3] = cols['sigma'][j,i]**-2
                muvwcounts[bpidx] += 1
                meanuvw[bpidx] += uvww[:3]

                corrdata = bufferview (corr_offset (bpidx, tidx, freqidx0),
                                       np.complex64, nchan)
                corrdata[:] = data[j,::step,i] # copy and convert

                flagdata = bufferview (flag_offset (bpidx, tidx, freqidx0),
                                       np.uint8, nchan)
                np.logical_not (flags[j,::step,i], flagdata)

                if flagdata.any ():
                    seenany[bpidx] = 1

        if not ms.iternext ():
            break

    dump ()

    tfinish = time.time ()
    elapsed = tfinish - tstart
    print ('   100%% (%d/%d) elapsed %s ETA 0s total %s   ' %
           (currec, nrecs, _sfmt (elapsed), _sfmt (elapsed)))
    unbufout.close ()

    # Finally, write out variables

    f.seek (vars_offset)
    savevariable (f, 'vispath', np.fromstring (b(vpath), dtype=np.byte))
    savevariable (f, 'basepols', pbps)
    savevariable (f, 'antnames', antnames)
    flaggedbps = pbps[np.where (seenany == 0)]
    savevariable (f, 'flaggedbps', flaggedbps)
    s = ' '.join ('%s=%s' % t for t in transpose_args.iteritems ())
    savevariable (f, 'transargs', np.fromstring (b(s), dtype=np.byte))

    wbad = np.where (muvwcounts == 0)
    muvwcounts[wbad] = 1
    meanuvw[:,0] /= muvwcounts # see _mir_transpose ()
    meanuvw[:,1] /= muvwcounts
    meanuvw[:,2] /= muvwcounts
    meanuvw[wbad] = 0
    meanuvw *= (freq0 + 0.5 * sdf * nf) / freq0
    savevariable (f, 'meanuvws', meanuvw)

    f.close ()
    ms.close ()
    return currec, nvis, data_size


# Data access

class TransposeFile (TransposeData):
    _overviewmode = False

    def __init__ (self, handle):
        self.handle = handle

        handle.seek (0)
        data = handle.read (header.size)
        (bom, version, nbp, nt, nchan, sfreq, sdf,
         time0, cadence, dofs, vofs) = header.unpack (data)

        if bom != BYTE_ORDER_MARKER:
            raise Exception ('wrong byte order')
        if version != FORMAT_VERSION:
            raise Exception ('unhandled format version')

        self.nbp = nbp
        self.nt = nt
        self.nchan = nchan
        self.data_offset = dofs
        self.slice_size = (8 * nchan + nchan + 4 * 8) * nt
        self.taxis = GridAxis (time0, cadence, nt)
        self.faxis = GridAxis (sfreq, sdf, nchan)

        self.vars = {}
        handle.seek (vofs)

        while True:
            data = handle.read (variable.size)
            if len (data) == 0:
                break

            magic, name = variable.unpack (data)
            if magic != VARIABLE_MAGIC:
                raise Exception ('corrupted variable table?')

            name = name.replace ('\0', '')
            self.vars[name] = np.load (handle)

        if 'basepols' not in self.vars:
            raise Exception ('no basepols var')
        self.pbps = self.vars['basepols']
        if self.pbps.size != self.nbp:
            raise Exception ('wrong-sized basepols var')

        self._curpbp = self._curgrid = None

        self._features = {}

        if 'delays' in self.vars:
            self._features['delays'] = True
        if 'meanuvws' in self.vars:
            self._features['meanuvw'] = True
        if 'flaggedbps' in self.vars:
            self._features['checkallflagged'] = True
            self._flaggedbpset = set (mtutil.pbp32ToBP (pbp) for pbp in self.vars['flaggedbps'])
        if 'antnames' in self.vars:
            self._features['antnames'] = True
            binnames = self.vars['antnames']
            nants = binnames.shape[0]
            self._antnames = antnames = [None] * nants

            for i in xrange (nants):
                antnames[i] = binnames[i].tostring ()


    def hasFeature (self, name):
        return self._features.get (name, False)


    def enableOverviewMode (self, t1, t2, f1, f2):
        self._overviewmode = True

        self.nt = 256
        self.nchan = 256
        self.data_offset = None
        self.slice_size = None
        self.taxis = GridAxis (t1, (t2 - t1) / (self.nt - 1), self.nt)
        self.faxis = GridAxis (f1, (f2 - f1) / (self.nchan - 1), self.nchan)


    def getBPs (self):
        for pbp in self.pbps:
            yield mtutil.pbp32ToBP (pbp)


    def getGrid (self, bp):
        pbp = mtutil.bpToPBP32 (bp)
        nt, nchan = self.nt, self.nchan

        if self._curpbp == pbp:
            return self._curgrid

        # Convert to Python int() to get automatic overflow handling
        # to fix same issues seen in _transpose().
        bpidx = int (self.pbps.searchsorted (pbp))

        if bpidx >= len (self.pbps) or self.pbps[bpidx] != pbp:
            raise Exception ('requesting basepol not present in the data')

        vgrid = VisGrid ()

        if self._overviewmode:
            vgrid.alloc (self.taxis, self.faxis)
            vgrid.data.fill (1)
            vgrid.data[0,0] = 0. #almost everything bright so we can see flags
            vgrid.flags.fill (1)
            vgrid.uvw.fill (0)
            vgrid.weights.fill (1)
        else:
            vgrid.taxis = self.taxis
            vgrid.faxis = self.faxis
            self.handle.seek (self.data_offset + bpidx * self.slice_size)
            data = self.handle.read (8 * nchan * nt)
            vgrid.data = np.fromstring (data, dtype=np.complex64).reshape (nt, nchan)
            data = self.handle.read (nchan * nt)
            vgrid.flags = np.fromstring (data, dtype=np.int8).reshape (nt, nchan)
            data = self.handle.read (4 * 8 * nt)
            data = np.fromstring (data, dtype=np.double).reshape (nt, 4)
            vgrid.uvw = data[:,:3]
            vgrid.weights = data[:,3]

        self._curpbp = pbp
        self._curgrid = vgrid
        return vgrid


    def getDelay (self, bp):
        from mirtask.util import POL_XX, POL_YY, POL_XY, POL_YX
        ant1, ant2, pol = bp2aap (bp)
        delays = self.vars['delays']

        if max (ant1, ant2) > delays.shape[0]:
            warn ('not enough antennas in delays array!')
            return 0

        if pol == POL_XX:
            pidx1, pidx2 = 0, 0
        elif pol == POL_YY:
            pidx1, pidx2 = 1, 1
        elif pol == POL_XY:
            pidx1, pidx2 = 0, 1
        elif pol == POL_YX:
            pidx1, pidx2 = 1, 0
        else:
            warn ('not sure what to do with this pol for delays')
            pidx1, pidx2 = 0, 0

        return delays[ant2-1,pidx2] - delays[ant1-1,pidx1]


    def getMeanUVW (self, bp):
        bpidx = int (self.pbps.searchsorted (mtutil.bpToPBP32 (bp)))
        return self.vars['meanuvws'][bpidx]


    def knownAllFlagged (self, bp):
        return bp in self._flaggedbpset

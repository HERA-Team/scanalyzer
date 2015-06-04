# -*- mode: python; coding: utf-8 -*-
# Copyright 2015 Peter Williams <peter@newton.cx> and collaborators.
# Licensed under the MIT License.

"""pwkit.scanalyzer.flag - storing and editing flag information for UV data

We consider a limited form of flags: they are regions in the
frequency/time/basepol space.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = b'FlagAPI FlagImplementation FlagTransform'.split ()


class FlagAPI (object):
    def add (self, tst, tend, fst, fend, bpmatch):
        """Create or modify a record. Returns its ident.

        If a record with the specified time/freq bounds already exists,
        it is modified and its ident is returned. Otherwise, a new record
        is added.

        *bpmatch* is either None, indicating that the record applies to
        all basepols, or a set of antpol 2-tuples.
        """
        raise NotImplementedError ()

    def get (self, ident):
        """Return (tst, tend, fst, fend, bpmatch)."""
        raise NotImplementedError ()

    def setBPMatch (self, ident, bpmatch):
        """Change the basepol matching component of a record."""
        raise NotImplementedError ()

    def resize (self, ident, tst, tend, fst, fend):
        """Change the freq/time shape of an existing region."""
        raise NotImplementedError ()

    def delete (self, ident):
        """Delete the record."""
        raise NotImplementedError ()

    def toggleBasepol (self, ident, bp):
        """Modify the bpmatch of the specified record by toggling
        whether the specified record applies to the specified
        antpols. A noop if the record has a broader basepol selection
        that causes it to apply to *bp* anyway."""
        raise NotImplementedError ()

    def applyToRegion (self, bp, taxis, faxis, flags):
        """*bp* is a basepol
        *flags* is array of bools of shape (ntimes, nfreqs)
        """
        raise NotImplementedError ()

    def getRecords (self, taxis, faxis):
        """Get list of flag records applying to the given region
        Returns iterable of (ident, tst, tend, fst, fend).
        """
        raise NotImplementedError ()

    def checkBasepol (self, bp):
        """Return whether the specified basepol is flagged given the
        current set of records."""
        raise NotImplementedError ()


TST, TEND, FST, FEND, BPMATCH = range (5)
DELETED = set () # right semantics for deleted records

def bpmatches (bpmatch, bp):
    if bpmatch is None:
        return True
    return bp in bpmatch

def normalrecord (rec):
    bpmatch = rec[BPMATCH]

    if bpmatch is DELETED:
        return False
    if bpmatch is not None and len (bpmatch) == 0:
        return False
    return True

def cmprecord (left, right):
    # The Python implementation of __cmp__ for tuples is good enough
    # for comparing the frequency/time components. We just need to be
    # crafty about the sets.

    lft = tuple (left[:4])
    rft = tuple (right[:4])

    result = cmp (lft, rft)
    if result:
        return result

    lbpmatch = left[BPMATCH]
    rbpmatch = right[BPMATCH]

    if lbpmatch is None:
        if rbpmatch is None:
            return 0
        return -1
    elif rbpmatch is None:
        return 1

    lbplen = len (lbpmatch)
    rbplen = len (rbpmatch)

    result = cmp (lbplen, rbplen)
    if result:
        return result

    return cmp (tuple (sorted (lbpmatch)), tuple (sorted (rbpmatch)))


class FlagImplementation (FlagAPI):
    def __init__ (self, path):
        self.records = []
        self.path = path


    def add (self, tst, tend, fst, fend, bpmatch):
        if (tst is None) ^ (tend is None):
            raise ValueError ('both or neither of tst and tend must be None')
        if (fst is None) ^ (fend is None):
            raise ValueError ('both or neither of fst and fend must be None')
        if not (bpmatch is None or isinstance (bpmatch, set)):
            raise ValueError ('bpmatch must be a set or None')

        cbounds = (tst, tend, fst, fend)

        for ident, (itst, itend, ifst, ifend, ibpmatch) in enumerate (self.records):
            if ibpmatch is DELETED:
                continue
            if (itst, itend, ifst, ifend) == cbounds:
                if bpmatch is None:
                    self.records[ident][BPMATCH] = None
                elif ibpmatch is None:
                    pass # noop, already an all-matcher
                else:
                    ibpmatch.update (bpmatch)
                return ident

        # Need to create new record.
        self.records.append ([tst, tend, fst, fend, bpmatch])
        return len (self.records) - 1


    def get (self, ident):
        if self.records[ident][BPMATCH] is DELETED:
            raise ValueError ('querying deleted record')
        return self.records[ident]


    def setBPMatch (self, ident, bpmatch):
        if self.records[ident][BPMATCH] is DELETED:
            raise ValueError ('modifying deleted record')
        if not (bpmatch is None or isinstance (bpmatch, set)):
            raise ValueError ('bpmatch must be a set or None')

        self.records[ident][BPMATCH] = bpmatch


    def resize (self, ident, tst, tend, fst, fend):
        if self.records[ident][BPMATCH] is DELETED:
            raise ValueError ('modifying deleted record')
        if (tst is None) ^ (tend is None):
            raise ValueError ('both or neither of tst and tend must be None')
        if (fst is None) ^ (fend is None):
            raise ValueError ('both or neither of fst and fend must be None')

        self.records[ident][TST:BPMATCH] = tst, tend, fst, fend


    def delete (self, ident):
        if self.records[ident][BPMATCH] is DELETED:
            raise ValueError ('modifying deleted record')
        self.records[ident][BPMATCH] = DELETED


    def toggleBasepol (self, ident, bp):
        """Returns whether the record was actually modified."""

        bpm = self.records[ident][BPMATCH]

        if bpm is DELETED:
            raise ValueError ('modifying deleted record')

        if bpm is None:
            # noop, matches all antpols
            return False

        if bp in bpm:
            bpm.remove (bp)
        else:
            bpm.add (bp)

        return True


    def applyToRegion (self, bp, taxis, faxis, flags):
        assert taxis.size == flags.shape[0]
        assert faxis.size == flags.shape[1]

        for tst, tend, fst, fend, bpmatch in self.records:
            if not bpmatches (bpmatch, bp):
                continue

            if tst is None and fst is None:
                # Unbounded in frequency and time. Total flaggitude.
                flags.fill (0)
            elif tst is None:
                # Unbounded in time.
                sl = faxis.overlapSlice (fst, fend)
                if sl is None:
                    continue

                flags[:,sl] = 0
            elif fst is None:
                # Unbounded in frequency.
                sl = taxis.overlapSlice (tst, tend)
                if sl is None:
                    continue

                flags[sl,:] = 0
            else:
                # Bounded in both.
                tsl = taxis.overlapSlice (tst, tend)
                if tsl is None:
                    continue

                fsl = faxis.overlapSlice (fst, fend)
                if fsl is None:
                    continue

                flags[tsl,fsl] = 0


    def getRecords (self, taxis, faxis):
        for ident, (tst, tend, fst, fend, bpmatch) in enumerate (self.records):
            if bpmatch is DELETED:
                continue
            if tst is not None and not taxis.overlaps (tst, tend):
                continue
            if fst is not None and not faxis.overlaps (fst, fend):
                continue

            yield ident, tst, tend, fst, fend


    def checkBasepol (self, bp):
        for tst, tend, fst, fend, bpmatch in self.records:
            if bpmatch is DELETED:
                continue
            if tst is not None or fst is not None:
                # only looking for all-freq/all-time flags
                continue
            if bpmatches (bpmatch, bp):
                return True

        return False


    def _save (self, stream):
        # We can currently write out our data in ARF multiflag
        # format and still save/restore all of the information we need.
        from .mtutil import jdToFull, fmtBP

        s = sorted ((rec for rec in self.records if normalrecord (rec)),
                     cmp=cmprecord)

        for tst, tend, fst, fend, bpmatch in s:
            record = []

            if tst is not None:
                record.append ('time(%s,%s)' % (jdToFull (tst), jdToFull (tend)))
            if fst is not None:
                record.append ('freq(%.18e,%.18e)' % (fst, fend))
            if bpmatch is not None:
                record.append ('basepol(%s)' % ','.join (fmtBP (x) for x in sorted (bpmatch)))

            print >>stream, ' '.join (record)


    def _load (self, stream):
        raise NotImplementedError


    def try_load (self):
        from ..io import try_open

        f = try_open (self.path)
        if f is not None:
            self._load (f)
            f.close ()


    def save (self):
        from os import rename
        self._save (open (self.path + '.new', 'w'))
        try:
            rename (self.path, self.path + '~')
        except OSError, e:
            if e.errno != 2:
                raise
        rename (self.path + '.new', self.path)


from .tdata import DataTransform

class FlagTransform (DataTransform):
    def __init__ (self, flagapi):
        self.flagapi = flagapi


    def getGrid (self, bp):
        vgrid = self.parent.getGrid (bp).shallowcopy ()
        vgrid.flags = vgrid.flags.copy ()
        self.flagapi.applyToRegion (bp, vgrid.taxis,
                                    vgrid.faxis, vgrid.flags)
        return vgrid


    _myfeatures = frozenset (('checkallflagged', ))

    def knownAllFlagged (self, bp):
        if self.parent.hasFeature ('checkallflagged'):
            if self.parent.knownAllFlagged (bp):
                return True

        return self.flagapi.checkBasepol (bp)

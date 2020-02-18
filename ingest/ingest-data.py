#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from os.path import isfile, isdir
from os import makedirs, path
from itertools import islice

import sys
sys.path.append('../') # this is to find libStealthy
from libStealthy.util import hqrand, DataFile, printProgress, window, progressBarIterable, setArgsFromConfig
from libStealthy.plots import selectSpan, checkOverlap

import configparser
import numpy as np
import shutil
import json

def genIPDChunks(count, num, length, offset=0):
    # count -> the number of IPDs
    # num -> number of reference strings
    # length -> length of a reference string
    if count < length:
        raise IndexError('total section overlap detected')
    starts = ((i * int(count // num)) + offset for i in range(num))
    for start in starts:
        end = start + length
        yield (start, end)

def main(args):
    # Read the configuration file
    try:
        args = setArgsFromConfig(args.testbedcfg, args, 'gencovers-ingest')
        try:
            args = setArgsFromConfig(args.testbedcfg, args)
        except:
            raise
    except FileNotFoundError:
        print('[!!!] Unable to find configuration file %s, aborting.' % args.testbedcfg)
        sys.exit(1)
    except KeyError:
        print('[~~~] No file specfic configs were found.')
    args.covlens = sorted(args.covlens)

    print('[###] Log file will be', args.logfile)
    print('[###] Data hierarchy rooted at', args.dataroot)
    
    if not isdir(args.dataroot):
        raise RuntimeError('[!!!] Dataroot %s is not a directory!' % args.dataroot)

    args.dataroot = path.join(args.dataroot, 'ingest')

#     if args.overrideData == 'yes':
#         shutil.rmtree(args.dataroot)
#         makedirs(args.dataroot)

    odir_for_plain = path.join(args.dataroot, 'plain/')
    odir_for_embed = path.join(args.dataroot, 'cover/')
    odir_for_ref = path.join(args.dataroot, 'reference/')
    dirs = (odir_for_plain, odir_for_embed, odir_for_ref)
    for d in dirs:
        if not isdir(d):
            makedirs(d)

    with open(args.logfile, 'a') as logstream:
        if args.parserID == 'perc':
            IPDs = []
            for filePath in args.ingest:
                try:
                    with open(filePath) as fp:
                        IPDs.extend(json.load(fp))
                except FileNotFoundError:
                    print('Unable to locate file: %s' % filePath)

        elif args.parserID == 'range':
            if len(args.ingest) != 1:
                print('Exactly one input file is allowed to be specified for range based ingest.')
                sys.exit()
            with open(args.ingest[0]) as fp:
                IPDs = json.load(fp)

        ipdLen = len(IPDs)

        refsForCl = []
        coversForCL = []
        plainForCL = []
        for cl in args.covlens:
            starts = []
            ends = []
            if args.parserID == 'perc':
                try:
                    refIPDCount = int(ipdLen * args.ref) #### 20% of the whole string
                    startEnd = tuple(genIPDChunks(refIPDCount, args.numrefs, args.reflen))
                    starts.extend((start for start, _ in startEnd)) ### extend means this list is not overwritten
                    ends.extend((end for _, end in startEnd))
                    refIPDs = (IPDs[start:end] for start, end, in startEnd)

                    coverIPDCount = int(ipdLen * args.cover) #### 40% of the whole string
                    if refIPDCount >= ends[-1]:
                        coverStartOffset = ends[-1] ## if we don't use up all of 20%, start as soon as possible
                    else:
                        coverStartOffset = refIPDCount
                    startEnd = tuple(genIPDChunks(coverIPDCount, args.numcovs, cl, offset=coverStartOffset))
                    starts.extend((start for start, _ in startEnd))
                    ends.extend((end for _, end in startEnd))
                    coverIPDs = (IPDs[start:end] for start, end, in startEnd)

                    plainIPDCount = int(ipdLen * args.plain) #### 40% of the whole string
                    if coverIPDCount + coverStartOffset >= ends[-1]:
                        plainStartOffset = ends[-1] ### start as soon as possible
                    else:
                        plainStartOffset = coverStartOffset + coverIPDCount
                    startEnd = tuple(genIPDChunks(plainIPDCount, args.numcovs, cl, offset=plainStartOffset))
                    starts.extend((start for start, _ in startEnd))
                    ends.extend((end for _, end in startEnd))
                    plainIPDs = (IPDs[start:end] for start, end, in startEnd)
                except IndexError:
                    print('Total section overlap detected for CL: %d' % cl)
                    continue

            elif args.parserID == 'range':
                if args.interactive:
                    if not any(args.ref):
                        title = 'Press left mouse button and drag select range for ref.'
                        args.ref = selectSpan(IPDs, title=title)
                    if not any(args.cover):
                        title = 'Press left mouse button and drag select range for cover.'
                        args.cover = selectSpan(IPDs, title=title)
                    if not any(args.plain):
                        title = 'Press left mouse button and drag select range for plain.'
                        args.plain = selectSpan(IPDs, title=title)
                refIPDs = [islice(IPDs, *args.ref)]
                coverIPDs = [islice(IPDs, *args.cover)]
                plainIPDs = [islice(IPDs, *args.plain)]
                args.numrefs = 1
                args.numcovs = 1
                starts = (start for start, _ in [args.ref, args.cover, args.plain])
                ends = (end for _, end in [args.ref, args.cover, args.plain])

            else:
                print("An action must be specified. Actions available: 'perc', 'range'")
                sys.exit()

            if args.checkOverlap:
                if args.pdf:
                    pdf = path.join(args.pdf, 'covLen%dOverlap.pdf' % cl)
                    checkOverlap(IPDs, starts, ends, title='Overlap for Cover Len: %d' % cl, pdf=pdf)
                else:
                    checkOverlap(IPDs, starts, ends, title='Overlap for Cover Len: %d' % cl)

            refsForCl.append(refIPDs)
            coversForCL.append(coverIPDs)
            plainForCL.append(plainIPDs)

        for cl, refIPDs, coverIPDs, plainIPDs in zip(args.covlens, refsForCl, coversForCL, plainForCL):
            refIPDs = progressBarIterable(refIPDs, prefix='Cover Length: %d Refs ' % cl, count=args.numrefs)
            coverIPDs = progressBarIterable(coverIPDs, prefix='Cover Length: %d Embed ' % cl, count=args.numcovs)
            plainIPDs = progressBarIterable(plainIPDs, prefix='Cover Length: %d Plain ' % cl, count=args.numcovs)

            ofile_ref_prefix = 'cover-for-ref-ingest-n%d' % args.reflen
            for i, ipds in enumerate(refIPDs):
                if ipds:
                    ofile = '%s%s-%03d.json' % (odir_for_ref, ofile_ref_prefix, i)
                    # if the file exist and we are not forcing an
                    # overwite, write to the logstream and skip the
                    # rest of the cycle.
                    if not args.force and isfile(ofile):
                        logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                        continue
                    with open(ofile, 'w') as fp:
                        json.dump(ipds, fp)

            ofile_cover_prefix = 'cover-for-embed-ingest-n%d' % cl
            for i, ipds in enumerate(coverIPDs):
                if ipds:
                    ofile = '%s%s-%03d.json' % (odir_for_embed, ofile_cover_prefix, i)
                    # if the file exist and we are not forcing an
                    # overwite, write to the logstream and skip the
                    # rest of the cycle.
                    if not args.force and isfile(ofile):
                        logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                        continue
                    with open(ofile, 'w') as fp:
                        json.dump(ipds, fp)

            ofile_plain_prefix = 'cover-plain-ingest-n%d' % cl
            for i, ipds in enumerate(plainIPDs):
                if ipds:
                    ofile = '%s%s-%03d.json' % (odir_for_plain, ofile_plain_prefix, i)
                    # if the file exist and we are not forcing an
                    # overwite, write to the logstream and skip the
                    # rest of the this cycle.
                    if not args.force and isfile(ofile):
                        logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                        continue
                    with open(ofile, 'w') as fp:
                        json.dump(ipds, fp)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate cover traffic from real IPDs for reference, plain traffic, and embedding.')

    parser.add_argument('--dataroot', metavar='DIR', type=str, default=None, help='Root directory for data hierarchy')
    parser.add_argument('--logfile', type=str, default='./covers.log', help='Name of log file for running covers')
    parser.add_argument('--testbedcfg', metavar='CFGFILE', type=str, default='../nextstep.cfg', help='Name of general testbed configuration file (default is ../nextstep.cfg)')
    parser.add_argument('--numrefs', type=int, default=None)
    parser.add_argument('--numcovs', type=int, default=None)
    parser.add_argument('--covlens', nargs='*', default=None)
    parser.add_argument('--reflen', type=int, default=None)
    parser.add_argument('-c', '--checkOverlap', action='store_true', default=False, help='Allow user to see where overlaps occur.  If in a non-gui enviroment (SSH), the -p/--pdf flag should also be set')
    parser.add_argument('-f', '--force', action='store_true', default=False, help='Force overwrite of existing covers/plains/refs.')
    parser.add_argument('-p', '--pdf', type=str, default='', help='The directory you want the PDF versons of the plots to go.  This is only used for checking overlaps')

    # Keep all files in HQ or delete all files to generate new ones
#     parser.add_argument('--overrideData', default='no', type=str, required=True, choices=['no', 'yes'], help="Delete previously generated files.")


    subParsers = parser.add_subparsers()

    rangeParser = subParsers.add_parser('range', help='Use absolute ranges to split up an input file into ref, cover and plain files')
    rangeParser.add_argument('--ref', type=int, nargs=2, default=(0,0))
    rangeParser.add_argument('--cover', type=int, nargs=2, default=(0,0))
    rangeParser.add_argument('--plain', type=int, nargs=2, default=(0,0))
    rangeParser.add_argument('-i', '--interactive', action='store_true', default=False, help='Use the interactive range finder to set the ranges.  Really slow with large datasets.')
    rangeParser.add_argument('ingest', type=str, nargs=1, help='The input json blobs.')
    rangeParser.set_defaults(parserID = 'range')

    percParser = subParsers.add_parser('perc')
    percParser.add_argument('--ref', type=float, default=0.2)
    percParser.add_argument('--cover', type=float, default=0.4)
    percParser.add_argument('--plain', type=float, default=0.4)
    percParser.add_argument('ingest', type=str, nargs='+', help='The input json blobs.')
    percParser.set_defaults(parserID = 'perc')

    args = parser.parse_args()

    main(args)

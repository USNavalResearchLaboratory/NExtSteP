#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from scipy.stats import ttest_ind, mannwhitneyu
from Crypto.Random import random as HQrandom
from math import modf, log, log10
from os.path import isfile, isdir
from os import makedirs, path

import sys
sys.path.append('../') # this is to find libStealthy
from libStealthy.util import rateF, hqrand, calc_mssage_length, calc_mssage_prob, calc_rate_constants, DataFile, SeedBank, progressBarIterable, setArgsFromConfig
from libStealthy.stego.fromIPDs import SCTimeReplayProbFromIPDsUEHQ, SCTimeReplayPercentileFromIPDsUEHQ, SCTimeReplayProbFromIPDsUEHQ_alphabet
from libStealthy.pRandom import PRNG

import configparser
import numpy as np
import warnings
import json
import math
import glob
import os

embeddedRates = ['const', 'log', 'cubrt', 'sqrt', '2rt3', 'linear']

def main(args):
    try:
        args = setArgsFromConfig(args.testbedcfg, args, 'embed')
        try:
            args = setArgsFromConfig(args.testbedcfg, args)
        except:
            raise
    except FileNotFoundError:
        print('[!!!] Unable to find configuration file %s, aborting.' % args.testbedcfg)
        sys.exit(1)
    except KeyError:
        print('[~~~] No file specfic configs were found.')

# # Various possible info. messages to print
#
#     print('[###] Log file will be', args.logfile)
#     print('[###] Data hierarchy rooted at', args.dataroot)
#     print('[###] Logs will be in', args.logdir)

    args.covlens = sorted(args.covlens)
    
    fastOutput = not args.humanReadable

    if not isdir(args.dataroot):
        raise RuntimeError('[!!!] Dataroot %s is not a directory!' % args.dataroot)
    if args.logdir is None:
        args.logdir = path.join(args.dataroot, 'logs')
    if not isdir(args.logdir):
         makedirs(args.logdir) 
    if 'ingest' in args.dist:
        args.dataroot = path.join(args.dataroot, 'ingest/')
    else:
        args.dataroot = path.join(args.dataroot, 'hq/')
    if args.bankPath is None:
        args.bankPath = path.join(args.dataroot, 'embedSeedbank.json')
    

    print('[###] Using', len(args.covlens), 'cover lengths:', args.covlens)


    if args.minidx < 0:
        raise RuntimeError('[!!!] Minimum index must be non-negative')

    # FUTURE: make the location an input parameter
    msg_ifile = '../urand-b' + str(args.nbits) + '.txt'


    # This is the number of bits to embed for most methods
    # For time-replay embedding, this will give the number of IPDs to change
    k0 = args.rateK 
    N0 = args.rateIPD 
    rateC = calc_rate_constants(N0, k0)
    
    normStr = '-k' + str(k0) + '-n' + str(N0)

    # How many bits should be embedded when using (# bits changed) as numerator
    MLs = [calc_mssage_length(cl, args.rate, rateC[args.rate]) for cl in args.covlens]
    print('[###] MLs are:', list(MLs))

    PRNGClass = PRNG.getRng(args.PRNGSource)
    #PRNGClass = HQRandom.StrongRandom()

    with open(path.join(args.logdir, args.logfile), 'a') as logstream:
        if not isfile(msg_ifile):
            logstream.write('FILE: Unable to find ifile (message to embed): ' + msg_ifile + '\n')
            raise RuntimeError('[!!!] Unable to find file with message to embed!')
        with open(msg_ifile, 'r') as istream:
            line = istream.readline()

        full_msg = "".join(["{0:04b}".format(int(str(c),16)) for c in line])
        # Set up inverseCDF
        if args.dist == 'sellkepareto' or args.dist == 'paretoa100b095':
            alpha = 100
            beta = 0.95
            ifileDistParamStr = 'sellke'
            ofileDistParamStr = 'sellkepareto'
        elif args.dist == 'narrow' or args.dist == 'paretoa100b10':
            alpha = 100
            beta = 10
            ifileDistParamStr = 'pareto-a100-b10'
            ofileDistParamStr = 'pareto-a100-b10'
        elif args.dist == 'pareto':
            alpha = args.paretoa
            beta = args.paretob
            if beta <= 0.0:
                raise ValueError('Beta must be positive!')
            ifileDistParamStr = 'pareto-a%d-b%s' % (alpha, str(beta).replace('.', 'p'))
            ofileDistParamStr = ifileDistParamStr
        elif args.dist == 'normal':
            distmean = args.mean
            sd = args.sd
            if sd <= 0.0:
                raise ValueError('Std. dev. must be positive!')
            ifileDistParamStr = 'normal-m%d-s%s' % (distmean, str(sd).replace('.', 'p'))
            ofileDistParamStr = ifileDistParamStr
        elif args.dist == 'ingest':
            ifileDistParamStr = 'ingest'
            ofileDistParamStr = 'ingest'
        elif args.dist == 'markov-ingest':
            ifileDistParamStr = 'markov-ingest'
            ofileDistParamStr = 'markov-ingest'

        if args.replayperc:
            if args.replayQ == None:
                raise RuntimeError('[!!!] Must provide one or more percentiles for replayperc!')

            for q in args.replayQ:
                if q < 0 or q > 100:
                    warnings.warn('Percentile q must be between 0 and 100!')
                idir_for_cover = path.join(args.dataroot, 'cover/')
                odir_for_embedded = path.join(args.dataroot, 'embedded/replayperc-Q' + str(int(q)) + normStr + '/')

                if not isdir(idir_for_cover):
                    raise RuntimeError('[!!!] Input directory not found! %s' % idir_for_cover)                    

                if not isdir(odir_for_embedded):
                    makedirs(odir_for_embedded)

                with SeedBank(args.bankPath, args.runID) as bank:
                    for msgLength, cl in zip(MLs, args.covlens):
                        if msgLength <= args.nbits:
                            msg = full_msg[:msgLength]
                        else:
                            raise RuntimeError('[!!!] Too few message bits!')

                        ifile_for_cover_prefix = 'cover-for-embed-' + ifileDistParamStr + '-n' + str(cl) + '-'
                        ofile_for_embedded_prefix = 'embedded-replayperc-' + 'Q' + str(int(q)) + '-m' + str(msgLength) + '-' + ofileDistParamStr + '-n' + str(cl) + '-'

                        inputFiles = DataFile.getValidDataFiles(idir_for_cover, ifile_for_cover_prefix)
                        if not args.reuseSeeds:
                            ifilesWithSeeds = bank.getNewSeeds(inputFiles, subRunID=cl, seedBits=PRNGClass.seedSize)
                        else:
                            ifilesWithSeeds = bank.getSeedsByID(cl)
                        prefix = '[###] Starting embedding for cl = %d' % cl
# # Same as following with more verbose reporting
#
#                         for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds)), prefix=prefix):
                        for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds))):
                            ifile, seed = ifileWithSeed
                            try:
                                df = DataFile(ifile)
                                coverIPDs = df.getData()
                                df.addMetaData('embedderSeed', seed)
                                df.addMetaData('PRNG', args.PRNGSource)
                                embedderObj = SCTimeReplayPercentileFromIPDsUEHQ(coverIPDs, PRNGClass, q, seed)
                            except FileNotFoundError:
                                logstream.write('FILE: Unable to find ifile: ' + ifile + '\n')

                            ofile = '%s%s%03d.json' % (odir_for_embedded, ofile_for_embedded_prefix, ifile_idx )
                            if isfile(ofile):
                                logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                                continue

                            stegoIPDs = embedderObj.generateIPDseq(msg, cl)

                            df.addMetaData('embeddingMethod', 'perc')
                            df.addMetaData('msglength', msgLength)
                            df.setData(stegoIPDs)
                            df.write(ofile, fastOutput)

        if args.donothing:
            # An embedder that does nothing
            # Produce the "embedded" files in the expected place,
            # but these are the same as the cover files.        
            
            # Force the rate to be sqrt so that we don't create
            # too many copies of data
            if args.rate != 'sqrt':
                raise RuntimeError('[!!!] Do-nothing embedder called with rate other thatn \'sqrt\'!')

            idir_for_cover = path.join(args.dataroot, 'cover/')
            odir_for_embedded = path.join(args.dataroot, 'embedded/donothing/')

            if not isdir(idir_for_cover):
                raise RuntimeError('[!!!] Input directory not found! %s' % idir_for_cover)                    

            if not isdir(odir_for_embedded):
                makedirs(odir_for_embedded)

            with SeedBank(args.bankPath, args.runID) as bank:
                for msgLength, cl in zip(MLs, args.covlens):
                    if msgLength <= args.nbits:
                        msg = full_msg[:msgLength]
                    else:
                        # Now raise error if not enough message bits
                        raise RuntimeError('[!!!] Too few message bits!')

                    ifile_for_cover_prefix = 'cover-for-embed-' + ifileDistParamStr + '-n' + str(cl) + '-'
                    ofile_for_embedded_prefix = 'embedded-donothing' + '-m' + str(msgLength) + '-' + ofileDistParamStr + '-n' + str(cl) + '-'

                    inputFiles = DataFile.getValidDataFiles(idir_for_cover, ifile_for_cover_prefix)
                    if not args.reuseSeeds:
                        ifilesWithSeeds = bank.getNewSeeds(inputFiles, subRunID=cl, seedBits=PRNGClass.seedSize)
                    else:
                        ifilesWithSeeds = bank.getSeedsByID(cl)
                    prefix = '[###] Starting embedding for cl = %d' % cl
# # Same as following with more verbose reporting
#
#                     for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds)), prefix=prefix):
                    for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds))):
                        ifile, seed = ifileWithSeed
                        try:
                            df = DataFile(ifile)
                            coverIPDs = df.getData()
                            df.addMetaData('embedderSeed', seed)
                            df.addMetaData('PRNG', args.PRNGSource)
                        except FileNotFoundError:
                            logstream.write('FILE: Unable to find ifile: ' + ifile + '\n')

                        ofile = '%s%s%03d.json' % (odir_for_embedded, ofile_for_embedded_prefix, ifile_idx )
                        if isfile(ofile):
                            logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                            continue
                        df.addMetaData('embeddingMethod', 'donothing')
                        df.addMetaData('msglength', msgLength)
                        # Simply use the coverIPDs as the data to write
                        # These should already be in df, but we'll use
                        # df.setData just to make sure/to parallel the 
                        # code for replayperc
                        df.setData(coverIPDs)
                        df.write(ofile, fastOutput)

        if args.replayprob:
            if args.replayQ == None:
                raise RuntimeError('[!!!] Must provide one or more percentiles for replayprob!')
            for q in args.replayQ:
                if q < 0 or q > 100:
                    warnings.warn('Percentile q must be between 0 and 100!')
                idir_for_cover = path.join(args.dataroot, 'cover/')
                odir_for_embedded = path.join(args.dataroot + 'embedded/replayprob-Q' + str(int(q)) + normStr + '/')

                if not isdir(idir_for_cover):
                    raise RuntimeError('[!!!] Input directory not found! %s' % idir_for_cover)

                if not isdir(odir_for_embedded):
                    makedirs(odir_for_embedded)

#                 C = args.rateK / args.rateIPD * rateF[args.rate](args.rateIPD)
                with SeedBank(args.bankPath, args.runID) as bank:
                    for cl in args.covlens:
                        if args.replayP == None:
#                             prob = C / rateF[args.rate](cl)
                            prob = calc_mssage_prob(cl, args.rate, rateC[args.rate])
                        else:
                            prob = args.replayP
                        
                        msgLength = cl
                        if msgLength <= args.nbits:
                            msg = full_msg[:msgLength]
                        else:
                            # Now raise error if not enough message bits
                            raise RuntimeError('[!!!] Too few message bits!')

                        ifile_for_cover_prefix = 'cover-for-embed-' + ifileDistParamStr + '-n' + str(cl) + '-'
                        ofile_for_embedded_prefix = 'embedded-replayprob-' + 'Q' + str(int(q)) + '-' + ofileDistParamStr + '-n' + str(cl) + '-'

                        inputFiles = DataFile.getValidDataFiles(idir_for_cover, ifile_for_cover_prefix)
                        if not args.reuseSeeds:
                            ifilesWithSeeds = bank.getNewSeeds(inputFiles, subRunID=cl, seedBits=PRNGClass.seedSize)
                        else:
                            ifilesWithSeeds = bank.getSeedsByID(cl)
                        prefix = 'Starting embedding for cl = %d' % cl
# # Same as following with more verbose reporting
#
#                         for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds)), prefix=prefix):
                        for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds))):
                            ifile, seed = ifileWithSeed

                            df = DataFile(ifile)

                            df.addMetaData('embedderSeed', seed)
                            df.addMetaData('PRNG', args.PRNGSource)
                            df.addMetaData('rate', args.rate)
                            coverIPDs = df.getData()
                            embedderObj=SCTimeReplayProbFromIPDsUEHQ(coverIPDs, PRNGClass, q, seed)

                            stegoIPDs, msgbits = embedderObj.generateIPDseq(msg, prob,cl)
                            ofile = '%s%sp%d-%03d.json' % (odir_for_embedded, ofile_for_embedded_prefix, int(10**5*prob), ifile_idx)
                            if isfile(ofile):
                                logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                                continue
                            df.addMetaData('embeddingMethod', 'prob')
                            df.addMetaData('msglength', msgLength)
                            df.addMetaData('msglengths_actually_embedded', msgbits)
                            df.setData(stegoIPDs)
                            df.write(ofile, fastOutput)

        if args.replayprob_alphabet:
            is_power_2 = lambda N: ((N & (N - 1)) == 0) and N != 0
            if is_power_2(len(args.replayQ)+1)==False:
                raise RuntimeError('[!!!] Number of percentiles provided must be of the form 2^k-1')
            if args.replayQ == None:
                raise RuntimeError('[!!!] Must provide one or more percentiles for replayprob_alphabet!')
            for q in args.replayQ:
                if q < 0 or q > 100:
                    warnings.warn('Percentile q must be between 0 and 100!')
            idir_for_cover = path.join(args.dataroot, 'cover/')

            qstring=''
            for q in args.replayQ:
                qstring=qstring+str(int(q))+'-'
            odir_for_embedded = path.join(args.dataroot, 'embedded/replayprob-Q' + qstring + normStr)

            if not isdir(idir_for_cover):
                raise RuntimeError('[!!!] Input directory not found! %s' % idir_for_cover)

            if not isdir(odir_for_embedded):
                makedirs(odir_for_embedded)

#             C = args.rateK / args.rateIPD * rateF[args.rate](args.rateIPD)
            with SeedBank(args.bankPath, runID=args.runID) as bank:
                for cl in args.covlens:
                    if args.replayP == None:
#                         prob = C / rateF[args.rate](cl)
                        prob = calc_mssage_prob(cl, args.rate, rateC[args.rate])
                    else:
                        prob = args.replayP
                    chunk_size=int(np.log2(len(args.replayQ)+1))
                    msgLength = cl*chunk_size
                    if msgLength <= args.nbits:
                        msg = full_msg[:msgLength]
                    else:
                        # Now raise error if not enough message bits
                        raise RuntimeError('[!!!] Too few message bits!')

                    ifile_for_cover_prefix = 'cover-for-embed-' + ifileDistParamStr + '-n' + str(cl) + '-'
                    ofile_for_embedded_prefix = 'embedded-replayprob-' + 'Q' + qstring + ofileDistParamStr + '-n' + str(cl) + '-'


                    inputFiles = DataFile.getValidDataFiles(idir_for_cover, ifile_for_cover_prefix)
                    if not args.reuseSeeds:
                        ifilesWithSeeds = bank.getNewSeeds(inputFiles, subRunID=cl, seedBits=PRNGClass.seedSize)
                    else:
                        ifilesWithSeeds = bank.getSeedsByID(cl)

                    prefix = 'Starting embedding for cl = %d' % cl
# # Same as following with more verbose reporting
#
#                     for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds)), prefix=prefix):
                    for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds))):
                        ifile, seed = ifileWithSeed

                        df = DataFile(ifile)

                        df.addMetaData('embedderSeed', seed)
                        df.addMetaData('PRNG', args.PRNGSource)
                        coverIPDs = df.getData()
                        embedderObj=SCTimeReplayProbFromIPDsUEHQ_alphabet(coverIPDs, PRNGClass, args.replayQ, seed)
                        stegoIPDs, msgbits = embedderObj.generateIPDseq(msg, prob,cl)
                        ofile = '%s%sp%d-%03d.json' % (odir_for_embedded, ofile_for_embedded_prefix, int(10**5*prob), ifile_idx)
                        if isfile(ofile):
                            logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                            continue
                        df.addMetaData('embeddingMethod', 'replayprob_alphabet')
                        df.addMetaData('msglength', msgLength)
                        df.addMetaData('msglengths_actually_embedded', msgbits)
                        df.setData(stegoIPDs)
                        df.write(ofile, fastOutput)

        if args.replayrepeat:
            if args.numrepeat < 0:
                warnings.warn('numrepeat must be non-negative')
            idir_for_cover = path.join(args.dataroot,'embed-hq/')
            odir_for_embedded = path.join(args.dataroot, 'embedded/replayrepeat-' + str(int(args.numrepeat)) + normStr + '/')
            if not isdir(idir_for_cover):
                raise RuntimeError('[!!!] Input directory not found! %s' % idir_for_cover)                

            if not isdir(odir_for_embedded):
                makedirs(odir_for_embedded)

            with SeedBank(args.bankPath, args.runID) as bank:
                for cl in args.covlens:
                    print('[###] Starting embedding for cl=', cl)
        
                    ifile_for_cover_prefix = 'cover-for-embed-' + ifileDistParamStr + '-n' + str(cl) + '-'
                    ofile_for_embedded_prefix = 'embedded-replayrepeat-'+ str(int(args.numrepeat)) + '-' + ofileDistParamStr + '-n' + str(cl) + '-'

                    inputFiles = DataFile.getValidDataFiles(idir_for_cover, ifile_for_cover_prefix)
                    if not args.reuseSeeds:
                        ifilesWithSeeds = bank.getNewSeeds(inputFiles, subRunID=cl, seedBits=PRNGClass.seedSize)
                    else:
                        ifilesWithSeeds = bank.getSeedsByID(cl)
                    
                    prefix = 'Starting embedding for cl = %d' % cl
# # Same as following with more verbose reporting
#
#                     for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds)), prefix=prefix):
                    for ifile_idx, ifileWithSeed in progressBarIterable(tuple(enumerate(ifilesWithSeeds))):
                        ifile, seed = ifileWithSeed

                        df = DataFile(ifile)

                        df.addMetaData('embedderSeed', seed)
                        df.addMetaData('PRNG', args.PRNGSource)
                        coverIPDs = df.getData()

                        for bitnum in range(len(coverIPDs)):
                            coverIPDs[bitnum]=coverIPDs[bitnum%args.numrepeat] 
                            
                        ofile = '%s%s%03d.json' % (odir_for_embedded, ofile_for_embedded_prefix, ifile_idx)
                        if isfile(ofile):
                            logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                            continue
                        df.setData(coverIPDs)
                        df.write(ofile, fastOutput)


if __name__ == "__main__":
    import argparse
    import sys
    RNGSources = PRNG.getListOfRNGSources()

    parser = argparse.ArgumentParser(description='Embed (simulated) encrypted messages into cover traffic to produce input data for testing detection methods.')
    parser.add_argument('--PRNGSource', choices=RNGSources, default=None, help='Select which PRNG you want to use.  Choices are: %s' % str(RNGSources))
    parser.add_argument('--dataroot', metavar='DIR', type=str, default=None, help='Root directory for data hierarchy')
    parser.add_argument('--logdir', metavar='DIR', type=str, default=None, help='Directory for log files')
    parser.add_argument('--logfile', type=str, default=None, help='Name of log file for running embeddings')
    parser.add_argument('--testbedcfg', metavar='CFGFILE', type=str, default='../nextstep.cfg', help='Name of general testbed configuration file (default is ../nextstep.cfg)')
    parser.add_argument('--fileEndingRegex', type=str, default=DataFile.fileNameEndngRegex, help='A regular expression which defines what the end of the filename should look like.  Defaults to %s' % DataFile.fileNameEndngRegex)
    parser.add_argument('--humanReadable', action='store_true', default=None, help='Use fast output mode for datafile.  This base64 encodes an array which means it isn\'t really human readable.')
    # Information about message data

    parser.add_argument('--nbits', metavar='N', type=int, default=None, help='Number of bits in message file; ignored if --msgfile is specified')
    parser.add_argument('--msgfile', type=str, default=None, help='File containing message to embed; if absent, use ./urand-b<N>.txt, where N is specified by --nbits option')

    # Indication of which distribution to use
    syntheticArgs = parser.add_argument_group('Synthetic Data Ingest', 'Set options related to the use of synthetic data.')
    syntheticArgs.add_argument('--dist', type=str, default=None, help='Cover distribution to use (sellkepareto and narrow are aliases for paretoa100b095 and paretoa100b10, respectively)', choices=['pareto', 'sellkepareto', 'narrow', 'paretoa100b10', 'paretoa100b095', 'normal', 'ingest', 'markov-ingest'])

    # More about the distribution
    syntheticArgs.add_argument('--paretoa', metavar='ALPHA', type=int, default=None, help='Alpha value for Pareto distribution (int)')
    # Process beta as a string to facilitate standardization of floating-point argument
    syntheticArgs.add_argument('--paretob', metavar='BETA', type=float, default=None, help='Beta value for Pareto distribution (float)')
    syntheticArgs.add_argument('--mean', metavar='MU', type=int, help='Mean value for normal distribution')
    syntheticArgs.add_argument('--sd', metavar='SIGMA', type=float, help='Standard deviation for normal distribution')

    # Indication of which ratio to use

    parser.add_argument('rate', type=str, help='Embedding density growth rate to use', choices=embeddedRates)
    parser.add_argument('--rateIPD', type=float, metavar='N', default=None, help='Num. cover IPDs at which embedding rates are equalized')
    parser.add_argument('--rateK', type=float, metavar='k', default=None, help='Num. bits that should be embedded into <rateIPD> IPDs')

    # Indication of which embedding method to consider
    embedding_group = parser.add_argument_group(title='Embedders')
    embedding_group = embedding_group.add_mutually_exclusive_group(required=True)
    # replay and sellke embeddings not yet support for HQ randomness
    embedding_group.add_argument('--replayperc', action='store_true', help='Use replay embedding skewed to use a specified set of percentiles')
    embedding_group.add_argument('--donothing', action='store_true', help='Embedding operation is simply the identity function (does not change IPD sequence); requires sqrt embedding rate')
    embedding_group.add_argument('--replayprob', action='store_true', help='Use probabalistic replay embedding skewed to use a specified set of percentiles')
    embedding_group.add_argument('--replayprob_alphabet', action='store_true', help='Use probabalistic replay embedding with any size of message alphabet, skewed to use a specified set of percentiles')
    embedding_group.add_argument('--replayrepeat', action='store_true', help='Repeat initial segment (of length --numrepeat) of cover string')
    # embedding_group.add_argument('--sellke8to3', action='store_true', help='Use Sellke 8-to-3 embedding')

    # Details of embeddings
    parser.add_argument('--res', type=int, default=None, help='Resolution for Sellke embedding (ignored for other embeddings; default is 10)')
    parser.add_argument('--replayP', metavar='P', type=float, default=None, help='Probabilty of changing each IPD (--replayprob only; ignored if a different embedding is specified). If provided, will overide a specified embeding rate')
    parser.add_argument('--replayQ', metavar='Q', type=int, nargs='*', default=None, help='Int(s) that indicate the percentile(s) to use for skewed replay embeddings (--replayperc and --replayprob only; ignored if a different embedding is specified)')

    parser.add_argument('--numrepeat', metavar='rep', type=int, default=None, help='Int telling how long initial segment of cover to use for replayrepeat (--replayrepeat only; ignored if a different embedding is specified)')
    parser.add_argument('--covlens', metavar='L', type=int, nargs='*', default=None, help='Integer(s) giving the cover length(s) to use (read from config file unless specified here)')

    parser.add_argument('--numcovs', metavar='N', type=int, default=100, help='Number of cover sequences to use (default is 100)')
    parser.add_argument('--minidx', metavar='IDX', type=int, default=0, help='Smallest index to use (use if some files previously used for embedding)')

    seedManagement = parser.add_argument_group(title='Seed Management')
    newSeedManagement = seedManagement.add_mutually_exclusive_group()
    newSeedManagement.add_argument('--newSeeds', default=None, action='store_true', help='Generate new seeds, --runID can be used to specify a custom name.')
    newSeedManagement.add_argument('--reuseSeeds', default=None, action='store_true', help='Use existing seeds, if set specify which set of seeds via --runID')
    seedManagement.add_argument('--runID', required='--reuseSeed' in sys.argv, type=str, default=None, help='The runID used to reffer to these seeds.  When generating new seeds if no runID is provided the current date and time will be used.  This argument is required when --reuseSeeds is set')
    seedManagement.add_argument('--bankPath', default=None, type=str, help='The path to the seed bank for the embedders.')

    args = parser.parse_args()

    if args.dist is None and args.ingestFiles is None:
        parser.print_usage()
        sys.exit()

    ########
    main(args)
    ########

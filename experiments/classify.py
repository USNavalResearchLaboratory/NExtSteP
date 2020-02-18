#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
from scipy.stats import ttest_ind, mannwhitneyu, chisquare
from multiprocessing import Pool, cpu_count
from os import makedirs, listdir, path
from os.path import isfile, isdir
from math import modf

import sys
sys.path.append('../') # this is to find libStealthy!
from libStealthy.util import rateF, calc_mssage_length, calc_mssage_prob, calc_rate_constants, window, DataFile, progressBarIterable, setArgsFromConfig
from libStealthy.classifiers import welcht, mwu, filler, sgfiller, cce, localMinCCE, Entropy, chisq
from libStealthy.srl import statisticNu0IPD, statisticSignedNu0IPD, calcProbLow
from libStealthy.stego.fromIPDs import SCTimeReplayFromIPDsUE
from libStealthy.exceptions import EncodingError
from libStealthy.pRandom import PRNG

import configparser
import numpy as np
import warnings
import random
import json

embeddedRates = ['const', 'log', 'cubrt', 'sqrt', '2rt3', 'linear']

def processDetection(ifileDetector):
    ifile, detector = ifileDetector
    try:
        df = DataFile(ifile)
        IPDs = df.getData()
        stat = detector.classify(IPDs)
        return stat
    except FileNotFoundError:
        print('[~~~] Unable to find traffic file:', ifile)


##################

def main(args):
    try:
        args = setArgsFromConfig(args.testbedcfg, args, 'classify')
        try:
            args = setArgsFromConfig(args.testbedcfg, args)
        except:
            raise
    except FileNotFoundError:
        print('[!!!] Unable to find configuration file %s, aborting.' % args.testbedcfg)
        sys.exit(1)
    except KeyError:
        print('[~~~] No file-specific configs were found.')

#     print('[###] Log file will be', args.logfile)
#     print('[###] Data hierarchy rooted at', args.dataroot)
#     print('[###] Logs will be in', args.logdir)

    args.covlens = sorted(args.covlens)

    if not isdir(args.dataroot):
        raise RuntimeError('[!!!] Dataroot %s is not a directory!' % args.dataroot)
    if args.logdir is None:
        args.logdir = path.join(args.dataroot, 'logs/')
    if not isdir(args.logdir):
         makedirs(args.logdir) 

    if args.fillerP is not None and not (0.0 < args.fillerP < 1.0):
        raise RuntimeError('[!!!] fillerP is not between 0.0 and 1.0!')
    
    if 'ingest' in args.dist: # No third way
        args.dataroot = path.join(args.dataroot, 'ingest/')
    else:
        args.dataroot = path.join(args.dataroot, 'hq/')
    
    print('[###] Using', len(args.covlens), 'cover lengths:', args.covlens)

    numRefIPDs = args.reflen
    rate = args.rate

    # This is the number of bits to embed for most methods
    # For time-replay embedding, this will give the number of IPDs to change

    k0 = args.rateK 
    N0 = args.rateIPD 
    rateC = calc_rate_constants(N0, k0)

    normStr = '-k' + str(k0) + '-n' + str(N0)


    # These should give the same results, but using the second format
    # to parallel other scripts
    #
#     embedding_prob=lambda N: (float(k0)/rateF[args.rate](N0)) * (rateF[args.rate](N) / float(N))
    embedding_prob=lambda N: calc_mssage_prob(N, rate, rateC[rate])
        
    if args.embedding == 'replayperc':
        embeddingTypeStr = 'replayperc-Q' + str(int(args.replayQ[0]))
        MLs = [calc_mssage_length(cl, args.rate, rateC[args.rate]) for cl in args.covlens]
    elif args.embedding == 'donothing':
        embeddingTypeStr = 'donothing'
        if args.rate != 'sqrt':
            raise RuntimeError('[!!!] Do-nothing embedding requires sqrt rate!')
        MLs = [calc_mssage_length(cl, args.rate, rateC[args.rate]) for cl in args.covlens]
    elif args.embedding == 'sellke8to3':
        embeddingTypeStr = 'sellke-L8-n3'
        MLs = [calc_mssage_length(cl, args.rate, rateC[args.rate]) for cl in args.covlens]
    # Embeddings with non-deterministic message length
    elif args.embedding =='replayprob':
        embeddingTypeStr='replayprob-Q'+str(int(args.replayQ[0]))
        Probs=[embedding_prob(cl) for cl in args.covlens]
    elif args.embedding == 'replayrepeat':
        embeddingTypeStr='replayrepeat-'+str(int(args.numrepeat))

    elif args.embedding =='replayprob_alphabet':
        qstring=''
        for q in args.replayQ:
            qstring=qstring+str(int(q))+'-'
        embeddingTypeStr='replayprob-Q'+qstring
        Probs=[embedding_prob(cl) for cl in args.covlens]

    else:
        raise RuntimeError('Unsupported embedding type!')


    if args.dist == 'sellkepareto' or args.dist == 'paretoa100b095':
        distParamStr = 'sellkepareto'
        ifileDistParamStr = 'sellke'
    elif args.dist == 'narrow' or args.dist == 'paretoa100b10':
        distParamStr = 'pareto-a100-b10'
        ifileDistParamStr = 'pareto-a100-b10'
    elif args.dist == 'pareto':
        alpha = args.paretoa
        beta = args.paretob
        if beta <= 0.0:
            raise ValueError('[!!!] Beta must be positive!')
        ifileDistParamStr = 'pareto-a%d-b%s' % (alpha, str(beta).replace('.', 'p'))
        distParamStr = ifileDistParamStr
    elif args.dist == 'normal':
        distmean = args.mean
        sd = args.sd
        if sd <= 0.0:
            raise ValueError('[!!!] Std. dev. must be positive!')
        ifileDistParamStr = 'normal-m%d-s%s' % (distmean, str(sd).replace('.', 'p'))
        distParamStr = ifileDistParamStr
    elif args.dist == 'ingest':
        distParamStr = 'ingest'
        ifileDistParamStr = 'ingest'
    elif args.dist == 'markov-ingest':
        distParamStr = 'markov-ingest'
        ifileDistParamStr = 'markov-ingest'


    if not isfile(path.join(args.logdir, args.logfile)):
        with open(path.join(args.logdir, args.logfile), 'w') as logstream:
            logstream.write('\n')

    with open(path.join(args.logdir, args.logfile), 'a') as logstream:
        # set up classifier to be used and provide classifier secfic
        # arguments.
        if args.welcht:
            detector = welcht()

        elif args.mwu:
            detector = mwu()

        elif args.chisq:
            detector = chisq(args.chisqBins)

        elif args.cce:
            if not args.cceLocalMin:
                detector = cce(args.cceSubSeqLen, args.cceBins)
            else:
                detector = localMinCCE(args.cceSubSeqLen,args.cceBins)

        elif args.filler:
            if args.fillerP is not None:
                inputRefFiles = tuple(True)
            detector = filler(args.fillerP)

        elif args.sgfiller:
            if args.fillerP is not None:
                inputRefFiles = tuple(True)
            detector = sgfiller(args.fillerP)


        # create our process pool for parallelization.
        with Pool(args.numProcs) as p:
            for i, cl in enumerate(args.covlens):
                # print("[###] Starting classify cl = %d\r" % cl, end='')
                idir_for_ref = args.dataroot + 'reference/'
                ifile_for_ref_prefix = 'cover-for-ref-' + ifileDistParamStr + '-n' + str(numRefIPDs) + '-'
                idir_for_plain = args.dataroot + 'plain/'
                ifile_plain_prefix = 'cover-plain-' + ifileDistParamStr + '-n' + str(cl) + '-'

                if args.embedding == 'replayperc':
                    ml=MLs[i]
                    idir_for_embedded = args.dataroot + 'embedded/' + embeddingTypeStr + normStr + '/'
                    ifile_for_embedded_prefix = 'embedded-' + embeddingTypeStr + '-m' + str(ml) + '-' + distParamStr + '-n' + str(cl) + '-'
                    ofile_middle = '-' + distParamStr + '-r' + str(numRefIPDs) + '-c' + str(cl) + '-m' + str(ml)
                    odir_for_classout = args.dataroot + 'classout/' + embeddingTypeStr + normStr + '/'
                    error_string= 'ML = '+str(ml)

                elif args.embedding == 'donothing':
                    ml=MLs[i]
                    idir_for_embedded = args.dataroot + 'embedded/' + embeddingTypeStr + '/'
                    ifile_for_embedded_prefix = 'embedded-' + embeddingTypeStr + '-m' + str(ml) + '-' + distParamStr + '-n' + str(cl) + '-'
                    ofile_middle = '-' + distParamStr + '-r' + str(numRefIPDs) + '-c' + str(cl) + '-m' + str(ml)
                    odir_for_classout = args.dataroot + 'classout/' + embeddingTypeStr + '/'
                    error_string= 'ML = '+str(ml)

                elif args.embedding == 'replayprob':
                    prob=Probs[i]
                    idir_for_embedded = args.dataroot + 'embedded/' + embeddingTypeStr + normStr + '/'
                    ifile_for_embedded_prefix = 'embedded-' + embeddingTypeStr + '-' + distParamStr +'-n' + str(cl)+'-p'+ str(int(10**5*prob)) +'-'
                    ofile_middle = '-' + distParamStr + '-r' + str(numRefIPDs) + '-c' + str(cl) + '-p' + str(int(10**5*prob))
                    odir_for_classout = args.dataroot + 'classout/' + embeddingTypeStr + normStr + '/'
                    error_string='p='+str(int(10**5*prob)/10**5)

                elif args.embedding =='replayprob_alphabet':
                    prob=Probs[i]
                    idir_for_embedded = args.dataroot + 'embedded/' + embeddingTypeStr + normStr + '/'
                    ifile_for_embedded_prefix = 'embedded-' + embeddingTypeStr + distParamStr +'-n' + str(cl)+'-p'+ str(int(10**5*prob)) +'-'
                    ofile_middle = '-' + distParamStr + '-r' + str(numRefIPDs) + '-c' + str(cl) + '-p' + str(int(10**5*prob))
                    odir_for_classout = args.dataroot + 'classout/' + embeddingTypeStr + normStr + '/'
                    error_string='p='+str(int(10**5*prob)/10**5)

                elif args.embedding =='replayrepeat':
                    idir_for_embedded = args.dataroot + 'embedded-hq/' + embeddingTypeStr + normStr + '/'
                    ifile_for_embedded_prefix = 'embedded-' + embeddingTypeStr + '-'+distParamStr +'-n' + str(cl)+'-'
                    ofile_middle = '-' + distParamStr + '-r' + str(numRefIPDs) + '-c' + str(cl) 
                    odir_for_classout = args.dataroot + 'classout-hq/' + embeddingTypeStr + normStr + '/'
                    error_string='numrepeats ='+str(int(args.numrepeat))

                else:
                    ml=MLs[i]
                    idir_for_embedded = args.dataroot + 'embedded/' + embeddingTypeStr + normStr + '/'
                    ifile_for_embedded_prefix = 'embedded-' + embeddingTypeStr + '-m' + str(ml) + '-' + distParamStr + '-n' + str(cl) + '-'
                    ofile_middle = '-' + distParamStr + '-r' + str(numRefIPDs) + '-c' + str(cl) + '-m' + str(ml)
                    odir_for_classout = args.dataroot + 'classout/' + embeddingTypeStr + normStr + '/'
                    error_string= 'ML = '+str(ml)

                if not isdir(odir_for_classout):
                    makedirs(odir_for_classout)

                classifier_outs = [[] for j in range(args.numrefs)]

                inputRefFiles = list(DataFile.getValidDataFiles(idir_for_ref, ifile_for_ref_prefix))[:args.numrefs]
                inputPlainFiles = list(DataFile.getValidDataFiles(idir_for_plain, ifile_plain_prefix))
                inputEmbeddedFiles = list(DataFile.getValidDataFiles(idir_for_embedded, ifile_for_embedded_prefix))

                ### Check whether we have some of each kind of required data files
                validRef = any(inputRefFiles)
                validPlain = any(inputPlainFiles)
                validEmbedded = any(inputEmbeddedFiles)                
                # We'll only continue if we have valid input files of each required types
                if validRef and validPlain and validEmbedded:
                    for j, ref_ifile in progressBarIterable(list(enumerate(inputRefFiles))):
# # Same as following with more verbose reporting
#
#                     for j, ref_ifile in progressBarIterable(list(enumerate(inputRefFiles)),prefix="[###] Starting classify cl = %d " % cl):
                        try:
                            refDF = DataFile(ref_ifile)
                            refIPDs = refDF.getData()
                            detector.setup(refIPDs)

                            # use a generator expression to create ifile/detector pairings
                            plainIfileGenerator = ((ifile, detector) for ifile in inputPlainFiles)
                            # do plain detection in parallel.
                            outs = p.map(processDetection, plainIfileGenerator)
                            # Indicate that this was just plain cover traffic
                            # Might use star notation to unpack the tuple
                            # classifier_outs[j].extend([(*stat, 0) for stat in outs])                            
                            
                            
                            classifier_outs[j].extend([list(stat) + [0] for stat in outs])

                            # use a generator expression to create ifile/detector pairings                            
                            embeddedIfileGenerator = ((ifile, detector) for ifile in inputEmbeddedFiles)
                            # do embedded detection in parallel.
                            outs = p.map(processDetection, embeddedIfileGenerator)
                            # Indicate that this traffic had embedded message
                            # Might use star notation to unpack the tuple                            
                            # classifier_outs[j].extend([(*stat, 1) for stat in outs])
                            classifier_outs[j].extend([list(stat) + [1] for stat in outs])
                        except FileNotFoundError:
                            logstream.write('FILE: Unable to find ifile: ' + ref_ifile + '\n')
                            print('[~~~] Unable to find reference file:', ref_ifile)

                    ofile = odir_for_classout + detector.detectorStr + ofile_middle + '.json'
                    classOutDF = DataFile(fast=False)
                    # classOutDF.addMetaData('refSeed', refDF.getMetaData('seed'))
                    # classOutDF.addMetaData('')
                    classOutDF.addMetaData('detectorStr', detector.detectorStr)
                    if detector.detectorStr in ['cce','localMinCCE'] and len(args.cceBins) == 1:
                        classOutDF.addMetaData('cceDepth',args.cceSubSeqLen)
                        classOutDF.addMetaData('cceBins',args.cceBins[0])
                        ofile = odir_for_classout + detector.detectorStr + '-d' + str(args.cceSubSeqLen) + '-n' + str(args.cceBins[0]) + ofile_middle + '.json'
                    elif detector.detectorStr in ['cce','localMinCCE'] and type(args.cceBins) is int:
                        classOutDF.addMetaData('cceDepth',args.cceSubSeqLen)
                        classOutDF.addMetaData('cceBins',args.cceBins)
                        ofile = odir_for_classout + detector.detectorStr + '-d' + str(args.cceSubSeqLen) + '-n' + str(args.cceBins) + ofile_middle + '.json'
                    elif detector.detectorStr in ['cce','localMinCCE']:
                        classOutDF.addMetaData('cceDepth',args.cceSubSeqLen)
                        classOutDF.addMetaData('cceBins','custom')
                        ofile = odir_for_classout + detector.detectorStr + '-d' + str(args.cceSubSeqLen) + '-cust' + ofile_middle + '.json'
                    
                    if detector.detectorStr in ['chisq']:
                        classOutDF.addMetaData('chisqBins',args.chisqBins)
                        ofile = odir_for_classout + detector.detectorStr + '-numBins' + str(args.chisqBins) + ofile_middle + '.json'

                    classOutDF.setData(classifier_outs)
                    classOutDF.write(ofile)
                else:
                    print('[~~~] Skipping', detector.detectorStr, 'for CL =', cl, 'and', error_string)
                    if not validRef:
                        logstream.write('FILE: Unable to find any reference file, starting with ' + ifile_for_ref_prefix + '\n')
                        print('[!!!] Unable to find any reference file, starting with:', ifile_for_ref_prefix)
                    if not validPlain:
                        logstream.write('FILE: Unable to find any plain-traffic file, starting with ' + ifile_plain_prefix + '\n')
                        print('[!!!] Unable to find any plain-traffic file, starting with:', ifile_plain_prefix)
                    if not validEmbedded:
                        logstream.write('FILE: Unable to find any embedded-traffic file, starting with ' + ifile_for_embedded_prefix + '\n')
                        print('[!!!] Unable to find any embedded-traffic file, starting with:', ifile_for_embedded_prefix)
            p.close()
            p.join()

if __name__ == "__main__":
    import argparse
    from pprint import pprint
    parser = argparse.ArgumentParser(description='Run classifiers for studies of possible square-root laws')

    parser.add_argument('--dataroot', metavar='DIR', type=str, default=None, help='Root directory for data hierarchy')
    parser.add_argument('--logdir', metavar='DIR', type=str, default=None, help='Directory for log files')
    parser.add_argument('--logfile', type=str, default=None, help='Name of log file for running classifiers')
    parser.add_argument('--testbedcfg', metavar='CFGFILE', type=str, default='../nextstep.cfg', help='Name of general testbed configuration file (default is ../nextstep.cfg)')
    parser.add_argument('--numProcs', type=int, default=None, help='How many processes should python use to run your classifiers.  Note anything over %d may decrease performance on this machine. Addtionally, increasing the process count will increase memory usage (because you are now running n pocesses in parallel all alocating their own memory).' % cpu_count())

    # Flags for different possible classifiers
    # Must pick exactly one of these
    classifier_group = parser.add_mutually_exclusive_group(required=True)
    classifier_group.add_argument('--mwu', action='store_true', help='Use Mann--Whitney U classifier')
    classifier_group.add_argument('--welcht', action='store_true', help='Use Welch\'s t test classifier') 
    classifier_group.add_argument('--chisq', action='store_true', help='Use chi-squared-test classifier')
    classifier_group.add_argument('--filler', action='store_true', help='Use Filler\'s frequency-based classifier.  Unless the --fillerP argument is given, the reference string is used to compute the probability of a 0 in the list.')
    classifier_group.add_argument('--sgfiller', action='store_true', help='Use signed version of Filler\'s frequency-based classifier.  Unless the --sgfillerP argument is given, the reference string is used to compute the probability of a 0 in the list.')
    classifier_group.add_argument('--cce', action='store_true', help='Use CCE based classifier')

    # Indication of which distribution to use
    parser.add_argument('dist', type=str, help='Cover distribution to use (sellkepareto and narrow are aliases for paretoa100b095 and paretoa100b10, respectively)', choices=['pareto', 'sellkepareto', 'narrow', 'paretoa100b10', 'paretoa100b095', 'normal', 'ingest', 'markov-ingest'])

    # More about the distribution
    parser.add_argument('--paretoa', metavar='ALPHA', type=int, default=None, help='Alpha value for Pareto distribution (int;)')
    parser.add_argument('--paretob', metavar='BETA', type=float, default=None, help='Beta value for Pareto distribution (float)')
    parser.add_argument('--mean', metavar='MU', type=int, help='Mean value for normal distribution')
    parser.add_argument('--sd', metavar='SIGMA', type=float, help='Standard deviation for normal distribution')

    # Specification of probLow to use for Filler's detector (if not using reference string)
    parser.add_argument('--fillerP', metavar='ProbLow', type=float, help='Underlying probability of a below-median value in the IPD sequence for use with Filler\'s detector')
    parser.add_argument('--sgfillerP', metavar='ProbLow', type=float, help='Underlying probability of a below-median value in the IPD sequence for use with signed version of Filler\'s detector')

    # CCE specific options
    CCEOptions = parser.add_argument_group(title='CCE Options')
    CCEOptions.add_argument('--cceBins', nargs='*', default = [5], help = "If a single int is supplied then it will be used as the number of bins.  If a list of numbers are supplied they will be used as the edges of the bins.  The bin count is also the branching factor of the tree.")
    CCEOptions.add_argument('--cceSubSeqLen', type=int, default = 50, help = 'The length of the sub-sequences made from the original sequence of IPDs')
    CCEOptions.add_argument('--cceLocalMin', action='store_true', help = 'Use the local minimum as an early stoping factor, this uses the simple CCE algorthim.')

    # Chi Square specific options
    ChiSqOptions = parser.add_argument_group(title='Chisq options')
    ChiSqOptions.add_argument('--chisqBins', type=int, default = 20, help = 'The number of bins to sort refIPDs into using Chisq analysis')

    # Rate function and constants
    parser.add_argument('rate', type=str, help='Embedding density growth rate to use', choices=embeddedRates)
    parser.add_argument('--rateIPD', type=float, metavar='N', default=None, help='Num. cover IPDs at which embedding rates are equalized')
    parser.add_argument('--rateK', type=float, metavar='k', default=None, help='Num. bits that should be embedded into <rateIPD> IPDs')

    # Indication of which embedding method to consider
    parser.add_argument('embedding', type=str, choices=['replay', 'replayperc', 'donothing', 'replayprob', 'sellke8to3', 'replayprob_alphabet', 'replayrepeat'], help='Embedding method to consider')

    # Details of embeddings
    parser.add_argument('--replayQ', metavar='Q', type=int, nargs='*', default=[90], help='Percentile (integer) to use for skewed replay embeddings (replayperc embedding; ignored if a different embedding is specified)')
    parser.add_argument('--numrepeat', metavar='rep', type=int, default=5, help='replayrepeat only; (ignored if a different embedding is specified) default 5')

    # Cover lengths to include
    parser.add_argument('--covlens', metavar='L', type=int, nargs='*', default=None, help='Integer(s) giving the cover length(s) to use (read from config file unless specified here)')

    # Numbers of covers and reference strings
    parser.add_argument('--numcovs', metavar='N', type=int, default=100, help='Number of cover sequences to use')
    parser.add_argument('--numrefs', metavar='N', type=int, default=None, help='Number of reference sequences to use')
    parser.add_argument('--refsspec', action='store_true', help='Use "special" reference IPD sequences with all classifiers (overrides --numrefs)')

    # Number of reference IPDs

    parser.add_argument('--reflen', type=int, default=None, help='Number of reference IPDs to use')

    ########

    args = parser.parse_args()

    
    main(args)

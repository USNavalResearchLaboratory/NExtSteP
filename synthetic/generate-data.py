#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
'''
Construct and write to disk files containing IPD sequences
These are the reference sequences for use by the detectors,
the cover sequences that will not have embeddings ("plain"),
and the cover sequences that will have embeddings ("cover").

Each output file has one JSON-encoded list of IPDs.

File names end with -nN-XXX.json, where N is the number of IPDs
and XXX is an index.

num_covers and num_refs in the code below indicate how many cover
and reference files should be generated.  (The number of covers is
for each type, plain and embed.)
'''

from Crypto.Random import random as HQrandom
from math import sqrt, log, log10
from os.path import isfile, isdir
from os import makedirs, path
from random import random
	
import sys
sys.path.append('../') # this is for libStealthy
from libStealthy.util import hqrand, DataFile, printProgress, window, equalAreaBins, setArgsFromConfig, SeedBank
from libStealthy.coverDist import Pareto, ParetoNarrow, ParetoSellke, MarkovIngest, Normal
from libStealthy.pRandom import PRNG

import multiprocessing
import configparser
import numpy as np
import random
import json
import time

def genParetoIPDsSpecial(numIPDs, seed, PRNGSource, alpha=100, beta=10, maxIPD=500):
    '''
    Generate a list of IPDs matching Pareto dist. as closely as possible

    Used to generate "special" reference sequences

    numIPDs specifies the number of IPDs that should be generated;
    these are returned in a list.

    alpha and beta are the parameters of the Pareto distribution being
    sampled.
    '''
    np.random = PRNG.getRng(PRNGSource)
    np.random.seed(seed)
    resFactor = 10 # 1 gives ms, 10 gives 10^{-4} s
    maxIPD = 500 # In ms, regardless of resFactor
    numRands = numIPDs #100000

    F = lambda x : 1 - (float(alpha)/float(x))**beta
    maxRand = F(maxIPD + 1.0/resFactor)
    minRand = 0.0

    points = np.linspace(minRand, maxRand, num=numRands, endpoint=True)

    Finv = lambda x : alpha * ((1 / (1.0 - x))**(1.0/beta))

    ipds = [int(resFactor*Finv(x)) for x in points]

    if len(ipds) != numRands:
        raise RuntimeError('[!!!] Wrong number of IPDs produced!')

    return ipds


def plotIPDs():
    '''
    Plot IPDs from a JSON-ecoded file
    '''
    import matplotlib.pyplot as plt

    betastr = '10'
    ifile_ref = 'cover-for-ref-pareto-a100-b' + betastr + '-n1000000-000.json'
    idir = '../../data/sqrtdata/reference/'
    titleStr = 'Pareto-Distributed IPDs (a=100, b=' + betastr + ')'
    ifile = idir + ifile_ref
    nsamples = 1000000

    resFactor = 10 # 1 gives ms, 10 gives 10^{-4} s

    maxIPD = 500 # In ms, regardless of resFactor

    _, ax = plt.subplots()

    bins = [i for i in range(resFactor*90, resFactor*(maxIPD+10))]

    with open(ifile, 'r') as istream:
        IPDs = json.loads(istream.readline())
    ax.hist(IPDs, bins, normed=1)
    ax.set_xlabel('IPD in 1e-4 seconds')
    ax.set_ylabel('Normalized frequency')

    plt.title(titleStr)

    fig_odir = '../../data/sqrtdata/samples/plots/'

    fig_ofile = fig_odir + 'n-' + str(nsamples) + 'cover-pareto-a100-b' + betastr + '.pdf'
    plt.savefig(fig_ofile)

    plt.show()


def genFile(ofile, numIPDs, dist, seed=None, fast=True):
    '''
    This funciton is used to generate the files in parallel.  Its
    passed to multiprocessing.Pool.apply_async which calles it with
    the argument tuple supplied to it.
    
    Args:
        ofile (str):
            The path to the output file.
        numIPDs (int):
            The number of IPDs you would like generated.
        distArgs (BaseCoverDist or tuple):
            An instance of a class implmenting BaseCoverDist.  If it
            is a tuple then it must contain the constructor for a 
            distrubuiton class and a tuple containing the arguments
            to be provided to the constructor.
        seed (int or bytes):
            The seed for the PRNG.
        fast (bool):
            Sould fast encoding be used.
    '''
    # if dist is tuple containing the constructor for a distrubution 
    # class and another tuple containing the arguments we want to 
    # provide to that constructor.  
    # For instance: 
    #     (Pareto, (alpha, beta, args.PRNGSource)) 
    # where Pareto is the constructor for the Pareto distrubution and 
    # the tuple is he arguments provided. We use tuple expansion to 
    # fill in the arguments. 
    if isinstance(dist, tuple):
        dist, args = dist 
        dist = dist(*args) 
    
    df = DataFile(fast=fast)
    df.mergeMetaData(dist.getMetaData())
    if seed is None:
        seed = HQrandom.getrandbits(256)


    IPDseq = dist.generateIPDseq(numIPDs, seed)
    df.setData(IPDseq)
    df.write(ofile)
    return (numIPDs, ofile, seed)

def main(args):    
    try:
        args = setArgsFromConfig(args.testbedcfg, args, 'gencovers')
        try:
            args = setArgsFromConfig(args.testbedcfg, args)
        except:
            raise
    except FileNotFoundError:
        print('[!!!] Unable to find configuration file %s, aborting.' % args.testbedcfg)
        sys.exit(1)
    except KeyError:
        print('[~~~] No file specfic configs were found.')

# # Various possible info. messages
#
#     print('[###] Log file will be', args.logfile)
    print('[###] Cover Dist will be', args.parserID)
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
    if 'ingest' in args.parserID:
        args.dataroot = path.join(args.dataroot, 'ingest')
    else:
        args.dataroot = path.join(args.dataroot, 'hq')

    if args.parserID == 'sellke':
        ofileDistParamStr = 'sellke'
        dist = (ParetoSellke, (args.PRNGSource,))
    elif args.parserID == 'narrow':
        ofileDistParamStr = 'pareto-a100-b10'
        dist = (ParetoNarrow, (args.PRNGSource,))
    elif args.parserID == 'pareto':
        alpha = args.paretoa
        beta = args.paretob
        if beta <= 0.0:
            raise ValueError('Beta must be positive!')
        ofileDistParamStr = 'pareto-a%d-b%s' % (alpha, str(beta).replace('.', 'p'))
        dist = (Pareto, (alpha, beta, args.PRNGSource))
    elif args.parserID == 'normal':
        distmean = args.mean
        sd = args.sd
        if sd <= 0.0:
            raise ValueError('Standard deviation must be positive!')
        ofileDistParamStr = 'normal-m%d-s%s' % (distmean, str(sd).replace('.', 'p'))
        dist = (Normal, (distmean, sd, args.PRNGSource,))
    elif args.parserID == 'markov-ingest':
        # Construct our corpus for the markov chain.
        corpus = []
        for corpusFile in args.corpusPath:
            with open(corpusFile) as fp:
                corpus.extend(json.load(fp))
        with SeedBank(args.bankPath, args.runID) as sb:
            if args.reuseSeeds:
                seed = next(sb.getSeedsByID('markovSeed'))[1]
            else:
                seed = next(sb.getNewSeeds('markovSeed','markovSeed'))[1]

            dist = MarkovIngest(corpus, args.stateSize, args.lookbackSize, args.bins, args.disjointLookback, args.PRNGSource, )
        ofileDistParamStr = 'markov-ingest'

    # setup path and create them if needed.
    odir_for_plain = path.join(args.dataroot, 'plain/')
    odir_for_embed = path.join(args.dataroot, 'cover/')
    odir_for_ref = path.join(args.dataroot, 'reference/')
    odir_for_refspec = path.join(args.dataroot, 'refspecial/')
    for d in (odir_for_plain, odir_for_embed, odir_for_ref, odir_for_refspec):
        if not isdir(d):
            makedirs(d)

    numRefIPDs = args.reflens

    if args.minidx < 0:
        raise ValueError('Minimum index must be non-negative!')

    ofile_for_ref_prefix = 'cover-for-ref-' + ofileDistParamStr + '-n' + str(numRefIPDs) + '-'
    ofile_for_refspec_prefix = 'refspecial-' + ofileDistParamStr + '-n' + str(numRefIPDs) + '-'

    # This is for actually generating data
    num_covers = args.numcovs
    num_refs = args.numrefs

    with open(path.join(args.logdir, args.logfile), 'a') as logstream:
        if args.refsspec:
            # Generate special cases, not realistic traffic
            
            if alpha != 100 or beta != 10:
                raise RuntimeError('[!!!] Special reference sequences currently assume alpha of 100 and beta of 10! alpha = ' + str(alpha) + ' and beta = ' + str(beta))

            paretoMedian = 100.0 / ((0.5 + 0.5 * (100.0 / 500.0)**(beta))**(1.0/beta))
            paretoMean = 100.0 * (5**beta - 5) * beta / float((5**beta - 1) * (beta - 1))

            ofile_idx = '000'
            ofile = odir_for_refspec + ofile_for_refspec_prefix + ofile_idx + '.json'
            if isfile(ofile):
                logstream.write('FILE: Skipping existing file ' + ofile + '\n')
            else:
                df = DataFile()
                df.addMetaData('refLen', numRefIPDs)
                df.addMetaData('PRNG', args.PRNGSource)
                df.addMetaData('paretoMedian', paretoMedian)
                medIPDseq = numRefIPDs * [int(round(10 * paretoMedian))]
                df.setData(medIPDseq)
                df.write(ofile)
                print('[###] For median reference string, (min, max) is', min(medIPDseq), max(medIPDseq))

            ofile_idx = '001'
            ofile = odir_for_refspec + ofile_for_refspec_prefix + ofile_idx + '.json'
            if isfile(ofile):
                logstream.write('FILE: Skipping existing file ' + ofile + '\n')
            else:
                df = DataFile()
                df.addMetaData('refLen', numRefIPDs)
                df.addMetaData('PRNG', args.PRNGSource)
                df.addMetaData('paretoMean', paretoMean)
                meanIPDseq = numRefIPDs * [int(round(10 * paretoMean))]
                df.setData(medIPDseq)
                df.write(ofile)
                print('For mean reference string, (min, max) is', min(meanIPDseq), max(meanIPDseq))

            ofile_idx = '002'
            ofile = odir_for_refspec + ofile_for_refspec_prefix + ofile_idx + '.json'
            if isfile(ofile):
                logstream.write('FILE: Skipping existing file ' + ofile + '\n')
            else:
                df = DataFile()
                seed = HQrandom.getrandbits(256)
                df.addMetaData('seed', seed)
                df.addMetaData('PRNG', args.PRNGSource)
                df.addMetaData('refLen', numRefIPDs)
                incIPDseq = genParetoIPDsSpecial(numRefIPDs, args.PRNGSource, seed)
                df.setData(medIPDseq)
                df.write(ofile)
                print('[###] For increasing reference string, (min, max) is', min(incIPDseq), max(incIPDseq))

            ofile_idx = '003'
            ofile = odir_for_refspec + ofile_for_refspec_prefix + ofile_idx + '.json'
            if isfile(ofile):
                logstream.write('FILE: Skipping existing file ' + ofile + '\n')
            else:
                df = DataFile()
                seed = HQrandom.getrandbits(256)
                df.addMetaData('seed', seed)
                df.addMetaData('PRNG', args.PRNGSource)
                df.addMetaData('refLen', numRefIPDs)
                decIPDseq = genParetoIPDsSpecial(numRefIPDs, args.PRNGSource, seed)
                decIPDseq.reverse()
                df.setData(decIPDseq)
                df.write(ofile)
                print('[###] For decreasing reference string, (min, max) is', min(decIPDseq), max(decIPDseq))

        else:
            # Not doing special cases
            
            # playing with where we create the SeedBank could yeild
            # perf improvements
            with SeedBank(args.bankPath, args.runID) as sb:
                # Here we setup a process pool to parallelize cover
                # creation.
                with multiprocessing.Pool(args.numProcs) as p:
                    jobs = []
                    if args.refs:
                        # We're creating reference sequences
                        print('[###] Doing', numRefIPDs, 'reference IPDs')
                        ofiles = []
                        for ofile_idx in range(num_refs):
                            ofile = '%s%s%03d.json' % (odir_for_ref, ofile_for_ref_prefix, ofile_idx)
                            if isfile(ofile):
                                logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                                continue
                            ofiles.append(ofile)
                        # Add job to the pool's queue to do and move on
                        # to the next job.
                        if args.reuseSeeds:
                            filesAndSeeds = sb.getSeedsByID('ref')
                        else:
                            filesAndSeeds = ((ofile, None) for ofile in ofiles)
                        for ofile, seed in filesAndSeeds:
                            job = p.apply_async(genFile, (ofile, numRefIPDs, dist, seed, fastOutput))
                            jobs.append(job)

                    else:  # Not doing reference strings
                        
                        print('[###] Using', len(args.covlens), 'cover lengths:', args.covlens)
                        for numCoverIPDs in args.covlens:
                            ofile_plain_prefix = 'cover-plain-' + ofileDistParamStr + '-n' + str(numCoverIPDs) + '-'
                            ofile_for_embed_prefix = 'cover-for-embed-' + ofileDistParamStr + '-n' + str(numCoverIPDs) + '-'

                            ofiles = []
                            for ofile_idx in range(num_covers):
                                ofile = '%s%s%03d.json' % (odir_for_plain, ofile_plain_prefix, ofile_idx)
                                if isfile(ofile):
                                    logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                                else:
                                    ofiles.append(ofile)
                            if args.reuseSeeds:
                                filesAndSeeds = sb.getSeedsByID(str(numCoverIPDs))
                            else:
                                filesAndSeeds = ((ofile, None) for ofile in ofiles)
                            for ofile, seed in filesAndSeeds:
                                # Add job to the pool's queue to do and
                                # move on to the next job.
                                job = p.apply_async(genFile, (ofile, numCoverIPDs, dist, seed, fastOutput))
                                jobs.append(job)
                                ofile = '%s%s%03d.json' % (odir_for_embed, ofile_for_embed_prefix, ofile_idx)

                            ofiles = []
                            for ofile_idx in range(num_covers):
                                ofile = '%s%s%03d.json' % (odir_for_embed, ofile_for_embed_prefix, ofile_idx)
                                if isfile(ofile):
                                    logstream.write('FILE: Skipping existing file ' + ofile + '\n')
                                else:
                                    ofiles.append(ofile)
                            if args.reuseSeeds:
                                filesAndSeeds = sb.getSeedsByID(str(numCoverIPDs))
                            else:
                                filesAndSeeds = ((ofile, None) for ofile in ofiles)
                            for ofile, seed in filesAndSeeds:
                                # Add job to the pool's queue to do and
                                # move on to the next job.
                                job = p.apply_async(genFile, (ofile, numCoverIPDs, dist, seed, fastOutput))
                                jobs.append(job)
                    p.close()            
                    # create progress bar
                    jobsReady = len([job for job in jobs if job.ready()])
                    totalJobs = len(jobs)
                    lenWritten = 0
                    while not jobsReady == totalJobs:
                        lenWritten = printProgress(jobsReady, totalJobs)
                        time.sleep(.2)
                        jobsReady = len([job for job in jobs if job.ready()])
                    # write whitespace to cover up the progress bar
                    # and reset the cursor
                    sys.stdout.write((' ' * lenWritten) + '\r')
                    p.join()
                    # Colate all seeds used and store them in theseedbank.
                    seedsUsed = {cl:{} for cl in args.covlens}
                    seedsUsed[numRefIPDs] = {}
                    for job in jobs:
                        cl, ofile, seed = job.get()
                        seedsUsed[cl][ofile] = seed
                    sb.bank['data'][sb.runID].update(seedsUsed)


if __name__ == "__main__":

    import argparse

    RNGSources = PRNG.getListOfRNGSources()

    parser = argparse.ArgumentParser(description='Generate cover traffic for reference, plain traffic, and embedding.')
    parser.add_argument('--PRNGSource', choices=RNGSources, default=None, help='Select which PRNG you want to use default is ChaCha20.  Choices are: %s' % str(RNGSources))
    parser.add_argument('--dataroot', metavar='DIR', type=str, default=None, help='Root directory for data hierarchy')
    parser.add_argument('--logdir', metavar='DIR', type=str, default=None, help='Directory for log files')
    parser.add_argument('--logfile', type=str, default=None, help='Name of log file for running covers')
    parser.add_argument('--testbedcfg', metavar='CFGFILE', type=str, default='../nextstep.cfg', help='Name of general testbed configuration file (default is ../nextstep.cfg)')
    parser.add_argument('--numProcs', type=int, default=multiprocessing.cpu_count(), help='The number of processes to use when generating files.  By default it uses as many processes as their are CPU cores')
    parser.add_argument('--humanReadable', action='store_true', default=None, help='Use human readable mode for the datafile.  Otherwise a faster non-humanreadable output is used')

    # Which types of files should we be generating?
    # For now, do either reference, special-case refs, or plain/cover
    parser.add_argument('--refs', action='store_true', help='Generate reference IPD sequences')
    parser.add_argument('--refsspec', action='store_true', help='Generate "special" reference IPD sequences (overrides --refs)')
    
    # Number of IPDs to generate
    parser.add_argument('--covlens', metavar='L', type=int, nargs='*', default=None, help='Integer(s) giving the cover length(s) to use (read from config file unless specified here)')
    parser.add_argument('--reflens', metavar='L', type=int, default=None, help='Integer giving the number of reference IPDs to produce in each sequence')

    # Numbers of cover/reference files
    parser.add_argument('--numcovs', metavar='N', type=int, default=None, help='Number of cover sequences to produce.')
    parser.add_argument('--numrefs', metavar='N', type=int, default=None, help='Number of reference sequences to produce.')
    parser.add_argument('--minidx', metavar='IDX', type=int, default=None, help='Smallest index to generate (use if some files previously generated)')
    
    seedManagement = parser.add_argument_group(title='Seed Management')
    newSeedManagement = seedManagement.add_mutually_exclusive_group()
    newSeedManagement.add_argument('--newSeeds', default=True, action='store_true', help='Generate new seeds, --runID can be used to specify a custom name.')
    newSeedManagement.add_argument('--reuseSeeds', default=False, action='store_true', help='Use existing seeds, if set specify which set of seeds via --runID')
    seedManagement.add_argument('--runID', required='--reuseSeed' in sys.argv, type=str, default=None, help='The runID used to reffer to these seeds.  When generating new seeds if no runID is provided the current date and time will be used.  This argument is required when --reuseSeeds is set')
    seedManagement.add_argument('--bankPath', default=None, type=str, help='The path to the seed bank for the embedders.')

    # Indication of which distribution to use
    subparsers = parser.add_subparsers(help='Cover distributions sub-parsers.  Use their name followed by -h to find more info about them. These should be specified after the rest of your flags are set.')
    Pareto.generateSubparser(subparsers)
    ParetoNarrow.generateSubparser(subparsers)
    ParetoNarrow.generateSubparser(subparsers, 'paretoa100b10')
    ParetoSellke.generateSubparser(subparsers)
    ParetoSellke.generateSubparser(subparsers, 'paretoa100b095')
    Normal.generateSubparser(subparsers)
    MarkovIngest.generateSubparser(subparsers)

    args = parser.parse_args()
    if hasattr(args, 'parserID'):
        main(args)
    else:
        print('Please select a cover distribution.')
        parser.print_help()

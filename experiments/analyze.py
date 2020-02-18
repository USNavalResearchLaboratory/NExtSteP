#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#
import matplotlib
matplotlib.use('PDF')
from math import sqrt, log, log10
from os.path import isfile, isdir
from os import makedirs, path
from numpy import mean

import sys
sys.path.append('../') # this is done to find libStealthy
from libStealthy.util import rateF, aur, pe, pe_pair, calc_mssage_length, calc_mssage_prob, calc_rate_constants, DataFile, setArgsFromConfig
from libStealthy.boot import basicBootCImean, studentizedBootCImean

import matplotlib.pyplot as plt
import configparser
import datetime
import math
import json

embeddedRates = ['empty', 'const', 'log', 'cubrt', 'sqrt', '2rt3', 'linear']
detectStrs = ['welcht', 'mwu', 'filler', 'sgfiller', 'cce', 'chisq']


def compute_AURs(ifile, detectorStr, useStat=False):
    aurList = []
    df = DataFile(ifile)
    # detectorStr = df.getMetaData('detectorStr')
    outcomes = df.getData()
    for outcomes_for_ref in outcomes:
        if detectorStr == 'filler':
            # These classout files have (nu, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
            aurList.append(aur(couts))
        elif detectorStr == 'sgfiller':
            # These classout files have (nu, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
            aurList.append(aur(couts, unrev=True))
        elif detectorStr == 'cce' or detectorStr == 'localMinCCE':
            # These classout files have (entropy, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
            aurList.append(aur(couts))
        elif useStat:
            # These classout files have (stat, p, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[2] for x in outcomes_for_ref]
            aurList.append(aur(couts))
        else: # Other detectors, using p value
            # These classout files have (stat, p, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[1])
            couts = [x[2] for x in outcomes_for_ref]
            aurList.append(aur(couts))
    return aurList


def compute_AURs_unrev(ifile, detectorStr, useStat=False):
    aurList = []
    df = DataFile(ifile)
    # detectorStr = df.getMetaData('detectorStr')
    outcomes = df.getData()
    for outcomes_for_ref in outcomes:
        if detectorStr == 'filler':
            # These classout files have (nu, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
            aurList.append(aur(couts, unrev=True))
        elif detectorStr == 'sgfiller':
            # These classout files have (nu, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
            aurList.append(aur(couts, unrev=True))
        elif detectorStr == 'cce' or detectorStr == 'localMinCCE':
            # These classout files have (entropyDelta, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
            aurList.append(aur(couts, unrev=True))
        elif useStat:
            # These classout files have (stat, p, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[2] for x in outcomes_for_ref]
            aurList.append(aur(couts, unrev=True))
        else: # Other detectors, using p value
            # These classout files have (stat, p, type)
            # where type = 0 for plain and type = 1 for embedded
            outcomes_for_ref.sort(key = lambda x : x[1])
            couts = [x[2] for x in outcomes_for_ref]
            aurList.append(aur(couts, unrev=True))
    return aurList


def compute_PEs(ifile, detectorStr, useStat=False):
    '''
    Compute P_E values
    '''
    valList = []
    df = DataFile(ifile)
    # detectorStr = df.getMetaData('detectorStr')
    outcomes = df.getData()
    for outcomes_for_ref in outcomes:
        # outcomes_for_ref.sort(key = lambda x : x[1])
        if detectorStr in ['filler', 'sgfiller', 'cce', 'localMinCCE']:
            # useStat is irrelevant here
            # Classout file has (nu, type)
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
            valList.append(pe(couts))
        elif useStat:
            # Classout has (stat, p, type), using stat
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[2] for x in outcomes_for_ref]
            valList.append(pe(couts))
        else:
            # Classout has (stat, p, type), using p
            outcomes_for_ref.sort(key = lambda x : x[1])
            couts = [x[2] for x in outcomes_for_ref]
            valList.append(pe(couts))
    return valList


def compute_01PE(ifile, detectorStr, useStat=False):
    '''
    Compute P_E 01 values
    
    In this version, a 0 before a 1 counts as an error, i.e.,
    a perfect detector would produce something in 1*0*.
    '''
    valList = []
    df = DataFile(ifile)
    # detectorStr = df.getMetaData('detectorStr')
    outcomes = df.getData()
    for outcomes_for_ref in outcomes:
        # outcomes_for_ref.sort(key = lambda x : x[1])
        if detectorStr in ['filler', 'sgfiller', 'cce', 'localMinCCE']:
            # useStat is irrelevant here
            # Classout file has (nu, type)
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
        elif useStat:
            # Classout has (stat, p, type), using stat
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[2] for x in outcomes_for_ref]
        else:
            # Classout has (stat, p, type), using p
            outcomes_for_ref.sort(key = lambda x : x[1])
            couts = [x[2] for x in outcomes_for_ref]
        valList.append(pe_pair(couts)[0])

    return valList


def compute_10PE(ifile, detectorStr, useStat=False):
    '''
    Compute P_E 10 values
    
    In this version, a 1 before a 0 counts as an error, i.e.,
    a perfect detector would produce something in 0*1*.
    '''
    valList = []
    df = DataFile(ifile)
    # detectorStr = df.getMetaData('detectorStr')
    outcomes = df.getData()
    for outcomes_for_ref in outcomes:
        # outcomes_for_ref.sort(key = lambda x : x[1])
        if detectorStr in ['filler', 'sgfiller', 'cce', 'localMinCCE']:
            # useStat is irrelevant here
            # Classout file has (nu, type)
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[1] for x in outcomes_for_ref]
        elif useStat:
            # Classout has (stat, p, type), using stat
            outcomes_for_ref.sort(key = lambda x : x[0])
            couts = [x[2] for x in outcomes_for_ref]
        else:
            # Classout has (stat, p, type), using p
            outcomes_for_ref.sort(key = lambda x : x[1])
            couts = [x[2] for x in outcomes_for_ref]
        valList.append(pe_pair(couts)[1])

    return valList


def main(args):

    try:
        args = setArgsFromConfig(args.testbedcfg, args, 'analyze')
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
#     print('[###] Log file will be', args.logfile)
    print('[###] Anaylsis type will be', args.type)
#     print('[###] Data hierarchy rooted at', args.dataroot)
#     print('[###] Logs will be in', args.logdir)

    args.covlens = sorted(args.covlens)

    if not isdir(args.dataroot):
        raise RuntimeError('[!!!] Dataroot %s is not a directory!' % args.dataroot)
    if args.logdir is None:
        args.logdir = path.join(args.dataroot, 'logs')
    if not isdir(args.logdir):
         makedirs(args.logdir) 

    if 'ingest' in  args.dist:
        args.dataroot = path.join(args.dataroot, 'ingest/')
    else:
        args.dataroot = path.join(args.dataroot, 'hq/')

#     print('[###] Using', len(args.covlens), 'cover lengths:', args.covlens)

    numRefIPDs = args.reflen

    valDict = {}
    valDictMean = {}

    fmtStrs = ['ko-', 'k^:', 'ks--', 'kd-.', 'kx-', 'k>:']
    pointMark = ['o', '^', 's', 'd', 'x', '>']
    ecolors = ['gray', 'gray', 'gray', 'gray', 'gray', 'gray']

    maxIPD = args.maxipd

    cipercent = args.cipercent

    # rateList is here to allow multiple rates to be plotted on the
    # same plot.  That's not supported at this time.
    rateList = []
    if args.linear:
        rateList.append('linear')
    if args.rt2over3:
        rateList.append('2rt3')
    if args.sqrt:
        rateList.append('sqrt')
    if args.cubrt:
        rateList.append('cubrt')
    if args.log:
        rateList.append('log')
    if args.const:
        rateList.append('const')

    # ratio is set to empty for empty embedding
    # The empty embedding and empty rate have been replaced by
    # the donothing embedding.
    if args.embedding == 'empty':
        raise RuntimeError('empty not supported at this time!')
        rateList = ['empty']

    if len(rateList) != 1:
        raise RuntimeError('[!!!] Requiring a unique embedding rate!')

    detectorStr = args.detector

    k0 = args.rateK
    N0 = args.rateIPD
    rateC = calc_rate_constants(N0, k0)

    normStr = '-k' + str(k0) + '-n' + str(N0)

    if args.embedding == 'replay':
        embeddingTypeStr = 'replay'
    elif args.embedding == 'replayperc':
        embeddingTypeStr = 'replayperc-Q' + str(int(args.replayQ[0]))
    elif args.embedding == 'donothing':
        embeddingTypeStr = 'donothing'
    elif args.embedding =='replayprob':
        embeddingTypeStr ='replayprob-Q' + str(int(args.replayQ[0]))
    elif args.embedding == 'replayprob_alphabet':
        qstring=''
        for q in args.replayQ:
            qstring=qstring+str(int(q))+'-'
        embeddingTypeStr='replayprob-Q'+qstring
    elif args.embedding=='replayrepeat':
        embeddingTypeStr ='replayrepeat-' + str(int(args.numrepeat))
    elif args.embedding == 'sellke8to3':
        embeddingTypeStr = 'sellke-L8-n3'
    elif args.embedding == 'empty':
        raise RuntimeError('empty not supported at this time!')
        embeddingTypeStr = 'empty'

    if args.dist == 'sellkepareto' or args.dist == 'paretoa100b095':
        distParamStr = 'sellkepareto'
    elif args.dist == 'narrow' or args.dist == 'paretoa100b10':
        distParamStr = 'pareto-a100-b10'
    elif args.dist == 'pareto':
        alpha = args.paretoa
        beta = args.paretob
        if beta <= 0.0:
            raise ValueError('Beta must be positive!')
        distParamStr = 'pareto-a%d-b%s' % (alpha, str(beta).replace('.', 'p'))
    elif args.dist == 'normal':
        distmean = args.mean
        sd = args.sd
        if sd <= 0.0:
            raise ValueError('Std. dev. must be positive!')
        distParamStr = 'normal-m%d-s%s' % (distmean, str(sd).replace('.', 'p'))
    elif args.dist == 'ingest':
        distParamStr = 'ingest'
    
    elif args.dist == 'markov-ingest':
        distParamStr = 'markov-ingest'


    if embeddingTypeStr == 'donothing':
        idir_for_classout = path.join(args.dataroot, 'classout/' + embeddingTypeStr + '/')
    else:
        idir_for_classout = path.join(args.dataroot, 'classout/' + embeddingTypeStr + normStr + '/')

    fig_odir = path.join(args.dataroot, 'analyzeout/plots%s/%s/' % (normStr, args.type))

    ifile_cce_infix = ''
    ofile_cce_params = ''
    if detectorStr in ['cce','localMinCCE'] and len(args.cceBins) == 1:
        ifile_cce_infix = '-d' + str(args.cceSubSeqLen) + '-n' + str(args.cceBins[0])
        ofile_cce_params = '-d' + str(args.cceSubSeqLen) + '-n' + str(args.cceBins[0])
    elif detectorStr in ['cce','localMinCCE'] and type(args.cceBins) is int:
        ifile_cce_infix = '-d' + str(args.cceSubSeqLen) + '-n' + str(args.cceBins)
        ofile_cce_params = '-d' + str(args.cceSubSeqLen) + '-n' + str(args.cceBins)
    elif detectorStr in ['cce','localMinCCE']:
        ifile_cce_infix = '-d' + str(args.cceSubSeqLen) + '-cust'
        ofile_cce_params = '-d' + str(args.cceSubSeqLen) + '-cust'
    if args.useStat and detectorStr in ['mwu', 'welcht']:
        fig_odir += 'useStat/'

    ifile_chisq_infix = ''
    ofile_chisq_params = ''
    if detectorStr in ['chisq']:
        ifile_chisq_infix = "-numBins" + str(args.chisqBins)
        ofile_chisq_params = "-numBins" + str(args.chisqBins)

    if not isdir(fig_odir):
        makedirs(fig_odir)

    badKeys = set()

    if not isfile(path.join(args.logdir, args.logfile)):
        with open(path.join(args.logdir, args.logfile), 'w') as logstream:
            logstream.write('\n')

    with open(path.join(args.logdir, args.logfile), 'a') as logstream:
        logstream.write(str(datetime.datetime.now()) + '\n')
        summ_str = 'Doing %s analysis for %s distribution and %s embedding and %s classifier.\n' % (args.type, distParamStr, embeddingTypeStr, detectorStr)
        logstream.write(summ_str)

        for rateStr in rateList:
            if args.embedding == 'empty' or rateStr == 'empty':
                raise RuntimeError('empty not supported at this time!')
                MLs = [0 for cl in args.covlens]
            elif args.embedding == 'replay' or args.embedding == 'replayperc' or args.embedding == 'donothing':
                if args.embedding == 'donothing' and rateStr != 'sqrt':
                    raise RuntimeError('[!!!] Do-nothing embedding requires sqrt rate!')
                MLs = [calc_mssage_length(cl, rateStr, rateC[rateStr]) for cl in args.covlens]
                print('[###] CLs:', args.covlens)
                print('[###] MLs:' , MLs)
            elif args.embedding == 'sellke8to3':
                MLs = [calc_mssage_length(cl, rateStr, rateC[rateStr]) for cl in args.covlens]
            elif args.embedding =='replayprob' or args.embedding=='replayprob_alphabet':
#                 embedding_prob=lambda N: (float(k0)/N0*rateF[rateStr](N0))/rateF[rateStr](N)
                embedding_prob=lambda N: calc_mssage_prob(N, rateStr, rateC[rateStr])

                Probs=[embedding_prob(cl) for cl in args.covlens]
            for i in range(len(args.covlens)):
                cl = args.covlens[i]
                if embeddingTypeStr == 'empty':
                    raise RuntimeError('empty not supported at this time!')
                    ifile_prefix = detectorStr + ifile_chisq_infix + ifile_cce_infix + '-' + distParamStr + '-r' + str(numRefIPDs) + '-c'+ str(cl)
                elif args.embedding== 'replayprob':
                    prob=Probs[i]
                    ifile_prefix = detectorStr + ifile_chisq_infix + ifile_cce_infix + '-' + distParamStr + '-r' + str(numRefIPDs) + '-c'+ str(cl)+'-p'+str(int(10**5*prob))
                elif args.embedding=='replayprob_alphabet':
                    prob=Probs[i]
                    ifile_prefix = detectorStr + ifile_chisq_infix + ifile_cce_infix + '-' + distParamStr + '-r' + str(numRefIPDs) + '-c'+ str(cl)+'-p'+str(int(10**5*prob))
                elif args.embedding=='replayrepeat':
                    ifile_prefix = detectorStr + ifile_chisq_infix + ifile_cce_infix + '-' + distParamStr + '-r' + str(numRefIPDs) + '-c'+ str(cl)
                else:
                    ml = MLs[i]
                    ifile_prefix = detectorStr + ifile_chisq_infix + ifile_cce_infix + '-' + distParamStr + '-r' + str(numRefIPDs) + '-c'+ str(cl) + '-m' + str(ml)
                ifile = idir_for_classout + ifile_prefix + '.json'
#                 print('[###] Opening', ifile)

                if not isfile(ifile):
                    logstream.write('FILE: Unable to find ifile: ' + ifile + '\n')
                    print('[~~~] Unable to find ifile:', ifile)
                    badKeys.add((rateStr, detectorStr, cl))
                    continue

                if args.type == 'aur':
#                     print('[###] Doing AUR analysis!')
                    valList = compute_AURs(ifile, detectorStr, args.useStat)
#                     print('Done with AUR analysis!')
                elif args.type == 'aurunrev':
#                     print('[###] Doing unreversed AUR analysis!')
                    valList = compute_AURs_unrev(ifile, detectorStr, args.useStat)
                elif args.type == 'pe':
#                     print('[###] Doing P_E analysis!')
                    valList = compute_PEs(ifile, detectorStr, args.useStat)
                elif args.type == 'pe01':
#                     print('[###] Doing 0--1 P_E analysis!')
                    valList = compute_01PE(ifile, detectorStr, args.useStat)
                elif args.type == 'pe10':
#                     print('[###] Doing 0--1 P_E analysis!')
                    valList = compute_10PE(ifile, detectorStr, args.useStat)
                valDict[(rateStr, detectorStr, cl)] = valList[:]
                valDictMean[(rateStr, detectorStr, cl)] = mean(valList)

#             print('===========')
        print('=======================')
        logstream.write('\n\n')
    ### End of stuff that might be logged

    # Set up a plot
    # For each detectorStr, extract x and y values from aurDictMean
    # Plot lines for all detectorStrs

    _, ax = plt.subplots()
    anyPlot = False

    if args.type == 'aur':
        ylabelStr = 'Mean AUR'
    if args.type == 'aurunrev':
        ylabelStr = 'Mean Unreversed AUR'
    elif args.type == 'pe':
        ylabelStr = 'Mean 1 - P_E'
    elif args.type == 'pe01':
        ylabelStr = 'Mean 1 - P_E_01'
    elif args.type == 'pe10':
        ylabelStr = 'Mean 1 - P_E_10'

    if args.plotstyle == 'meanwhisker':
        if args.cimethod == 'basic':
            ylabelStr += ' (bars give basic-bootstrap ' + str(cipercent) + '% CIs)'
        elif args.cimethod == 'minmax':
            ylabelStr += ' (bars give min. and max.)'
        elif args.cimethod == 'none':
            pass
        elif args.cimethod == 'studentized':
            ylabelStr += ' (bars give studentized-bootstrap ' + str(cipercent) + '% CIs)'
    elif args.plotstyle == 'meanvals':
        ylabelStr += ' (dots show all values)'

    for j in range(len(rateList)):
        rateStr = rateList[j]
        NsToPlot = []
        valsToPlot = []
        valLists = []
        valLowerCI = []
        valUpperCI = []


        for i in range(len(args.covlens)):
            cl = args.covlens[i]
            if (rateStr, detectorStr, cl) in badKeys:
                continue

            NsToPlot.append(cl)
            valMean = valDictMean[(rateStr, detectorStr, cl)]
            valsToPlot.append(valMean)
            valLists.append(valDict[(rateStr, detectorStr, cl)])

            if args.cimethod == 'minmax':
                ### Set up the error bars to show the entire range of values
                valLowerCI.append(valMean - sorted(valDict[(rateStr, detectorStr, cl)])[0])
                valUpperCI.append(sorted(valDict[(rateStr, detectorStr, cl)])[-1] - valMean)

            elif args.cimethod == 'basic':
                lowerCI, upperCI = basicBootCImean(valDict[(rateStr, detectorStr, cl)],alpha=(cipercent/100.0))
                valLowerCI.append(valMean - lowerCI)
                valUpperCI.append(upperCI - valMean)

            elif args.cimethod == 'studentized':
                lowerCI, upperCI = studentizedBootCImean(valDict[(rateStr, detectorStr, cl)],alpha=(cipercent/100.0))
                valLowerCI.append(valMean - lowerCI)
                valUpperCI.append(upperCI - valMean)

            elif args.cimethod == 'none':
                # We shouldn't be using these values anyway
                valLowerCI.append(None)
                valUpperCI.append(None)

        seqLists = []
        if len(valLists) > 0 and len(valLists[0]) > 0:
            for i in range(len(valLists[0])):
                seq = []
                for k in range(len(valLists)):
                    seq.append(valLists[k][i])
                seqLists.append(seq)

        if args.plotstyle == 'meanwhisker': # Traditional (means with whiskers)
            if len(NsToPlot) == 0:
                continue
            anyPlot = True
            ax.errorbar(NsToPlot, valsToPlot, fmt=fmtStrs[j], label=rateStr, yerr=[valLowerCI, valUpperCI], ecolor=ecolors[j])
            ax.set_ylabel(ylabelStr)

        elif args.plotstyle == 'meanvals': # Dots plus mean
            if len(NsToPlot) == 0:
                continue
            anyPlot = True
            for k in range(len(NsToPlot)):
                xval = NsToPlot[k]
                yvals = valLists[k]
                ax.scatter([xval] * len(yvals), yvals, facecolors='none', marker=pointMark[j], edgecolors=ecolors[j], alpha=0.5)
            ax.set_ylabel(ylabelStr)
            ax.plot(NsToPlot, valsToPlot, fmtStrs[j], label=rateStr)

        elif args.plotstyle == 'allcurves': # Curves plus mean
            if len(NsToPlot) == 0:
                continue
            anyPlot = True
            for k in range(len(seqLists)):
                ax.plot(NsToPlot, seqLists[k], alpha=0.5)
            ax.set_ylabel(ylabelStr)
            ax.plot(NsToPlot, valsToPlot, fmtStrs[j], label=rateStr)

    if anyPlot:
        if args.type == 'aurunrev':
            ax.set_ylim([0.0, 1.01])
        else:
            ax.set_ylim([.5, 1.01])
        ax.set_xlim([0,maxIPD])
        ax.set_xlabel('Length of cover IPD sequence')

        titleStr = 'Results for '
        if embeddingTypeStr == 'replay':
            titleStr += 'replay embedding'
        elif args.embedding== 'replayperc':
            titleStr += 'replay perc. emb. (Q = ' + str(int(args.replayQ[0])) + ')'
        elif args.embedding== 'donothing':
            titleStr += 'do-nothing embedding'
        elif args.embedding == 'replayprob':
            titleStr += 'replay prob. emb. (Q = ' + str(int(args.replayQ[0])) + ')'
        elif args.embedding =='replayprob_alphabet':
            qstring=''
            for q in args.replayQ:
                qstring=qstring+str(int(q))+'-'
            titleStr+= 'replay prob. emb. (Q= '+ str(args.replayQ)+ ')'
        elif args.embedding == 'replayrepeat':
            titleStr += 'replay repeat (length repeated seq = ' + str(int(args.numrepeat)) + ')'
        elif embeddingTypeStr == 'sellke-L8-n3':
            titleStr += '8-to-3 embedding'
        elif embeddingTypeStr == 'empty':
            raise RuntimeError('empty not supported at this time!')
            titleStr += 'empty embedding'
        else:
            print('[~~~] embeddingTypeStr unmatched:', embeddingTypeStr)

        if embeddingTypeStr != 'empty':
            if distParamStr == 'sellkepareto':
                titleStr += ' into Pareto (a100, b0.95)'
            elif distParamStr == 'pareto-a100-b10':
                titleStr += ' into Pareto (a100, b10)'

        titleStr += ' \n with ' + detectorStr + ofile_chisq_params + ' classifier using reflen ' + str(numRefIPDs)

        plt.title(titleStr)

        fig_ofile = fig_odir + detectorStr + ofile_chisq_params + ofile_cce_params + '-' + embeddingTypeStr + '-' + distParamStr + '-' + args.type + '-' + args.plotstyle + '-' + rateList[0] + '-r' + str(numRefIPDs) + '-current.pdf'
        print('[###] Writing figure to:', fig_ofile)
        plt.savefig(fig_ofile)

        #write mean values (one per cover length) to files
        mean_ofile= fig_odir + detectorStr + '-' + embeddingTypeStr + '-' + distParamStr + '-mean_' + args.type + '-' + rateList[0] + '-r' + str(numRefIPDs) + '.json'
        with open(mean_ofile, 'w') as ostream:
            out = {
                'args.covlens':args.covlens
                ,'vals':valsToPlot
            }
            json.dump(out, ostream)
            
    else:
        print('[###] Nothing to plot!')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Produce analysis plots (multiple embedding ratios per plot) for studies of possible square-root laws')

    parser.add_argument('--dataroot', metavar='DIR', type=str, default=None, help='Root directory for data hierarchy')
    parser.add_argument('--outroot', metavar='DIR', type=str, default='./', help='Root directory for output hierarchy for writing plots, etc. (default is ./)')
    parser.add_argument('--logdir', metavar='DIR', type=str, default=None, help='Directory for log files')
    parser.add_argument('--logfile', type=str, default=None, help='Name of log file for analysis')
    parser.add_argument('--testbedcfg', metavar='CFGFILE', type=str, default='../nextstep.cfg', help='Name of general testbed configuration file (default is ../nextstep.cfg)')

    # Flags for different types of analysis plots
    parser.add_argument('--type', type=str, choices=['aur', 'pe', 'pe01', 'pe10', 'aurunrev'], default='aur', help='Type of analysis to run (AUR [default], 1-P_E, 1-P_E_01, 1-P_E_10, or unreversed AUR (nonstandard)')

    # Plot types
    parser.add_argument('--plotstyle', type=str, choices=['meanwhisker', 'meanvals', 'allcurves'], default='meanwhisker', help='Type of aur plot to show (mean with whiskers [default], mean with all values, or all---individual and mean---curves)')
    parser.add_argument('--cimethod', type=str, choices=['basic', 'minmax', 'none', 'studentized'], default='basic', help='Method for generating confidence intervals; basic bootstrapping (default; tested with data), min./max. values, none, and studentized bootstrapping (not tested with data)')
    parser.add_argument('--cipercent', type=float, default=95, help='Float percentage for confidence interval (default is 95)')

    # Options for presentation of analysis
    parser.add_argument('--maxipd', type=int, help='Largest IPD value to include in horizontal axis of analysis plots.', default=None)

    # Flags for different possible ratios
#     parser.add_argument('--empty', action='store_true', help='Include empty ratio')
# # --empty is unsupported by other parts of the testbed
# # Code is still present here, but ensure it's not used by removing the flag.
    parser.add_argument('--const', action='store_true', help='Include constant ratio')
    parser.add_argument('--log', action='store_true', help='Include logarithmic ratio')

    parser.add_argument('--cubrt', action='store_true', help='Include cubrt ratio')
    parser.add_argument('--sqrt', action='store_true', help='Include sqrt ratio')
    parser.add_argument('--rt2over3', action='store_true', help='Include two-thirds-root ratio')
    parser.add_argument('--linear', action='store_true', help='Include linear ratio')

    parser.add_argument('--rateIPD', type=float, metavar='N', default=None, help='Num. cover IPDs at which embedding rates are equalized')
    parser.add_argument('--rateK', type=float, metavar='k', default=None, help='Num. bits that should be embedded into <rateIPD> IPDs')

    # Indication of which distribution to use
    parser.add_argument('dist', type=str, help='Cover distribution to use (sellkepareto and narrow are aliases for paretoa100b095 and paretoa100b10, respectively)', choices=['pareto', 'sellkepareto', 'narrow', 'paretoa100b10', 'paretoa100b095', 'normal', 'ingest', 'markov-ingest'])

    # More about the distribution
    parser.add_argument('--paretoa', metavar='ALPHA', type=int, default=None, help='Alpha value for Pareto distribution (int)')
    # Process beta as a string to facilitate standardization of floating-point argument
    parser.add_argument('--paretob', metavar='BETA', type=float, default=None, help='Beta value for Pareto distribution (float)')
    parser.add_argument('--mean', metavar='MU', type=int, help='Mean value for normal distribution')
    parser.add_argument('--sd', metavar='SIGMA', type=float, help='Standard deviation for normal distribution')

    # Indication of which classifier to use
    parser.add_argument('detector', type=str, help='Classifier to use', choices=['welcht', 'mwu', 'filler', 'sgfiller', 'cce', 'localMinCCE', 'chisq'])
    parser.add_argument('--useStat', action='store_true', help='Sort by test statistic instead of p value')


    # CCE specific options
    CCEOptions = parser.add_argument_group(title='CCE Options')
    CCEOptions.add_argument('--cceBins', nargs='*', default = [5], help = "If a single int is supplied then it will be used as the number of bins.  If a list of numbers are supplied they will be used as the edges of the bins.  The bin count is also the branching factor of the tree.")
    CCEOptions.add_argument('--cceSubSeqLen', type=int, default = 50, help = 'The length of the sub-sequences made from the original sequence of IPDs')

    # Chi Square specific options
    ChiSqOptions = parser.add_argument_group(title='Chisq options')
    ChiSqOptions.add_argument('--chisqBins', type=int, default = 20, help = 'The number of bins to sort refIPDs into using Chisq analysis')

    # Indication of which embedding method to consider
    parser.add_argument('embedding', type=str, choices=['replay', 'replayperc', 'donothing', 'sellke8to3', 'replayprob', 'replayprob_alphabet', 'replayrepeat'], help='Embedding method to consider')

    # Details of embeddings
    parser.add_argument('--replayQ', metavar='Q', type=int, nargs='*', default=[90], help='Percentile (integer) to use for skewed replay embeddings (replayperc embedding; ignored if a different embedding is specified)')
    parser.add_argument('--numrepeat', metavar='rep', type=int, default=5, help='replayrepeat only; (ignored if a different embedding is specified) default 5')


    # Cover lengths to include
    parser.add_argument('--covlens', metavar='L', type=int, nargs='*', default=None, help='Integer(s) giving the cover length(s) to use (read from config file unless specified here)')

    # Number of reference IPDs
    parser.add_argument('--reflen', type=int, default=None, help='Number of reference IPDs to use')



    args = parser.parse_args()
    main(args)

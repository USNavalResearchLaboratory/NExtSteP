#!/bin/bash
#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#

# The number and length of the reference/plain/cover sequences
# may need to be adjusted based on how much traffic is available.
# This can be done on the command line or via the .cfg file.
# It seems useful to have rateIPD be no larger than the smallest
# cover length used.
# Similarly, the maxipd value for the analyze.py script may need
# to be much smaller for ingested data than for synthetic data.

logdir='../../data/logs'
logfile='ingest.log'

# The inputfile should be a JSON-encoded list of IPD values
inputfile='../../data/inputdata/nextstep-input.json'

echo "`date`: ------- INGEST DATA -------" >> $logdir/$logfile
echo "`date`: ------- INGEST DATA -------"
echo "`date`: genconvers perc $inputfile" >> $logdir/$logfile
echo "`date`: genconvers perc $inputfile"
# --overrideData not supported in this release
### if you are starting for the first time and the data directory does not exist yet, then --overrideData no
# python ingest-data.py perc --overrideData yes $inputfile
python ingest-data.py perc $inputfile
echo

#######

cd ../experiments/

########################### Do-nothing embedding ####################################

echo "`date`: ------- EMBED (INGEST) -------" >> $logdir/$logfile
echo "`date`: ------- EMBED (INGEST) -------"
for rate in 'sqrt'; do
	for emb in 'donothing'; do
		echo "`date`: embed $rate --dist ingest --$emb" >> $logdir/$logfile
		echo "`date`: embed $rate --dist ingest --$emb"
		python embed.py $rate --dist ingest --$emb
	done
done
echo

echo "`date`: ------- CLASSIFY (MWU, FILLER, WELCHT) -------" >> $logdir/$logfile
echo "`date`: ------- CLASSIFY (MWU, FILLER, WELCHT) -------"
for rate in 'sqrt'; do
	for det in 'mwu' 'welcht' 'filler'; do
		for emb in 'donothing'; do
			echo "`date`: classify ingest $rate $emb --$det" >> $logdir/$logfile
			echo "`date`: classify ingest $rate $emb --$det" 
			python classify.py ingest $rate $emb --$det
		done
	done
done
echo

echo "`date`: ------- ANALYZE (MWU, WELCHT, FILLER) -------" >> $logdir/$logfile
echo "`date`: ------- ANALYZE (MWU, WELCHT, FILLER) -------"
for rate in 'sqrt'; do
	for emb in 'donothing'; do
		for det in 'welcht' 'mwu' 'filler'; do
			echo "`date`: analyze ingest $det $emb --$rate --plotstyle allcurves" >> $logdir/$logfile
			echo "`date`: analyze ingest $det $emb --$rate --plotstyle allcurves" 
			python analyze.py ingest $det $emb --$rate --plotstyle allcurves	
		done
	done
done

echo "`date`: ------- CLASSIFY (FULL CCE) -------" >> $logdir/$logfile
echo "`date`: ------- CLASSIFY (FULL CCE) -------"
for rate in 'sqrt'; do
	for emb in 'donothing'; do
#         for depth in 5 15 50; do
		for depth in 15; do
			echo "`date`: classify ingest $rate $emb --cce --cceSubSeqLen $depth --cceBins 5" >> $logdir/$logfile
			echo "`date`: classify ingest $rate $emb --cce --cceSubSeqLen $depth --cceBins 5"
			python classify.py ingest $rate $emb --cce --cceSubSeqLen $depth --cceBins 5
		done
		
# 		for depth in 5 10 25; do
		for depth in 15; do
			for nbins in 3 10; do
				echo "`date`: classify ingest $rate $emb --cce --cceSubSeqLen $depth --cceBins $nbins" >> $logdir/$logfile
				echo "`date`: classify ingest $rate $emb --cce --cceSubSeqLen $depth --cceBins $nbins" 
				python classify.py ingest $rate $emb --cce --cceSubSeqLen $depth --cceBins $nbins
			done
		done
	done
done
echo

echo "`date`: ------- ANALYZE (FULL CCE) -------" >> $logdir/$logfile
echo "`date`: ------- ANALYZE (FULL CCE) -------"
for rate in 'sqrt'; do
	for emb in 'donothing'; do
		for det in 'cce'; do
# 			for depth in 5 15 50; do
			for depth in 15; do
				echo "`date`: analyze ingest $det $emb --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins 5" >> $logdir/$logfile
				echo "`date`: analyze ingest $det $emb --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins 5" 
				python analyze.py ingest $det $emb --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins 5
			done
		
# 			for depth in 5 10 25; do
			for depth in 15; do
				for nbins in 3 10; do
					echo "`date`: analyze ingest $det $emb --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins $nbins" >> $logdir/$logfile
					echo "`date`: analyze ingest $det $emb --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins $nbins" 
					python analyze.py ingest $det $emb --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins $nbins
				done
			done
		done
	done
done

echo "`date`: ------- CLASSIFY (CHISQUARED) -------" >> $logdir/$logfile
echo "`date`: ------- CLASSIFY (CHISQUARED) -------"
for rate in 'sqrt'; do
	for emb in 'donothing'; do
		for numBins in 10 25 100; do
			echo "`date`: classify ingest $rate $emb --chisq --chisqBins $numBins" >> $logdir/$logfile
			echo "`date`: classify ingest $rate $emb --chisq --chisqBins $numBins" 
			python classify.py ingest $rate $emb --chisq --chisqBins $numBins
		done
	done
done
echo


echo "`date`: ------- ANALYZE (CHISQUARED) -------" >> $logdir/$logfile
echo "`date`: ------- ANALYZE (CHISQUARED) -------"
for rate in 'sqrt'; do
	for emb in 'donothing'; do
		for numBins in 10 25 100; do
			echo "`date`: analyze ingest chisq $emb --$rate --plotstyle allcurves --chisqBins $numBins" >> $logdir/$logfile
			echo "`date`: analyze ingest chisq $emb --$rate --plotstyle allcurves --chisqBins $numBins" 
			python analyze.py ingest chisq $emb --$rate --plotstyle allcurves --chisqBins $numBins
		done
	done
done

################################################## With Embeddings #############################

echo "`date`: ------- EMBED (INGEST) -------" >> $logdir/$logfile
echo "`date`: ------- EMBED (INGEST) -------"
for rate in 'const' 'linear' 'sqrt' 'cubrt' '2rt3'; do
	for emb in 'replayperc' 'replayprob'; do
		for Q in 85 50; do
			echo "`date`: embed $rate --dist ingest --$emb --replayQ $Q" >> $logdir/$logfile
			echo "`date`: embed $rate --dist ingest --$emb --replayQ $Q"
			python embed.py $rate --dist ingest --$emb --replayQ $Q
		done
	done
done
echo

echo "`date`: ------- CLASSIFY (MWU, FILLER, WELCHT) -------" >> $logdir/$logfile
echo "`date`: ------- CLASSIFY (MWU, FILLER, WELCHT) -------"
for rate in 'const' 'linear' 'sqrt' 'cubrt' '2rt3'; do
	for Q in 85 50; do
		for det in 'mwu' 'welcht' 'filler'; do
			for emb in 'replayperc' 'replayprob'; do
				echo "`date`: classify ingest $rate $emb --replayQ $Q --$det" >> $logdir/$logfile
				echo "`date`: classify ingest $rate $emb --replayQ $Q --$det" 
				python classify.py ingest $rate $emb --replayQ $Q --$det
			done
		done
	done
done
echo

echo "`date`: ------- ANALYZE (MWU, WELCHT, FILLER) -------" >> $logdir/$logfile
echo "`date`: ------- ANALYZE (MWU, WELCHT, FILLER) -------"
for rate in 'const' 'linear' 'sqrt' 'cubrt' 'rt2over3'; do
	for Q in 85 50; do
		for emb in 'replayperc' 'replayprob'; do
			for det in 'welcht' 'mwu' 'filler'; do
				echo "`date`: analyze ingest $det $emb --replayQ $Q --$rate --plotstyle allcurves" >> $logdir/$logfile
				echo "`date`: analyze ingest $det $emb --replayQ $Q --$rate --plotstyle allcurves" 
				python analyze.py ingest $det $emb --replayQ $Q --$rate --plotstyle allcurves
			done
		done
	done
done

echo "`date`: ------- CLASSIFY (FULL CCE) -------" >> $logdir/$logfile
echo "`date`: ------- CLASSIFY (FULL CCE) -------"
for rate in 'const' 'linear' 'sqrt' 'cubrt' '2rt3'; do
	for Q in 85 50; do
		for emb in 'replayperc' 'replayprob'; do
# 			for depth in 5 15 50; do
			for depth in 15; do
				echo "`date`: classify ingest $rate $emb --replayQ $Q --cce --cceSubSeqLen $depth --cceBins 5" >> $logdir/$logfile
				echo "`date`: classify ingest $rate $emb --replayQ $Q --cce --cceSubSeqLen $depth --cceBins 5"
				python classify.py ingest $rate $emb --replayQ $Q --cce --cceSubSeqLen $depth --cceBins 5
			done
			
			for depth in 15; do
				for nbins in 3 10; do
					echo "`date`: classify ingest $rate $emb --replayQ $Q --cce --cceSubSeqLen $depth --cceBins $nbins" >> $logdir/$logfile
					echo "`date`: classify ingest $rate $emb --replayQ $Q --cce --cceSubSeqLen $depth --cceBins $nbins" 
					python classify.py ingest $rate $emb --replayQ $Q --cce --cceSubSeqLen $depth --cceBins $nbins
				done
			done
		done
	done
done
echo

echo "`date`: ------- ANALYZE (FULL CCE) -------" >> $logdir/$logfile
echo "`date`: ------- ANALYZE (FULL CCE) -------"
for rate in 'const' 'linear' 'sqrt' 'cubrt' 'rt2over3'; do
	for Q in 85 50; do
		for emb in 'replayperc' 'replayprob'; do
			for det in 'cce'; do
# 				for depth in 5 15 50; do
				for depth in 15; do
					echo "`date`: analyze ingest $det  $emb --replayQ $Q --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins 5" >> $logdir/$logfile
					echo "`date`: analyze ingest $det  $emb --replayQ $Q --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins 5" 
					python analyze.py ingest $det  $emb --replayQ $Q --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins 5
				done
			
# 				for depth in 5 10 25; do
				for depth in 15; do
					for nbins in 3 10; do
						echo "`date`: analyze ingest $det  $emb --replayQ $Q --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins $nbins" >> $logdir/$logfile
						echo "`date`: analyze ingest $det  $emb --replayQ $Q --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins $nbins" 
						python analyze.py ingest $det  $emb --replayQ $Q --$rate --plotstyle allcurves --cceSubSeqLen $depth --cceBins $nbins
					done
				done
			done
		done
	done
done

echo "`date`: ------- CLASSIFY (CHISQUARED) -------" >> $logdir/$logfile
echo "`date`: ------- CLASSIFY (CHISQUARED) -------"
for rate in 'const' 'linear' 'sqrt' 'cubrt' '2rt3'; do
	for Q in 85 50; do
		for emb in 'replayperc' 'replayprob'; do
			for numBins in 10 25 100; do
				echo "`date`: classify ingest $rate $emb --replayQ $Q --chisq --chisqBins $numBins" >> $logdir/$logfile
				echo "`date`: classify ingest $rate $emb --replayQ $Q --chisq --chisqBins $numBins" 
				python classify.py ingest $rate $emb --replayQ $Q --chisq --chisqBins $numBins
			done
		done
	done
done
echo

echo "`date`: ------- ANALYZE (CHISQUARED) -------" >> $logdir/$logfile
echo "`date`: ------- ANALYZE (CHISQUARED) -------"
for rate in 'const' 'linear' 'sqrt' 'cubrt' 'rt2over3'; do
	for Q in 85 50; do
		for emb in 'replayperc' 'replayprob'; do
			for numBins in 10 25 100; do
				echo "`date`: analyze ingest chisq $emb --replayQ $Q --$rate --plotstyle allcurves --chisqBins $numBins" >> $logdir/$logfile
				echo "`date`: analyze ingest chisq $emb --replayQ $Q --$rate --plotstyle allcurves --chisqBins $numBins" 
				python analyze.py ingest chisq $emb --replayQ $Q --$rate --plotstyle allcurves --chisqBins $numBins
			done
		done
	done
done


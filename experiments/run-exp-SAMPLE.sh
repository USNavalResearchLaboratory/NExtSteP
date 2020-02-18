#!/bin/bash
#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#

# Do embeddings

# Here, just do one rate constant
# Use rK = 200 for Q = 85; rK = 3200 for Q = 50
for rK in 200; do
    # Iterate through embedding rates
    for rate in 'linear' 'cubrt' 'sqrt' '2rt3' 'const'; do
        echo "`date`: starting embedding for $rate"
        # Iterate through embedding algorithms
        for emb in 'replayperc' 'replayprob'; do
            # Iterate through Q parameter for replayperc embedding
            for Q in 85; do
                # Iterate through (sythetic) traffic distributions
                # Each distribution needs to have been used to generate
                # reference/plain/cover traffic
                for dist in 'narrow' 'normal'; do
                    echo "`date`: embed $dist $rate $emb Q$Q rK$rK"                        
                    echo "`date`: embed $dist $rate $emb Q$Q rK$rK" >> exp-script.log  
                    python embed.py $rate --testbedcfg ../nextstep.cfg --rateK $rK --dist $dist --$emb --replayQ $Q
                done
            done
        done
    done
done

# Now do classification and analysis

# All the parameters here need to have been covered in the embedding
# process above.  One could move the embedding into this loop if desired.
#
# Rate constant
for rK in 200; do
    # Traffic distribution
    for dist in 'narrow' 'normal'; do
        for emb in 'replayperc' 'replayprob'; do
            for Q in 85; do
                for rate in 'linear' '2rt3' 'sqrt' 'cubrt' 'const'; do
                    # The analysis script uses 'rt2over3' instead of '2rt3'
                    # Issue with starting a flag with a numeral
                    if [ $rate == '2rt3' ]; then
                        analyzeRate='rt2over3'
                    else
                        analyzeRate=$rate
                    fi
                    echo "`date`: starting classification for $rate"
                    # Iterate over reference lengths
                    for reflen in 100000; do
                        # Run detectors.
                        # In each case, run corresponding analysis script to
                        # produce a plot as soon as the classify script finishes
                        # 
                        # Run the first set of detectors.
                        # These do not take any additional parameters
                        for det in 'welcht' 'mwu' 'filler'; do
                            echo "`date`: classify $dist $rate $emb Q$Q $det rK$rK rL$reflen"
                            echo "`date`: classify $dist $rate $emb Q$Q $det rK$rK rL$reflen" >> exp-script.log  
                            python classify.py $dist $rate $emb --$det --testbedcfg ../nextstep.cfg --replayQ $Q --rateK $rK --reflen $reflen

                            # Iterate through metrics on detectors
                            # Each produces a different plot
                            for type in 'aur' 'aurunrev' 'pe' 'pe01' 'pe10'; do
                                echo "`date`: analyze $dist $rate $emb Q$Q rK$rK $det $type"
                                echo "`date`: analyze $dist $rate $emb Q$Q rK$rK $det $type" >> exp-script.log
                                # This shows the curves relative to each
                                # reference string plus a summary curve
                                #
                                # --useStat uses the test statistic instead
                                # of the p value, when the latter is computed,
                                # to sort the outputs.  This is ignored by tests
                                # that do not produce a p value.
                                python analyze.py --testbedcfg ../nextstep.cfg --reflen $reflen --rateK $rK $dist $det $emb --replayQ $Q --$analyzeRate --plotstyle allcurves --type $type --useStat
                            done

                        done

                        # Now do chi-squared detector, which also takes nbins as
                        # a parameter
                        for nbins in 15; do
                            echo "`date`: classify $dist $rate $emb Q$Q chisq rK$rK rL$reflen bins$nbins" >> exp-script.log
                            echo "`date`: classify $dist $rate $emb Q$Q chisq rK$rK rL$reflen bins$nbins" 
                            python classify.py --testbedcfg ../nextstep.cfg --chisq --reflen $reflen --rateK $rK $dist $rate $emb --replayQ $Q  --chisqBins $nbins

                            for type in 'aur' 'aurunrev' 'pe' 'pe01' 'pe10'; do
                                echo "`date`: analyze $dist $rate $emb Q$Q rK$rK chisq $type"
                                echo "`date`: analyze $dist $rate $emb Q$Q rK$rK chisq $type" >> exp-script.log
                                python analyze.py --testbedcfg ../nextstep.cfg --reflen $reflen --rateK $rK $dist chisq $emb --chisqBins $nbins --replayQ $Q --$analyzeRate --plotstyle allcurves --type $type --useStat
                            done

                        done

                        # Iterate through the CCE tree depth and branching
                        # factor (nbins).
                        #
                        # Using large values of both depth and nbins can be
                        # problematic.
                        for depth in 5; do
                            for nbins in 15; do
                                echo "`date`: classify rL$reflen rK$rK $dist $rate $emb Q$Q cceSSL$depth bins$nbins" >> exp-script.log
                                echo "`date`: classify rL$reflen rK$rK $dist $rate $emb Q$Q cceSSL$depth bins$nbins"
                                python classify.py --reflen $reflen --testbedcfg ../nextstep.cfg --rateK $rK $dist $rate $emb --replayQ $Q --cce --cceSubSeqLen $depth --cceBins $nbins

                                for type in 'aur' 'aurunrev' 'pe' 'pe01' 'pe10'; do
                                    echo "`date`: analyze $dist $rate $emb Q$Q rK$rK cce bins$nbins depth$depth $type"
                                    echo "`date`: analyze $dist $rate $emb Q$Q rK$rK cce bins$nbins depth$depth $type" >> exp-script.log
                                    python analyze.py --testbedcfg ../nextstep.cfg --reflen $reflen --rateK $rK $dist cce $emb --replayQ $Q --cceSubSeqLen $depth --cceBins $nbins --$analyzeRate --plotstyle allcurves --type $type --useStat
                                done
                            done
                        done

                    done
                done
            done
        done
    done
done

# Now the "donothing" embedder, which just produces "embedded" IPDs that
# are the same as the "cover" IPDs.  This requires the "sqrt" embedding rate
# to avoid generating lots of copies of the same data.
#
# "donothing" embedder is useful to get a reference for how much a detector
# will distinguish "cover" IPDs from "plain" IPDs, even before an embedding
# is done.

for rate in 'sqrt'; do
    echo "`date`: starting donothing embedding for $rate"
    for emb in 'donothing'; do
        for dist in 'narrow' 'normal'; do
            echo "`date`: embed $dist $rate $emb"                        
            echo "`date`: embed $dist $rate $emb" >> exp-script.log  
            python embed.py $rate --testbedcfg ../nextstep.cfg --dist $dist --$emb
        done
    done  
done

for dist in 'narrow' 'normal'; do
    for emb in 'donothing'; do
        for rate in 'sqrt'; do
            if [ $rate == '2rt3' ]; then
                analyzeRate='rt2over3'
            else
                analyzeRate=$rate
            fi
            echo "`date`: starting classification for $rate"
            for reflen in 100000; do

                for det in 'welcht' 'mwu' 'filler'; do
                    echo "`date`: classify $dist $rate $emb $det rL$reflen"
                    echo "`date`: classify $dist $rate $emb $det rL$reflen" >> exp-script.log  
                    python classify.py $dist $rate $emb --$det --testbedcfg ../nextstep.cfg --reflen $reflen

                    for type in 'aur' 'aurunrev' 'pe' 'pe01' 'pe10'; do
                        echo "`date`: analyze $dist $rate $emb $det $type"
                        echo "`date`: analyze $dist $rate $emb $det $type" >> exp-script.log
                        python analyze.py --testbedcfg ../nextstep.cfg --reflen $reflen $dist $det $emb --sqrt --plotstyle allcurves --type $type --useStat
                    done

                done

                for nbins in 15; do
                    echo "`date`: classify $dist $rate $emb Q$Q chisq rK$rK rL$reflen bins$nbins" >> exp-script.log
                    echo "`date`: classify $dist $rate $emb Q$Q chisq rK$rK rL$reflen bins$nbins" 
                    python classify.py --testbedcfg ../nextstep.cfg --chisq --reflen $reflen $dist $rate $emb --chisqBins $nbins

                    for type in 'aur' 'aurunrev' 'pe' 'pe01' 'pe10'; do
                        echo "`date`: analysze $dist $rate $emb Q$Q rK$rK chisq $type"
                        echo "`date`: analyze $dist $rate $emb Q$Q rK$rK chisq $type" >> exp-script.log
                        python analyze.py --testbedcfg ../nextstep.cfg --reflen $reflen $dist chisq $emb --chisqBins $nbins --sqrt --plotstyle allcurves --type $type --useStat
                    done

                done

                for depth in 5; do
                    for nbins in 15; do

                        echo "`date`: classify rL$reflen rK$rK $dist $rate $emb Q$Q cceSSL$depth bins$nbins" >> exp-script.log
                        echo "`date`: classify rL$reflen rK$rK $dist $rate $emb Q$Q cceSSL$depth bins$nbins"
                        python classify.py --reflen $reflen --testbedcfg ../nextstep.cfg $dist $rate $emb --cce --cceSubSeqLen $depth --cceBins $nbins

                        for type in 'aur' 'aurunrev' 'pe' 'pe01' 'pe10'; do
                            echo "`date`: analyze $dist $rate $emb Q$Q rK$rK cce bins$nbins depth$depth $type"
                            echo "`date`: analyze $dist $rate $emb Q$Q rK$rK cce bins$nbins depth$depth $type" >> exp-script.log
                            python analyze.py --testbedcfg ../nextstep.cfg --reflen $reflen $dist cce $emb --cceSubSeqLen $depth --cceBins $nbins --sqrt --plotstyle allcurves --type $type --useStat
                        done
                    done
                done

            done
        done
    done
done

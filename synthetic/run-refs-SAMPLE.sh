#!/bin/bash
#
# U.S. Naval Research Laboratory
# Center for High Assurance Computer Systems
#

# Generate references for both narrow Pareto and normal distributions
# Can vary reference lengths here
# The reference lengths that are used by the experimental scripts need
# to have been generated here.  Can also rely on the default value given by
# the .cfg file instead of explicitly specifying the reference length.
for rl in 100000; do
    python generate-data.py --testbedcfg ../nextstep.cfg --refs --reflens $rl narrow
    python generate-data.py --testbedcfg ../nextstep.cfg --refs --reflens $rl normal
done

# Generate cover/plain sequences for both narrow Pareto and normal distributions
python generate-data.py narrow
python generate-data.py normal

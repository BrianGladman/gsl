#! /bin/bash -xv
#
# Diehard is a battery of tests for random number generators.
# It is available from http://stat.fsu.edu/pub/diehard/
#
# Note that 'diehard' needs to read its file "operm.cov" in the current
# directory, so you'll need to copy that file here if you plan to run this
# script in the current directory.
#
# The symptom of not doing this is that the OPERM5 test fails completely.
#

#for RNG in mt19937 ranlux ranlux389 cmrg mrg taus tt800
for RNG in gfsr4
do
GSL_RNG_TYPE=$RNG ./rng-dump > tmp.$RNG
./diehard > tmp.$RNG.log <<EOF
tmp.$RNG
16
EOF
rm -f tmp.$RNG
done

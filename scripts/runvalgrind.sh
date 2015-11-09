#!/bin/sh
#
# Run all tests under valgrind/memcheck
#
# First, compile the library with:
# > ./configure --disable-shared CFLAGS="-g -Wall"
# > make ; make check
# to compile the library and tests without optimization

outfile="valgrind.out"
rm -f $outfile

# valgrind exit code to check for errors
vcode="100"
vflags="--leak-check=full --error-exitcode=${vcode}"

npass=0
nfail=0
ntot=0

for testprog in $(ls */test); do
  echo "Running valgrind on ${testprog}"
  valgrind ${vflags} ${testprog} >> $outfile 2>&1
  rc=$?
  if [[ $rc == ${vcode} ]]; then
    echo "FAIL: $testprog"
    nfail=$((nfail+1))
  else
    npass=$((npass+1))
  fi
  ntot=$((ntot+1))
done

echo "FAIL: ${nfail}/${ntot}"
echo "PASS: ${npass}/${ntot}"

echo "output file is $outfile"

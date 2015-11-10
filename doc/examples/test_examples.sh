#!/bin/sh

nfail=0
npass=0
ntot=0

function dotest
{
  prog=$1
  file_out=$2
  file_err=$3
  args="$4"

  tmpout=$(mktemp)
  tmperr=$(mktemp)

  echo "testing $prog"
  ntot=$((ntot+1))

  eval ./${prog} ${args} 1> $tmpout 2> $tmperr

  # test stdout output
  str=$(/bin/diff $tmpout $file_out)
  if [ -n "$str" ]; then
    echo "FAIL(stdout): $prog"
    echo "difference in $file_out:"
    echo $str
    nfail=$((nfail+1))
  elif [ -n "$file_err" ]; then
    # test stderr output
    str=$(/bin/diff $tmperr $file_err)
    if [ -n "$str" ]; then
      echo "FAIL(stderr): $prog"
      echo "difference in $file_err:"
      echo $str
      nfail=$((nfail+1))
    else
      npass=$((npass+1))
    fi
  elif [ -s $tmperr ]; then
    echo "FAIL(stderr): nonzero output but no file for comparison"
    nfail=$((nfail+1))
  else
    npass=$((npass+1))
  fi

  rm -f $tmpout
  rm -f $tmperr
}

dotest blas blas.txt "" ""
dotest bspline bspline.txt bspline.err ""
dotest cblas cblas.txt "" ""
dotest cdf cdf.txt "" ""
dotest cheb cheb.txt "" ""
dotest combination combination.txt "" ""
dotest const const.txt "" ""
dotest diff diff.txt "" ""
dotest dwt dwt.txt "" "ecg.dat"
dotest eigen eigen.txt "" ""
dotest eigen_nonsymm eigen_nonsymm.txt "" ""
dotest fft fft.txt "" ""
dotest fftmr fftmr.txt "" ""
dotest fftreal fftreal.txt "" ""
dotest fitreg fitreg.txt fitreg.err ""
dotest fitting fitting.txt "" ""
dotest fitting2 fitting2.txt "" "19 < exp.dat"
dotest histogram2d histogram2d.txt "" ""
dotest ieee ieee.txt "" ""
dotest ieeeround ieeeround.txt "" ""
dotest integration integration.txt "" ""
dotest interp interp.txt "" ""
dotest interp2d interp2d.txt "" ""
dotest interp_compare interp_compare.txt "" ""
dotest interpp interpp.txt "" ""
dotest intro intro.txt "" ""
dotest largefit largefit.txt largefit.err ""
dotest largefit largefit2.txt largefit2.err "1"
dotest linalglu linalglu.txt "" ""
dotest matrix matrix.txt matrix.err ""
dotest matrixw matrixw.txt "" ""

echo "FAIL: ${nfail}/${ntot}"
echo "PASS: ${npass}/${ntot}"

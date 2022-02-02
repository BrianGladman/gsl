#!/bin/sh

exdir="doc/examples"

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

  eval ${exdir}/${prog} ${args} 1> $tmpout 2> $tmperr

  # test stdout output
  str=$(/bin/diff $tmpout ${exdir}/${file_out})
  if [ -n "$str" ]; then
    echo "FAIL(stdout): $prog"
    echo "difference in $file_out:"
    echo $str
    nfail=$((nfail+1))
  elif [ -n "$file_err" ]; then
    # test stderr output
    str=$(/bin/diff $tmperr ${exdir}/${file_err})
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
dotest bspline_deriv bspline_deriv.txt "" ""
dotest bspline_gram bspline_gram1.txt "" ""
dotest bspline_gram bspline_gram2.txt "" "0 0.6"
dotest bspline_interp bspline_interp.txt "" ""
dotest bspline_lsbreak bspline_lsbreak.txt bspline_lsbreak.err ""
dotest bspline_lsend bspline_lsend.txt bspline_lsend.err ""
dotest bspline_lsorder bspline_lsorder.txt "" ""
dotest bspline_per bspline_per.txt bspline_per.err ""
dotest bspline_proj bspline_proj.txt "" ""
dotest cblas cblas.txt "" ""
dotest cdf cdf.txt "" ""
dotest cheb cheb.txt "" ""
dotest combination combination.txt "" ""
dotest const const.txt "" ""
dotest diff diff.txt "" ""
dotest dwt dwt.txt "" "${exdir}/ecg.dat"
dotest eigen eigen.txt "" ""
dotest eigen_nonsymm eigen_nonsymm.txt "" ""
dotest fft fft.txt "" ""
dotest fftmr fftmr.txt "" ""
dotest fftreal fftreal.txt "" ""
dotest filt_edge filt_edge.txt "" ""
dotest fitreg fitreg.txt fitreg.err ""
dotest fitreg2 fitreg2.txt fitreg2.err ""
dotest fitting fitting.txt "" ""
dotest fitting2 fitting2.txt "" "19 < ${exdir}/exp.dat"
dotest gaussfilt gaussfilt.txt "" ""
dotest gaussfilt2 gaussfilt2.txt "" ""
dotest histogram2d histogram2d.txt "" ""
dotest ieee ieee.txt "" ""
dotest ieeeround ieeeround.txt "" ""
dotest impulse impulse.txt "" ""
dotest integration integration.txt "" ""
dotest integration2 integration2a.txt "" "10 5"
dotest integration2 integration2b.txt "" "10 6"
dotest interp interp.txt "" ""
dotest interp2d interp2d.txt "" ""
dotest interp_compare interp_compare.txt "" ""
dotest interpp interpp.txt "" ""
dotest intro intro.txt "" ""
dotest largefit largefit.txt largefit.err ""
dotest largefit largefit2.txt largefit2.err "1e-5"
dotest linalglu linalglu.txt "" ""
#dotest matrix matrix.txt matrix.err ""
dotest matrixw matrixw.txt "" ""
dotest min min.txt "" ""
dotest monte monte.txt "" ""
dotest movstat1 movstat1.txt "" ""
dotest movstat2 movstat2.txt "" ""
dotest movstat3 movstat3.txt "" ""
dotest multiset multiset.txt "" ""
dotest nlfit nlfit.txt nlfit.err ""
dotest nlfit2 nlfit2.txt nlfit2.err ""
dotest nlfit2b nlfit2b.txt nlfit2b.err ""
dotest nlfit3 nlfit3.txt nlfit3.err ""
# Times will be different for nlfit4
#dotest nlfit4 nlfit4.txt nlfit4.err ""
dotest ode-initval ode-initval.txt "" ""
dotest permseq permseq.txt "" ""
dotest permshuffle permshuffle.txt "" ""
dotest poisson poisson.txt poisson.err ""
dotest polyroots polyroots.txt "" ""
dotest qrng qrng.txt "" ""
dotest randpoisson randpoisson.txt "" ""
dotest randwalk randwalk.txt "" ""
dotest rng rng.txt "" ""
dotest rngunif rngunif.txt "" ""
dotest robfit robfit.txt "" "100"
dotest rootnewt rootnewt.txt "" ""
dotest roots roots.txt "" ""
dotest rquantile rquantile.txt "" ""
dotest rstat rstat.txt "" ""
dotest siman siman.txt "" ""
dotest sortsmall sortsmall.txt "" ""
dotest specfun specfun.txt "" ""
dotest specfun_e specfun_e.txt "" ""
dotest specfun2 specfun2.txt "" ""
dotest spmatrix spmatrix.txt "" ""
dotest stat stat.txt "" ""
dotest statsort statsort.txt "" ""
dotest sum sum.txt "" ""
dotest vectorview vectorview.txt "" ""

export GSL_RNG_TYPE=mrg
export GSL_RNG_SEED=123
dotest rngunif rngunif2.txt rngunif2.err ""
unset GSL_RNG_TYPE
unset GSL_RNG_SEED

export GSL_RNG_SEED=123
dotest randpoisson randpoisson2.txt randpoisson2.err ""
unset GSL_RNG_SEED

# write test.dat, perform test, and delete
$exdir/ntuplew
dotest ntupler ntuple.txt "" ""
rm -f test.dat

# test vector read/write
$exdir/vectorw
dotest vectorr vectorr.txt "" ""
rm -f test.dat

echo "FAIL: ${nfail}/${ntot}"
echo "PASS: ${npass}/${ntot}"

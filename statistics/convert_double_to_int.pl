#!/usr/local/bin/perl
$/ = undef ;
$_ = <> ;

s/a dataset/an integer dataset/g ;
s/datasets/integer datasets/g ;

s/const double data\[\]/const int data[]/g ;
s/const double data1\[\]/const int data1[]/g ;
s/const double data2\[\]/const int data2[]/g ;

s/double\s+(\w+)\s+=\s+data/int $1 = data/g ;

s/double\s+gsl_stats_max/int\ngsl_stats_max/mg ;
s/double\s+gsl_stats_min/int\ngsl_stats_min/mg ;

s/gsl_stats_/gsl_stats_i/g ;

s/sum_of_squares/isum_of_squares/g;
s/\#include \"/\#include \"i/g ;

print ;


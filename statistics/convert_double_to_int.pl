#!/usr/local/bin/perl
$/ = undef ;
$_ = <> ;

s/a dataset/an integer dataset/g ;
s/datasets/integer datasets/g ;

s/const double data\[\]/const int data[]/g ;
s/const double data1\[\]/const int data1[]/g ;
s/const double data2\[\]/const int data2[]/g;
s/const double sorted_data\[\]/const int sorted_data[]/g ;
s/double data\[\]/int data[]/g;
s/compare_doubles/compare_ints/g;
s/const double \*/const int */g;
s/sizeof\(double\)/sizeof(int)/g;

s/double\s+(\w+)\s+=\s+data/int $1 = data/g ;

s/double\s+gsl_stats_max/int\ngsl_stats_max/mg ;
s/double\s+gsl_stats_min/int\ngsl_stats_min/mg ;

s/gsl_stats_/gsl_stats_int_/g ;

s/sum_of_squares/int_sum_of_squares/g;
s/\#include \"(\w+)\.h\"/\#include \"int_$1.h\"/g ;

s/gsl_statistics.h/gsl_statistics_int.h/g ;
s/_GSL_STATISTICS_H/_GSL_STATISTICS_INT_H/g;

print ;


#!/usr/local/bin/perl
$/ = undef ;
$_ = <> ;

s/a dataset/an integer dataset/g ;
s/datasets/integer datasets/g ;

s/double data/int data/g ;
s/double sorted_data/int sorted_data/g ;
s/compare_doubles/compare_ints/g;
s/const double \*/const int \*/g;
s/sizeof\(double\)/sizeof(int)/g;

s/double\s+(\w+)\s+=\s+data/int $1 = data/g ;

s#double /\* BASE #int /\* BASE #g ;

s/gsl_stats_/gsl_stats_int_/g ;

s/sum_of_squares/int_sum_of_squares/g;
s/\#include \"(\w+)\.h\"/\#include \"int_$1.h\"/g ;

s/gsl_statistics.h/gsl_statistics_int.h/g ;
s/_GSL_STATISTICS_H/_GSL_STATISTICS_INT_H/g;

print "/* This is an automatically generated file created by\n" ;
print "   convert_double_to_int.pl -- do not edit                */\n\n";
print ;


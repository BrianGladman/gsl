#!/usr/bin/perl

open(FILE, "../gsl_version.h");
while (<FILE>) {
    next unless /define GSL_VERSION/;
    ($v) = /"(.*)"/;
    last ;
}
close (FILE);


while (<>) {
    s/\@VERSION\@/$v/g;
    print;
}

#!/usr/bin/perl

while (<>) {
    if (/^extern/ && !/^extern inline/) {
        $a = $_; $a =~ s/^extern/__declspec(dllexport)/;
        $b = $_; $b =~ s/^extern/__declspec(dllimport)/;
        print "#ifdef GSL_EXPORTS\n";
        print "$a";
        print "#elif defined(GSL_IMPORTS)\n";
        print "$b";
        print "#else\n";
        print "$_";
        print "#endif\n";
    } else {
        print;
    }
}

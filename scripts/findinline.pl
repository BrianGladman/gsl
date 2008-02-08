#!/usr/bin/perl
## find all the inline functions
## perl scripts/findinline.pl gsl/*.h | sort | uniq

while (<>) {
    if (/^(extern\s+)inline/) {
        if (s/\(.*//) {
            s/^(extern\s+)inline.*(gsl_|GSL_)/$2/;
            print ;
        } else {
            while (!s/\(.*//) { $_ = <> ; }
            print;
        }
    }
}

#!/usr/bin/perl
$/=';' ;

while (<>) {
    s/\/\*/\001/g;
    s/\*\//\002/g;
    s/\001[^\002]*//g;
    s/[^\001]*\002//g;
    s/^#.*//mg;
    s/__BEGIN_DECLS//g;
    s/__END_DECLS//g;
    s/\s+/ /g;
    s/^\s+//;
    s/;$//;
    s/\s+\)/\)/g;
    s/(\w)\(/$1 \(/g;
    next if !/\(.*\)\s*$/ ;
    s/\b_(gsl_\w+_view)/\1/g;
    s/(\w+)(\[|,|\))/\@var{$1}$2/g;
    s/\s*\,(\w)/, $1/g;
    s/\@var\{\*/* \@var\{/g ;
    s/\*\@var\{/* \@var\{/g ;
    s/\*(\S)/* $1/g ;
    s/\@var{void}/void/g;
    s/\(\s+/(/g;
    next if /^([a-z]+ )+gsl_(block|vector|matrix|stats|sort|permute|sort_vector|permute_vector)_(u?int|u?short|u?long|u?char|float|long_double|complex_float|complex_double|complex_long_double)/;
    print "\@deftypefun $_\n";
    print "\@end deftypefun\n\n" ;
}

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
    next if !/\(/ ;
    s/(\w+)(,|\))/\@var{$1}$2/g;
    s/\@var\{\*/*\@var\{/g ;
    print "\@deftypefun $_\n";
    print "\@end deftypefun\n\n" ;
}

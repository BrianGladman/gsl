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
    ($function) = /(\w+) \(/;
    next if $function !~ /^gsl_/;
    next if $function =~ /blas_raw/;
    $function =~ s/_impl$//;
    $function =~ s/_e$//;

    $base = $function;
    $base =~ s/_long_double_/_/;
    $base =~ s/_float_/_/;
    $base =~ s/_ulong_/_/;
    $base =~ s/_long_/_/;
    $base =~ s/_uint_/_/;
    $base =~ s/_int_/_/;
    $base =~ s/_ushort_/_/;
    $base =~ s/_short_/_/;
    $base =~ s/_uchar_/_/;
    $base =~ s/_char_/_/;
    $base =~ s/_long_double//;
    $base =~ s/_float//;
    $base =~ s/_ulong//;
    $base =~ s/_long//;
    $base =~ s/_uint//;
    $base =~ s/_int//;
    $base =~ s/_ushort//;
    $base =~ s/_short//;
    $base =~ s/_uchar//;
    $base =~ s/_char//;

    if ($base ne $function) {
        $secondary{$base} = $function;
    } else {        
        $primary{$function} = 1;
    }
}

for (sort keys %primary) {
    print $_,"\n" ;
}

for (sort keys %secondary) {
    print $secondary{$_},"\n" if !$primary{$_} ;
}

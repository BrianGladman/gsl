#!/usr/bin/perl

while (<>) {
    next unless /^\@def/ ;
    ($function) = /(\w+)\s+\(/;
    $function =~ s/_e$//;
    $function =~ s/_impl$//;

    $primary{$function} = 1;
}

for (sort keys %primary) {
    print $_,"\n" ;
}

#for (sort keys %secondary) {
#    print $secondary{$_},"\n" if !$primary{$_} ;
#}

while(<>){
    next unless /^\s*const\s+gsl_rng_type\s+\*\s*(\w+)/;
    print "$1\n" ;
}

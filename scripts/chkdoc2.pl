while(<>){
    next unless /^\@deffn\s+\S+\s+(\S+)/;
    print "$1\n" ;
}

#!/usr/bin/perl

# Use this script to add a constructor for the vector and matrix views objects.
# Some compilers will complain if there is no constructor for an object
# which contains a const element.
#
# usage:  perl ./scripts/fixcplusplus.pl -i.bak vector/gsl*.h matrix/gsl*.h

while (<>) {
    if (/\} (gsl_(vector|matrix)_.*_const_view)/) {
        print "#ifdef __cplusplus\n";
        print "   $1() { } ;\n";
        print "#endif\n" ;
    }
    print ;
}

#! /bin/sh

# Run this to generate all the auto-generated files needed by the GNU
# configure program

libtoolize
aclocal
automake --add-missing
autoconf

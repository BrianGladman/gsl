#! /bin/sh

echo 'This script will run aclocal, automake --add-missing and autoconf.'
echo 'You should have the automake and autoconf packages installed'
echo 'on your system.'

aclocal
automake --add-missing
autoconf

echo 'Now it will run configure with all the options you passed to ' $0
echo 'Running configure' $*

./configure $*

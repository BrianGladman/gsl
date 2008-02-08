#!/usr/bin/perl

$footer=q(<hr>The GNU Scientific Library - a free numerical library licensed under the GNU GPL<br>Back to the <a href="/software/gsl/">GNU Scientific Library Homepage</a>);

while (<>) {
    s/href=\"([^\"]+)\.html#\1\"/href=\"$1\.html\"/g;
    s/\s+(<\/body><\/html>)/$footer$1/;
    s/\s+(<\/pre>)/$1/;
    print;
}

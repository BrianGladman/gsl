#!/usr/local/bin/perl

print <<EOF ;
List of Known Problems
======================

The errors below have been automatically extracted from the output of
"make check" and are known to the developers.  

If you find a bug which is not on this list please report it to the
mailing list gsl-discuss\@sourceware.cygnus.com. Thank you.

EOF

print "-" x 79, "\n\n" ;

chomp(@lines = <>) ;

@n = () ;

for ($i = 0 ; $i < @lines ; $i++) {
    $_ = $lines[$i] ;
    if ((/^\S+: / || /^!/)
	&& !/PASS:/ 
	&& !/mdate-sh/
	&& !/cp: .\/libgsl.a.c: No such file or directory/
	&& !/sh internal 2K buffer overflow/ 
	&& !/cvs server: Updating \S+$/ 
	&& !/Entering directory/
	&& !/Leaving directory/
	&& !/Nothing to be done/) {
	$n[$i] = 2 ;
    }
}

for ($i = 0 ; $i < @lines ; $i++) {
    if ($n[$i]) {
	$c = 1 ;
    } else {
	$c = 0 if $lines[$i] =~ /^Making/ ;
	$n[$i] = 1 if $c > 0 ;
	$c-- if $c > 0 ;
    }
}

for ($i = @lines - 1 ; $i >= 0 ; $i--) {
    $c = 1 if $n[$i] ;
    next if $n[$i] ;
    $n[$i] = 1 if $c > 0 ;
    $c-- if $c > 0 ;
    $c = 0 if $lines[$i] =~ /^Making/ ;
}

for ($i = 0 ; $i < @lines ; $i++) {
    $_ = $lines[$i] ;
    if ($n[$i] > 1) {
	print "*** $_\n" ;
	$prev = 1 ;
    } elsif ($n[$i] == 1) {
	print "    $_\n" ;
	$prev = 1 ;
    } else {
	print "--------\n" if $prev == 1 ;
	$prev = 0 ;
    }
}


#! /usr/bin/perl

while (<DATA>) { print ; } ;
chomp(@lines = <>) ;

@n = () ;

for ($i = 0 ; $i < @lines ; $i++) {
    $_ = $lines[$i] ;
    if ((/^\S+: / || /^!/)
	&& !/PASS:/ 
	&& !/mdate-sh/
	&& !/cp: .\/libgsl.a.c: No such file or directory/
	&& !/config.status/
	&& !/sh internal 2K buffer overflow/ 
	&& !/cvs server: Updating \S+$/ 
	&& !/Entering directory/
	&& !/Leaving directory/
        && !/is up to date/
	&& !/Nothing to be done/) {
	$n[$i] = 2 ;
    }
}

for ($i = 0 ; $i < @lines ; $i++) {
    if ($n[$i]) {
	$c = 7 ;
    } else {
	$c = 0 if $lines[$i] =~ /^Making/ ;
	$n[$i] = 1 if $c > 0 ;
	$c-- if $c > 0 ;
    }
}

for ($i = @lines - 1 ; $i >= 0 ; $i--) {
    $c = 7 if $n[$i] ;
    next if $n[$i] ;
    $n[$i] = 1 if $c > 0 ;
    $c-- if $c > 0 ;
    $c = 0 if $lines[$i] =~ /^Making/ ;
}

for ($i = 0 ; $i < @lines ; $i++) {
    $_ = $lines[$i] ;
    $dir = $_ if /^Making/ ;
    if ($n[$i] > 1) {
	print "*** $_\n" ;
	$prev = 1 ;
    } elsif (/^Running/) {
        print "\n\n" ;
        print "=" x 75, "\n";
        print $_, "\n\n";
    } elsif ($n[$i] == 1) {
        print "$dir\n" if $dir ne $prevdir ;
        $prevdir = $dir ;
	print "    $_\n" ;
	$prev = 1 ;
    } else {
	print "--------\n" if $prev == 1 ;
	$prev = 0 ;
    }
}

######################################################################
__END__
List of Known Problems
======================

The errors below have been automatically extracted from the output of
"make check" and are known to the developers.  The output is in three
sections,

  1) make 
  2) make check  (using extended-precision floating point registers)
  3) make check  (using strict IEEE double-precision arithmetic)

If you find a bug which is not on this list please report it to the
mailing list bug-gsl@gnu.org. Thank you.

p.s. If you want to send us the output of 'make check' please extract
only the relevant parts, because the complete output is huge (over 10
megabytes).

-------------------------------------------------------------------------------

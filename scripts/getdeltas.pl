#!/usr/bin/perl -w
# -*- Mode: perl; indent-tabs-mode: nil -*-
#
# getdeltas.pl: Script to extract and calculate deltas between expected and actual GSL test results.
#
# Usage: getdeltas.pl must be run from the root build directory.
#        One or more modules can be specified on the commandline; default is all.
#        Alternatively, modules can be excluded by preceding the list with -x.
#        Test results are expected to be in $testout (see filename variables below) in each module.
#        Output is written to $gsldeltas with a timestamp extension.
#        The PASS and FAIL counts are slightly skewed by check-TESTS and built-in summary lines.
#
# copyright (c) 2001 Henry Sobotka <sobotka@axess.com>
# include GNU General Public License
#
use strict;

# Filenames
#
my $testout   = "test.out";
my $gsldeltas = "gsldeltas";
my $versionh  = "gsl_version.h";

# GSL subdirectories with tests
#
my @mods = qw(sys err complex block vector matrix permutation sort ieee-utils blas linalg eigen
	      specfunc dht qrng rng randist fft poly fit multifit statistics siman sum integration
	      interpolation histogram ode-initval roots multiroots min multimin monte ntuple diff);

# Globals
#
my $header = '';
my %deltas = ();
my %dmods  = ();
my @mlist  = ();
my $passes = 0;
my $fails  = 0;

# Get version and format output header
#
sub setheader {

    my $vstr = "#define GSL_VERSION ";

    open(VERSIONH, "<$versionh") or die "Can't open $versionh: $!\n";
    my @version = grep /$vstr/, <VERSIONH>;
    close VERSIONH;

    my $ver = substr($version[0], length($vstr));
    $ver =~ s/"//g;
    chomp $ver;

    my ($sec, $min, $hour, $mday, $mon, $year) = (localtime)[0..5];
    $year += 1900;
    $mon += 1;
    my $timestamp = sprintf "%d/%02d/%02d %02d:%02d:%02d", $year, $mon, $mday, $hour, $min, $sec;
    $header = "GNU Scientific Library $ver test result deltas - $^O - $timestamp";
}

# Find deltas in $testout files
#
sub parsetestout {

    my @dirs = ();
    my $dir = shift;

    if ($dir eq 'all') {
	push @mlist, "ALL"; @dirs = @mods;
    }
    else {
	do { push @dirs, $dir; } while (defined($dir = shift));

	if ($ARGV[0] eq "-x") {
	    @mlist = grep !/-x/, @ARGV;
	    unshift @mlist, "ALL EXCEPT";
	}
	else { @mlist = @dirs; }
    }

    foreach my $d (@dirs) {

	my $c = 0;
	my $file = $d.'/'.$testout;
	if (!-e $file) { print "ERROR: Can't find $file\n"; next; }
	open(TESTOUT, "<$file") or die "Can't open $file: $!\n";

	while (defined (my $line = <TESTOUT>)) {

	    ++$passes if ($line =~ /^PASS:/);
	    ++$fails if ($line =~ /^FAIL:/);
	    next unless $line =~ /obs/;

	    my $vals = substr ($line, (rindex($line, "(") + 1));
	    my ($obs, $exp);

	    if ($line =~ /theoretical/) {
		$vals =~ /(.+) vs (.+)\)/;
		$obs = $1; $exp = $2;
	    }
	    elsif ($line =~ /observed/) {
		$vals =~ /(.+) observed vs (.+) exp/;
		$obs = $1; $exp = $2;
	    }
	    elsif ($line =~ /obs vs/) {
		$vals =~ /(.+) obs vs (.+) exp/;
		$obs = $1; $exp = $2;
	    }
	    elsif ($line =~ /vs plain/) {
		$vals =~ /obs (.+) vs plain (.+)\)/;
		$obs = $1; $exp = $2;
	    }
	    elsif ($line =~ /vs exp/) {
		$vals =~ /obs (.+) vs exp (.+)\)/;
		$obs = $1; $exp = $2;
	    }

	    if (defined($obs) && defined($exp)) {
		if ($obs != $exp) {
		    my $diff = $obs - $exp;
		    push @{$deltas{$diff}}, "[$d] $line";
		    ++$c;
		}
	    }
	}
	close TESTOUT;
	$dmods{$d} = $c unless ($c == 0);
    }
}

# Write deltas to file
#
sub writedeltas {

    my $ext = time;
    open (DELTAS, ">$gsldeltas-$ext") or die "Can't open $gsldeltas-$ext: $!\n";
    print DELTAS "$header\n";
    print DELTAS "=" x length($header);
    print DELTAS "\nModules: @mlist\n\n";

    my $c = 0;
    foreach my $d (sort {$b <=> $a} keys (%deltas)) {
	print DELTAS "DELTA $d\n";
	foreach my $t (@{$deltas{$d}}) {
	    print DELTAS "   $t";
	    ++$c;
	}
	print DELTAS "\n";
    }
    print DELTAS "\n", '-' x 24, " SUMMARY ", '-' x 24;
    printf DELTAS "\n\n%8d deltas TOTAL (counted %d PASS, %d FAIL)\n", $c, $passes, $fails;
    foreach my $k (sort keys (%dmods)) {
	printf DELTAS "%8d in %s\n", $dmods{$k}, $k;
    }
    close DELTAS;
    print "Wrote $c entries to $gsldeltas-$ext\n";
}

sub getmodules {

    if (defined ($ARGV[0])) {

	if ($ARGV[0] eq "-x") {
	    my @dirs = ();
	    foreach my $m (@mods) {
		push @dirs, $m unless (grep /$m/, @ARGV);
	    }
	    parsetestout @dirs;
	}
	else { parsetestout @ARGV; }
    }
    else { parsetestout 'all'; }
}


# main
#
&setheader;
&getmodules;
&writedeltas;
exit;

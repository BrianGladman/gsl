#!/usr/bin/perl

push(@ARGV,"Makefile") ;

while (<>) {
    $line = "" ;

    if (s/\s+\\$/ /) {
	chomp;
	$line = $_ ;
	while (<>) {
	    chomp ;
	    if (s/\s+\\$/ /) {
		$line .= $_ ;
	    } else {
		$line .= $_ ;
		last ;
	    }
	}
    } else {
	$line = $_ ;
    }

    if ($line =~ /SUBDIRS\s*=\s*(.*)/) {
	@dirs = split(' ',$1) ;
	for (@dirs) {
	    push (@ARGV, "$_/Makefile") ;
	}
    }

    if ($line =~ /(libgsl\w+)_a_OBJECTS\s*=\s*(.*)/) {
	$lib = $1 ;
	($dir) = split('/', $ARGV) ;
	@files = split(' ',$2) ;
	@files = map("$dir/$_", @files) ;
	$cmd = "\t" . "\$(AR) q libgsl.a " . join(" ",@files) ;
	push(@cmds,$cmd) ;
	push(@deps,"$dir/$lib.a") ;
    }
}

open(TMP,">Makefile.am.tmp") ;
open(FILE,"<Makefile.am") ;
select(TMP) ;
while (<FILE>) {
    last if /^libgsl.a:/ ;
    print ;
}
close (FILE) ;

print "libgsl.a: ", join(" ",@deps), "\n" ;
print "\t", "\@rm -f libgsl.a\n" ;
print join("\n", @cmds) ;
print "\n" ;
print "\t", "\$(RANLIB) libgsl.a\n" ;
close(TMP) ;

rename("Makefile.am.tmp","Makefile.am");

#!/usr/bin/perl
while (<>) {
    if (/^\@c\{\$(.*)\$\}\s*$/) {
	$expected = "\@math{$1}" ;
	$n = length($expected) ;
	$next_line = <> ;
	print $_ if (substr($next_line, 0, $n) ne $expected) ;
	print $next_line ;
	next ;
    }

    $in_example = 1 if /^\@example/ ;
    $in_example = 0 if /^\@end example/ ;

    s/ \\over /\//g if $in_example ;
	
    print $_ ;
}

#!/usr/bin/perl

if (@ARGV < 2) {
    die "addgpl.pl authors [files]\n";
};

$authors = shift(@ARGV) ;
@DATA = <DATA>;

for $file (@ARGV) {
    next unless $file =~ /\.c$|\.h$/ ;

    open(FILE,"<$file") ;
    @lines = <FILE> ;
    close(FILE);

    next if grep(/copyright/i, @lines);  # skip any existing copyrights
    next if !grep(/^\#|^\s*[a-z]/, @lines) ;  # doesn't contain any code
    
    print "$file $authors\n" ;
    next if $authors =~ /--test/i;

    $header = &notice($authors, $file) ;
    rename("$file", "$file.bak") if ! -e "$file.bak" ;
    open(OUT,">$file");
    print OUT $header, "\n" ;
    print OUT "\n" ;
    print OUT @lines ;
    close(OUT);
}

sub notice {
    my ($authors,$file) = @_ ;
    
    use Text::Wrap;
    $Text::Wrap::columns=78;

    my @header = @DATA;

    for (@header) {
        if ( s/\@AUTHORS\@/$authors/ ) {
            $_ = wrap ('', '', $_);
        } else {
            s/\@FILE\@/$file/;
        }
    }

    my $header = join('',@header);
    $header =~ s/^/\/* /;
    $header =~ s/\n/\n * /mg;
    $header =~ s/\* $/\*\//;

    return $header ;
}
   

__END__
@FILE@

Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 @AUTHORS@

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

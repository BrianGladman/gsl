while (<>) {
    s/\s*$//;
    if (/^\@example/ .. /^\@end example/) {
        push (@lines, $_);
    } else {
        #print if length($_) > 55;
        if (@lines) {
            my $long = grep(length($_) > 55, @lines);
            if ($long) {
                map(s/^\@example/\@smallexample/, @lines);
                map(s/^\@end example/\@end smallexample/, @lines);
                for ($i = 0 ; $i < @lines; $i++) { $lines[$i] .= " "x55;  substr($lines[$i],55,1) = "|" ; $lines[$i] =~ s/\s*$//;} ;
                print join("\n", @lines), "\n" ;
            };
            undef @lines ;
        } ;
        #print;
    }

}

for $file (@ARGV) {
    open (FILE, "<$file") ;
    $c = 0;
    $l = 0;
    while (<FILE>) {
        $c++ if /\@.index/ ;
        $l++ ;
    }
    $score{$file} = $c/$l ;
}

@files =  sort {$score{$b} <=> $score{$a}} keys %score ;

for (@files) {
    printf ("%30s %f\n", $_, $score{$_}) ;
}

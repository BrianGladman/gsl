@groupa = ( .0421, .0941, .1064, .0242, .1331, .0773, .0243, .0815,
	   .1186, .0356, .0728, .0999, .0614, .0479 ) ;

@groupb = ( .1081, .0986, .1566, .1961, .1125, .1942, .1079, .1021,
	   .1583, .1673, .1675, .1856, .1688, .1512 ) ;

@igroupa = ( 17 , 18 , 16 , 18 , 12 , 20 , 18 , 20 , 20 , 22 , 20 , 10
	    , 8 , 12 , 16 , 16 , 18 , 20 , 18 , 21 ) ;

@igroupb = ( 19 , 20 , 22 , 24 , 10 , 25 , 20 , 22 , 21 , 23 , 20 , 10
	    , 12 , 14 , 12 , 20 , 22 , 24 , 23 , 17 ) ;

print "mean groupa = ", &mean(@groupa), "\n" ;
print "mean groupb = ", &mean(@groupb), "\n" ;
print "mean igroupa = ", &mean(@igroupa), "\n" ;
print "mean igroupb = ", &mean(@igroupb), "\n" ;

print "variance groupa = ", &variance(@groupa), "\n" ;
print "variance groupb = ", &variance(@groupb), "\n" ;
print "variance igroupa = ", &variance(@igroupa), "\n" ;
print "variance igroupb = ", &variance(@igroupb), "\n" ;

print "variance_est groupa = ", &variance_est(@groupa), "\n" ;
print "variance_est groupb = ", &variance_est(@groupb), "\n" ;
print "variance_est igroupa = ", &variance_est(@igroupa), "\n" ;
print "variance_est igroupb = ", &variance_est(@igroupb), "\n" ;

print "sd groupa = ", &sd(@groupa), "\n" ;
print "sd groupb = ", &sd(@groupb), "\n" ;
print "sd igroupa = ", &sd(@igroupa), "\n" ;
print "sd igroupb = ", &sd(@igroupb), "\n" ;

print "sd_est groupa = ", &sd_est(@groupa), "\n" ;
print "sd_est groupb = ", &sd_est(@groupb), "\n" ;
print "sd_est igroupa = ", &sd_est(@igroupa), "\n" ;
print "sd_est igroupb = ", &sd_est(@igroupb), "\n" ;

print "pv groupa groupb = ", &pv(\@groupa,\@groupb), "\n" ;
print "pv igroupa igroupb = ", &pv(\@igroupa,\@igroupb), "\n" ;

print "t groupa groupb = ", &t(\@groupa,\@groupb), "\n" ;
print "t igroupa igroupb = ", &t(\@igroupa,\@igroupb), "\n" ;



sub mean {
    my @x = @_ ;
    my $sum ;
    for $x (@x) { $sum += $x ; } ;
    return $sum / scalar(@x) ;
}

sub variance {
    my @x = @_ ;
    my $sum ;
    my $mean = &mean(@x) ;
    for $x (@x) { $sum += ($x-$mean)*($x-$mean) ; } ;
    return $sum / scalar(@x) ;
}

sub variance_est {
    my @x = @_ ;
    my $sum ;
    my $mean = &mean(@x) ;
    for $x (@x) { $sum += ($x-$mean)*($x-$mean) ; } ;
    return $sum / (scalar(@x)-1) ;
}

sub sd {
    my @x = @_ ;
    return sqrt(&variance(@x)) ;
}   

sub sd_est {
    my @x = @_ ;
    return sqrt(&variance_est(@x)) ;
}   

sub pv {
    my ($x, $y) = @_ ;
    my $v1 = &variance_est(@$x) ;
    my $v2 = &variance_est(@$y) ;
    my $n1 = scalar(@$x) ;
    my $n2 = scalar(@$y) ;
    return ((($n1 - 1)*$v1)+(($n2-1)*$v2)) / ($n1 + $n2 - 2);
}
	
sub t {
    my ($x, $y) = @_ ;
    my $mean1 = &mean(@$x);
    my $mean2 = &mean(@$y);
    my $sd1 = &sd_est(@$x);
    my $sd2 = &sd_est(@$y);
    my $pv = &pv($x,$y); 
    my $n1 = scalar(@$x) ;
    my $n2 = scalar(@$y) ;
    
    my $t = ($mean1-$mean2)/(sqrt($pv*((1.0/$n1)+(1.0/$n2))));
  
    return $t;
}

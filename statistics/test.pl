@weight = (.0000, .0000, .0000, 3.000, .0000, 1.000, 1.000, 1.000,
           0.000, .5000, 7.000, 5.000, 4.000, 0.123) ;

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

print "absdev groupa = ", &absdev(@groupa), "\n" ;
print "absdev groupb = ", &absdev(@groupb), "\n" ;
print "absdev igroupa = ", &absdev(@igroupa), "\n" ;
print "absdev igroupb = ", &absdev(@igroupb), "\n" ;

print "skew groupa = ", &skew(@groupa), "\n" ;
print "skew groupb = ", &skew(@groupb), "\n" ;
print "skew igroupa = ", &skew(@igroupa), "\n" ;
print "skew igroupb = ", &skew(@igroupb), "\n" ;

print "kurt groupa = ", &kurt(@groupa), "\n" ;
print "kurt groupb = ", &kurt(@groupb), "\n" ;
print "kurt igroupa = ", &kurt(@igroupa), "\n" ;
print "kurt igroupb = ", &kurt(@igroupb), "\n" ;

print "median groupa = ", &median(@groupa), "\n" ;
print "median groupb = ", &median(@groupb), "\n" ;
print "median igroupa = ", &median(@igroupa), "\n" ;
print "median igroupb = ", &median(@igroupb), "\n" ;

print "median_less_one groupa = ", &median_less_one(@groupa), "\n" ;
print "median_less_one groupb = ", &median_less_one(@groupb), "\n" ;
print "median_less_one igroupa = ", &median_less_one(@igroupa), "\n" ;
print "median_less_one igroupb = ", &median_less_one(@igroupb), "\n" ;

print "cov groupa groupb = ", &cov(\@groupa,\@groupb), "\n" ;
print "cov igroupa igroupb = ", &cov(\@igroupa,\@igroupb), "\n" ;

print "pv groupa groupb = ", &pv(\@groupa,\@groupb), "\n" ;
print "pv igroupa igroupb = ", &pv(\@igroupa,\@igroupb), "\n" ;

print "t groupa groupb = ", &t(\@groupa,\@groupb), "\n" ;
print "t igroupa igroupb = ", &t(\@igroupa,\@igroupb), "\n" ;

print "weighted mean groupa = ", &weighted_mean(\@weight, \@groupa), "\n" ;
print "weighted variance groupa = ", &weighted_variance(\@weight, \@groupa), "\n" ;
print "weighted variance_est groupa = ", &weighted_variance_est(\@weight, \@groupa), "\n" ;
print "weighted sd groupa = ", &weighted_sd(\@weight, \@groupa), "\n" ;
print "weighted sd_est groupa = ", &weighted_sd_est(\@weight, \@groupa), "\n" ;
print "weighted absdev groupa = ", &weighted_absdev(\@weight, \@groupa), "\n" ;
print "weighted skew groupa = ", &weighted_skew(\@weight, \@groupa), "\n" ;
print "weighted kurt groupa = ", &weighted_kurt(\@weight, \@groupa), "\n" ;

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

sub absdev {
    my @x = @_ ;
    my $sum ;
    my $mean = &mean(@x) ;
    for $x (@x) { $sum += abs($x-$mean) ; } ;
    return $sum / scalar(@x) ;
}

sub skew {
    my @x = @_ ;
    my $sum ;
    my $mean = &mean(@x) ;
    my $sd = &sd_est(@x) ;
    for $x (@x) { $d = ($x - $mean)/$sd ; $sum += $d*$d*$d ; } ;
    return $sum / scalar(@x) ;
}


sub kurt {
    my @x = @_ ;
    my $sum ;
    my $mean = &mean(@x) ;
    my $sd = &sd_est(@x) ;
    for $x (@x) { $d = ($x - $mean)/$sd ; $sum += $d*$d*$d*$d ; } ;
    return $sum / scalar(@x) - 3 ;
}

sub weighted_mean {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    my $sum; my  $den ;
    for (my $i = 0 ; $i < @x ; $i++) { 
        $sum += $w[$i] * $x[$i] ; 
        $den += $w[$i] ;
    } ;
    return $sum / $den ;
}

sub weighted_variance {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    my $sum ; my $den;
    my $mean = &weighted_mean(\@w, \@x) ;
    for (my $i = 0 ; $i < @x ; $i++) { 
        $sum += $w[$i] * ($x[$i]-$mean)*($x[$i]-$mean) ; 
        $den += $w[$i] ;
    } ;
    return $sum / $den ;
}

sub weighted_variance_est {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    my $sum ; my $sum_w; my $sum_w2 ;
    my $mean = &weighted_mean(\@w,\@x) ;
    for (my $i = 0 ; $i < @x ; $i++) { 
        $sum += $w[$i] * ($x[$i]-$mean)*($x[$i]-$mean) ; 
        $sum_w += $w[$i] ;
        $sum_w2 += $w[$i] * $w[$i] ;
    } ;
    
    return ($sum_w / ($sum_w * $sum_w - $sum_w2)) * $sum ;
}

sub weighted_sd {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    return sqrt(&weighted_variance(\@w,\@x)) ;
}   

sub weighted_sd_est {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    return sqrt(&weighted_variance_est(\@w,\@x)) ;
}   

sub weighted_absdev {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    my $sum ; my $den;
    my $mean = &weighted_mean(\@w, \@x) ;
    for (my $i = 0 ; $i < @x ; $i++) { 
        $sum += $w[$i] * abs($x[$i]-$mean) ; 
        $den += $w[$i] ;
    } ;
    return $sum / $den;
}

sub weighted_skew {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    my $sum ; my $den ;
    my $mean = &weighted_mean(\@w, \@x) ;
    my $sd = &weighted_sd_est(\@w, \@x) ;
    for (my $i = 0 ; $i < @x ; $i++) { 
        $d = ($x[$i] - $mean)/$sd ; 
        $sum += $w[$i] * $d*$d*$d ; 
        $den += $w[$i] ;
    } ;
    return $sum / $den ;
}

sub weighted_kurt {
    my ($w_ref, $x_ref) = @_ ; my @w = @$w_ref ; my @x = @$x_ref ;
    my $sum ; my $den ;
    my $mean = &weighted_mean(\@w,\@x) ;
    my $sd = &weighted_sd_est(\@w,\@x) ;
    for (my $i = 0 ; $i < @x ; $i++) { 
        $d = ($x[$i] - $mean)/$sd ; 
        $sum += $w[$i]*$d*$d*$d*$d ; 
        $den += $w[$i] ;} ;
    return( $sum / $den) - 3 ;
}


sub median {
    my @x = @_ ;
    @x = sort {$a <=> $b} @x ;
    print "x= ",join(",",@x),"\n" ;
    my $n = scalar(@x);
    my $median  ;
    if ($n % 2 == 1) {
	$median = $x[$n/2] ;
    } else {
	$median = 0.5*($x[$n/2 - 1] + $x[$n/2]) ;
    }
    return $median ;
}

sub median_less_one {
    my @x = @_ ;
    @x = sort @x ;
    pop(@x) ;
    return &median(@x) ;
}

sub cov {
    my ($x, $y) = @_ ;
    my $mean1 = &mean(@$x);
    my $mean2 = &mean(@$y);
    my $n1 = scalar(@$x) ;
    my $n2 = scalar(@$y) ;
    my $cov = 0;
    for (my $i = 0 ; $i < $n1 ; $i++) {
        print "i = $i d1 = ", $x->[$i] - $mean1, " d2 = ", $y->[$i] - $mean2, "\n" ;
        $cov += ($x->[$i] - $mean1) * ($y->[$i] - $mean2);
    }
    return $cov / ($n1 - 1); 
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


sub do_math {
    local($_) = $_[0];
    return $_ if !$htmlhelp;
    print STDERR "math on $_\n";

    &unprotect_texi;
    &unprotect_html;

    s/\@(\{|\})/$1/g;
    s/(\\le|<=)\b/&le;/g;
    s/(\\ge|>=)\b/&ge;/g;
    s/\\ne\b/&ne;/g;
    s/(\\lt|<)\b/&lt;/g;
    s/(\\gt|>)\b/&gt;/g;

    s/\_(\w)\b/<sub>$1<\/sub>/g;
    s/\^(\w)\b/<sup>$1<\/sup>/g;

    s/\_(\\?[a-z]+)\b/<sub>$1<\/sub>/ig;
    s/\^(\\?[a-z]+)\b/<sup>$1<\/sup>/ig;

    s/(\d)\^(-?\d+)/$1<sup>$2<\/sup>/g;

    s/\_\{([^\{\}]+)\}/<sub>$1<\/sub>/g;
    s/\^\{([^\{\}]+)\}/<sup>$1<\/sup>/g;
    s/\_(\([^\(\)]+\))/<sub>$1<\/sub>/g;
    s/\^(\([^\(\)]+\))/<sup>$1<\/sup>/g;

    s/\\infty\b/&infin;/g;
    s/\\approx\b/&cong;/g;
    s/\\times\b/&times;/g;
    s/\\int/&int;/g;
    s/\\(Re|Im)/$1/g;

    s/\\(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma|tau|upsilon|phi|chi|psi|omega)\b/&$1;/g;
    s/\\(Alpha|Beta|Gamma|Delta|Epsilon|Zeta|Eta|Theta|Iota|Kappa|Lambda|Mu|Nu|Xi|Omicron|Pi|Rho|Sigma|Tau|Upsilon|Phi|Chi|Psi|Omega)\b/&$1;/g;
    s/\\(ln|exp|log|sin|cos|tan|sinh|cosh|tanh|sec|cosec|cot|sech|cosech|csch|coth)\b/$1/g;
    s/\\(arcsin|arccos|arctan|arcsinh|arccosh|arctanh|arcsec|arccsc|arccot|arcsech|arccsch|arccoth)\b/$1/g;


    print STDERR "gives $_\n";
    return $_;
}

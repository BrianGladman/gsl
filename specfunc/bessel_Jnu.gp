besselJnu(nu,x) =
{
        (x/2)^nu * suminf(k=0, (-(x^2)/4)^k/(k! * gamma(k+nu+1)))
}
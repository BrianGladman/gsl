## Copyright (C) 2010 Pedro Gonnet <pedro.gonnet@gmail.com>
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{int}, @var{err}, @var{nr_points}] =} cquad (@var{f}, @var{a}, @var{b}, @var{tol})
## @deftypefnx {Function File} {[@var{int}, @var{err}, @var{nr_points}] =} cquad (@var{f}, @var{a}, @var{b}, @var{tol}, @var{sing})
## Numerically evaluates an integral using the doubly-adaptive
## quadrature described by P. Gonnet in @cite{"Increasing the
## Reliability of Adaptive Quadrature Using Explicit Interpolants",
## ACM Transactions on Mathematical Software, in Press, 2010}.
## The algorithm uses Clenshaw-Curtis quadrature rules of increasing
## degree in each interval and bisects the interval if either the
## function does not appear to be smooth or a rule of maximum
## degree has been reached. The error estimate is computed from the
## L2-norm of the difference between two successive interpolations
## of the integrand over the nodes of the respective quadrature rules.
## 
## For example,
##
## @example
##    int = cquad ( f , a , b , 1.0e-6 );
## @end example
##
## @noindent computes the integral of a function @var{f} in the interval
## [@var{a},@var{b}] to the relative precision of six
## decimal digits.
## The integrand @var{f} should accept a vector argument and return a vector
## result containing the integrand evaluated at each element of the
## argument, for example
##
## @example
##    f = @@(x) x .* sin ( 1 ./ x ) .* sqrt ( abs ( 1 - x ) );
## @end example
##
## If the integrand has known singularieites or discontinuities
## in any of its derivatives inside the interval,
## as does the above example at x=1, these can be specified in
## the additional argument @var{sing} as follows
##
## @example
##    int = cquad ( f , a , b , 1.0e-6 , [ 1 ] );
## @end example
##
## The two additional output variables @var{err} and @var{nr_points}
## return an estimate of the absolute integration error and
## the number of points at which the integrand was evaluated
## respectively.
## If the adaptive integration did not converge, the value of
## @var{err} will be larger than the requested tolerance. It is
## therefore recommended to verify this value for difficult
## integrands.
##
## If either @var{a} or @var{b} are @code{+/-Inf}, @code{cquad}
## integrates @var{f} by substituting the variable of integration
## with @code{x=tan(pi/2*u)}.
##
## @code{cquad} is capable of dealing with non-numerical
## values of the integrand such as @code{NaN}, @code{Inf}
## or @code{-Inf}, as the above example at x=0.
## If the integral diverges and @code{cquad} detects this, 
## a warning is issued and @code{Inf} or @code{-Inf} is returned.
##
## Note that @code{cquad} is a general purpose quadrature algorithm
## and as such may be less efficient for smooth or otherwise
## well-behaved integrand than other methods such as
## @code{quadgk} or @code{trapz}.
##
## @seealso{triplequad, dblquad, quadgk, quadl, quadv, trapz}
## @end deftypefn

## Author: Pedro G. Gonnet <pedro.gonnet@gmail.com>
## Keywords: Quadrature, Integration

function [ int , err , nr_points ] = cquad ( f , a , b , tol , sing )

    % declare persistent variables
    persistent n xi bee ...
        V Vinv Vcond T_left T_right w U
    
    % have the persistent variables been declared already?
    if ( ~exist ('U') || isempty (U) )
    
        % the nodes of the different rules and the coefficients
        % of their newton polynomials
        n = [4,8,16,32];
        xi{1} = sin(pi*(-n(1):2:n(1))/(2*n(1)))';
        xi{2} = sin(pi*(-n(2):2:n(2))/(2*n(2)))';
        xi{3} = sin(pi*(-n(3):2:n(3))/(2*n(3)))';
        xi{4} = sin(pi*(-n(4):2:n(4))/(2*n(4)))';
        bee{1} = [0., .233284737407921723637836578544e-1, 0., -.831479419283098085685277496071e-1, 0.,0.0541462136776153483932540272848 ]';
        bee{2} = [0., .883654308363339862264532494396e-4, 0., .238811521522368331303214066075e-3, 0., .135365534194038370983135068211e-2, 0., -.520710690438660595086839959882e-2, 0.,0.00341659266223572272892690737979 ]';
        bee{3} = [0., .379785635776247894184454273159e-7, 0., .655473977795402040043020497901e-7, 0., .103479954638984842226816620692e-6, 0., .173700624961660596894381303819e-6, 0., .337719613424065357737699682062e-6, 0., .877423283550614343733565759649e-6, 0., .515657204371051131603503471028e-5, 0.,-.203244736027387801432055290742e-4, 0.,0.0000134265158311651777460545854542 ]';
        bee{4} = [0., .703046511513775683031092069125e-13, 0., .110617117381148770138566741591e-12, 0., .146334657087392356024202217074e-12, 0., .184948444492681259461791759092e-12, 0., .231429962470609662207589449428e-12, 0., .291520062115989014852816512412e-12, 0., .373653379768759953435020783965e-12, 0., .491840460397998449461993671859e-12, 0., .671514395653454630785723660045e-12, 0., .963162916392726862525650710866e-12, 0., .147853378943890691325323722031e-11, 0., .250420145651013003355273649380e-11, 0., .495516257435784759806147914867e-11, 0., .130927034711760871648068641267e-10, 0., .779528640561654357271364483150e-10, 0., -.309866395328070487426324760520e-9, 0., .205572320292667201732878773151e-9]';

        % compute the Vandermonde-like Legendre matrices
        xil = (xi{4} - 1) / 2; xir = (xi{4} + 1) / 2;
        Vx = ones ( n(4)+1 ); Vl = ones ( n(4)+1 ); Vr = ones ( n(4)+1 );
        Vx(:,2) = xi{4}; Vl(:,2) = xil; Vr(:,2) = xir;
        for i=3:n(4)+1
            Vx(:,i) = (2*i-3) / (i-1) * xi{4} .* Vx(:,i-1) - (i-2) / (i-1) * Vx(:,i-2);
            Vl(:,i) = (2*i-3) / (i-1) * xil .* Vl(:,i-1) - (i-2) / (i-1) * Vl(:,i-2);
            Vr(:,i) = (2*i-3) / (i-1) * xir .* Vr(:,i-1) - (i-2) / (i-1) * Vr(:,i-2);
        endfor;
        for i=1:n(4)+1
            Vx(:,i) = Vx(:,i) * sqrt(4*i-2)/2;
            Vl(:,i) = Vl(:,i) * sqrt(4*i-2)/2;
            Vr(:,i) = Vr(:,i) * sqrt(4*i-2)/2;
        endfor;
        V{4} = Vx; Vinv{4} = inv ( V{4} );
        V{3} = V{4}(1:2:end,1:n(3)+1); Vinv{3} = inv ( V{3} );
        V{2} = V{4}(1:4:end,1:n(2)+1); Vinv{2} = inv ( V{2} );
        V{1} = V{4}(1:8:end,1:n(1)+1); Vinv{1} = inv ( V{1} );
        Vcond = [ cond(V{1}) , cond(V{2}) , cond(V{3}) , cond(V{4}) ];

        % compute the shift matrices
        T_left = Vinv{4} * Vl; 
        T_right = Vinv{4} * Vr;

        % compute the integration vector
        w = [ sqrt(2) , zeros(1,n(4)) ] / 2; % legendre
        
        % set-up the downdate matrix
        k = [0:n(4)]';
        U = diag(sqrt((k+1).^2 ./ (2*k+1) ./ (2*k+3))) + diag(sqrt(k(3:end).^2 ./ (4*k(3:end).^2-1)),2);
        
    endif; % if exist('n')
        
    % create the original datatype
    ivals = struct( ...
        'a', [], 'b',[], ...
        'c', [], ...
        'fx', [], ...
        'int', [], ...
        'err', [], ...
        'depth', [], ...
        'rdepth', [], ...
        'ndiv', [] );
    
    % onvert function given as a string to a function handle (from quadgk)
    if ( ischar (f) )
      f = @(x) feval (f, x);
    endif;
    
    % get the interval(s)
    if ( ~exist ('sing') || isempty (sing) ), iv = [ a ; b ];
    else iv = [ a ; sort(sing(:)) ; b ]; endif;
    
    % if a or b are +/- Inf, transform the interval
    % NOTE: there are probably better ways of doing this, e.g. what
    % quadgk or Waldvogel does!
    if ( isinf (a) || isinf (b) )
        f = @(x) f ( tan ( pi / 2 * x ) ) .* ( 1 + tan ( pi * x / 2 ).^2 ) * pi / 2;
        if ( isinf (a) ), a = -1;
        else a = 2 * atan (a) / pi; endif;
        if ( isinf (b) ), b = 1;
        else b = 2 * atan (b) / pi; endif;
        iv = [ a ; 2 * atan( iv(2:end-1) ) / pi ; b ];
    endif;

    % compute the first interval(s)
    for k=1:length(iv)-1
    
        % evaluate f in the interval
        m = (iv(k) + iv(k+1)) / 2; h = (iv(k+1) - iv(k)) / 2;
        points = m + h * xi{4};
        fx = f (points);
        
        % check for nans and clean them up
        nans = find ( ~isfinite ( fx ) );
        fx(nans) = 0;
        
        % compute the coefficients of the two highest degree interpolations
        ivals(k).c = zeros(n(4)+1,4);
        ivals(k).c(1:n(4)+1,4) = Vinv{4} * fx;
        ivals(k).c(1:n(3)+1,3) = Vinv{3} * fx([1:2:n(4)+1]);
        
        % re-instate the nans
        fx(nans) = NaN;
        
        % store the data
        ivals(k).fx = fx;
        ivals(k).a = iv(k); ivals(k).b = iv(k+1);
        ivals(k).depth = 4;
        ivals(k).ndiv = 0;
        ivals(k).rdepth = 1;
        
        % compute the integral
        ivals(k).int = 2 * h * w * ivals(k).c(:,4);
        
        % compute the error estimate
        c_diff = norm(ivals(k).c(:,4) - ivals(k).c(:,3));
        ivals(k).err = 2 * h * c_diff;
        if ( c_diff / norm(ivals(k).c(:,4)) > 0.1 ),  ivals(k).err = max( ivals(k).err , 2 * h * norm(ivals(k).c(:,4)) ); endif;
        
    endfor;
    
    % init some globals
    int = sum( [ ivals(:).int ] );
    err = sum( [ ivals(:).err ] );
    nr_ivals = length(ivals);
    int_final = 0; err_final = 0; err_excess = 0;
    nr_points = nr_ivals * (n(4) + 1);
    ndiv_max = 20;
    [ dummy , i_max ] = max( [ ivals(:).err ] );
    
    % do we even need to go this way?
    if ( err < int * tol ), return; endif;
    
    % main loop
    while ( true )
    
        % get some global data
        iv = ivals(i_max);
        a = iv.a; b = iv.b;
        depth = iv.depth; 
        split = 0;
        
        % check the depth of the interval. if it is less than 4,
        % then try to increase the degree. if this fails or the
        % interval is at depth 4, split the interval.
        
        % interval of depth 1
        if ( depth < 4 )
        
            % adjust the depth
            iv.depth = depth+1;
            depth = depth + 1;
            
            % get the new function values
            points = (a+b)/2 + (b-a)*xi{depth}/2;
            fx(1:2:n(depth)+1) = iv.fx;
            fx(2:2:n(depth)) = f (points(2:2:n(depth)));
            fx = fx(1:n(depth)+1);
            nr_points = nr_points + n(depth)-n(depth-1);
            
            % purge nans and compute the coefficients
            nans = find( ~isfinite ( fx ) );
            fx(nans) = 0;
            c_new = Vinv{depth} * fx;
            
            % downdate to remove any nans
            if ( length(nans) > 0 )
                b_new = bee{depth}; n_new = n(depth);
                for i=nans'
                    b_new(1:end-1) = (U(1:n(depth)+1,1:n(depth)+1) - diag(ones(n(depth),1)*xi{depth}(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1;
                endfor;
            endif;
            fx(nans) = NaN;
            
            % store the function values
            iv.fx = fx;
            
            % compute the error estimate
            nc = norm(c_new);
            iv.c(1:n(depth)+1,depth) = c_new;
            c_diff = norm(iv.c(:,depth-1) - iv.c(:,depth));
            iv.err = (b-a) * c_diff;
            
            % compute the new integral
            iv.int = (b-a) * w(1) * c_new(1);
            
            % split this interval prematurely?
            if ( nc > 0 && c_diff / nc > 0.1 ), split = 1; endif;
            
        % interval of any other depth
        else
            split = 1;
        endif;
        
        % can we safely ignore this interval?
        % discard if the nodes are below machine resolution or the
        % error estimate falls below the local tolerance
        if ( points(2) <= points(1) || points(end) <= points(end-1) || ...
            iv.err < abs(iv.int) * eps * Vcond(iv.depth) )
            err_final = err_final + iv.err;
            int_final = int_final + iv.int;
            ivals(i_max) = ivals(nr_ivals);
            nr_ivals = nr_ivals - 1;
            
        % if we have to split the interval...
        elseif ( split )
        
            % get the center
            m = (a + b) / 2;
            
            % construct the left interval
            ivl = struct ();
            ivl.a = a; ivl.b = m;
            ivl.depth = 1; ivl.rdepth = iv.rdepth + 1;
            ivl.c = zeros ( n(4)+1 , 4 );
            fx = [ iv.fx(1) ; f( (a+m)/2 + (m-a)*xi{1}(2:end-1)/2 ) ; iv.fx((end+1)/2) ];
            nr_points = nr_points + n(1)-1;
            nans = find ( ~isfinite( fx ) );
            fx(nans) = 0;
            c_new = Vinv{1} * fx;
            if ( length(nans) > 0 )
                b_new = bee{1}; n_new = n(1);
                for i=nans
                    b_new(1:end-1) = (U(1:n(1)+1,1:n(1)+1) - diag(ones(n(1),1)*xi{1}(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1;
                endfor;
            endif;
            fx(nans) = NaN;
            ivl.fx = fx;
            ivl.c(1:n(1)+1,1) = c_new; nc = norm(c_new);
            c_diff = norm(ivl.c(:,1) - T_left * iv.c(:,depth))
            ivl.err = (m - a) * c_diff;
            ivl.int = (m - a) * ivl.c(1,1) * w(1);
            
            % check for divergence
            ivl.ndiv = iv.ndiv + (abs(iv.c(1,1)) > 0 && ivl.c(1,1) / iv.c(1,1) > 2);
            if ( ivl.ndiv > ndiv_max && 2*ivl.ndiv > ivl.rdepth ), warning (sprintf ('possibly divergent integral in the interval [%e,%e]! (h=%e)',a,m,m-a)); int = sign(int) * Inf; return; endif;
            
            % construct the right interval
            ivr = struct ();
            ivr.a = m; ivr.b = b;
            ivr.depth = 1; ivr.rdepth = iv.rdepth + 1;
            ivr.c = zeros(n(4)+1,4);
            fx = [ iv.fx((end+1)/2) ; f((m+b)/2+(b-m)*xi{1}(2:end-1)/2) ; iv.fx(end) ];
            nr_points = nr_points + n(1)-1;
            nans = find ( ~isfinite( fx ) );
            fx(nans) = 0;
            c_new = Vinv{1} * fx;
            if ( length(nans) > 0 )
                b_new = bee{1}; n_new = n(1);
                for i=nans
                    b_new(1:end-1) = (U(1:n(1)+1,1:n(1)+1) - diag(ones(n(1),1)*xi{1}(i),1)) \ b_new(2:end);
                    b_new(end) = 0; c_new = c_new - c_new(n_new+1) / b_new(n_new+1) * b_new(1:end-1);
                    n_new = n_new - 1;
                endfor;
            endif;
            fx(nans) = NaN;
            ivr.fx = fx;
            ivr.c(1:n(1)+1,1) = c_new; nc = norm(c_new);
            c_diff = norm(ivr.c(:,1) - T_right * iv.c(:,depth));
            ivr.err = (b - m) * c_diff;
            ivr.int = (b - m) * ivr.c(1,1) * w(1);
            
            % check for divergence
            ivr.ndiv = iv.ndiv + (abs(iv.c(1,1)) > 0 && ivr.c(1,1) / iv.c(1,1) > 2);
            if ( ivr.ndiv > ndiv_max && 2*ivr.ndiv > ivr.rdepth ), warning (sprintf ('possibly divergent integral in the interval [%e,%e]! (h=%e)',m,b,b-m)); int = sign(int) * Inf; return; endif;
            
            % store the intervals
            ivals(i_max) = ivl;
            nr_ivals = nr_ivals + 1;
            ivals(nr_ivals) = ivr;
            
        % otherwise, just store the interval
        else
        
            ivals(i_max) = iv;
            
        endif;
        
        % compute the running err and new max
        [ dummy , i_max ] = max ( [ ivals(1:nr_ivals).err ] );
        err = err_final + sum( [ ivals(1:nr_ivals).err ] );
        int = int_final + sum( [ ivals(1:nr_ivals).int ] );
        
        % nuke smallest element if stack is larger than 200
        if ( nr_ivals > 200 )
            [ dummy , i_min ] = min ( [ ivals(1:nr_ivals).err ] );
            err_final = err_final + ivals(i_min).err;
            int_final = int_final + ivals(i_min).int;
            ivals(i_min) = ivals(nr_ivals);
            if ( i_max == nr_ivals ), i_max = i_min; endif;
            nr_ivals = nr_ivals - 1;
        endif;
        
        % get up and leave?
        if ( err == 0 || err < abs(int) * tol || (err_final > abs(int) * tol && err - err_final < abs(int) * tol) || nr_ivals == 0 )
            break;
        endif;
        
    endwhile; % main loop
    
    % clean-up and tabulate
    err = err + err_excess;
    
end


%!assert (cquad (@sin,-pi,pi, 1e-6), 0, 1e-6)
%!assert (cquad (inline('sin'),-pi,pi, 1e-6), 0, 1e-6)
%!assert (cquad ('sin',-pi,pi, 1e-6), 0, 1e-6)

%!assert (cquad (@sin,-pi,0, 1e-6), -2, 2e-6)
%!assert (cquad (@sin,0,pi, 1e-6), 2, 2e-6)
%!assert (cquad (@(x) 1./sqrt(x), 0, 1,1e-6), 2, 2e-6)
%!assert (cquad (@(x) abs (1 - x.^2), 0, 2, 1e-6, [ 1 ]), 2, 2e-6)
%!assert (cquad (@(x) 1./(sqrt(x).*(x+1)), 0, Inf, 1e-6), pi, 3.1e-6)

%!assert (cquad (@(x) exp(-x .^ 2), -Inf, Inf, 1e-6), sqrt(pi), 1e-6)

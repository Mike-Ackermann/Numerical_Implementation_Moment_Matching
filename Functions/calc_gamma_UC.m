
function gamma_sig = calc_gamma_UC(s,n,der_order)
%s is interpolation point
%n is order of system
%der_order is order of derivative.
%it is assumed that s is within some delta of the unit circle, so that s^n
%does not blow up.


%if calulating M0, just put s^(0:n) in a vector.
if der_order == 0
    exponents = 0:n;
    gamma_sig = (s.^exponents).';
else
    gamma_sig = ones(n+1,1);
    
    for i = 1:der_order %get coefficients for multiplication
        gamma_sig = polyder(gamma_sig)';
    end
    gamma_sig = flip(gamma_sig); %get into right order
    gamma_sig = [zeros(der_order,1);gamma_sig];
    
    for i = der_order+2:n+1 % calculate gamma_n(sigma)
        gamma_sig(i) = gamma_sig(i)*(s^(i-1-der_order));
    end
    
end
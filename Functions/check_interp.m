function [val,LS] = check_interp(s,U,n,tau)
% NEED TO CHANGE TO MATCH PAPER

%Checks interpolation conditions (Modified conditions from Thm 13 in
%Burohman, Besselink, Scherpen 2020

%%%%% OUTPUTS %%%%%%
%val is a boolean for if we can interpolate (pass uniqueness criterion)
%LS is boolean for if we can only compute a Least Squares solution

%%%%% INPUTS %%%%%%
%U is the orthogonal subspace for the current data window
%n is order of system
%s is the interpolation point

gamma_sig = calc_gamma(s,n,0);
gamma_sig = gamma_sig/norm(gamma_sig);
z = [zeros(n+1,1);gamma_sig];
b = [gamma_sig;zeros(n+1,1)];

v = z-(U*(U'*z));
b_perp = b-(U*(U'*b));

%%%%%%%% UNIQUENESS CONDITION %%%%%%%%%%%%%
% check if the norm of the component of z orthogonal to U is not less than
% tau.  "Not less than" is used to allow tau = nan to effectively skip this
% check.
val = ~(norm(v) < tau);

%%%%%%%% EXISTENCE CONDITION %%%%%%%%%%%%%

%in_range_U = norm(b_perp) < tau;
% b_perp is the part of b orthogonal to range(U).  If b_perp is not
% (numerically) linearly dependent on v we say that b in not in range([U
% b]).  This also captures the case where b is in range(U) since 
% b_perp >= ((v*(v'*b_perp))/(norm(v)^2) >= 0

% if tau = nan, then LS is always 1 (use \, not special solver)
LS = ~(b_perp - ((v*(v'*b_perp))/(norm(v)^2)) < tau);




[~,rank_U] = size(U);
if val
    % if [U z b] is full rank, then we solve a LS problem
    LS = rank([U z b]) == rank_U+2;
else
    % if z is not sufficiently linearly independent of U we cannot compute
    % transfer function value.
    return
end


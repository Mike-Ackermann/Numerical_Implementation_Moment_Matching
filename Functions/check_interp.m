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

val = 0;
LS = 0;

gamma_sig = calc_gamma(s,n,0);
gamma_sig = gamma_sig/norm(gamma_sig);
z = [zeros(n+1,1);gamma_sig];
b = [gamma_sig;zeros(n+1,1)];


%calculate nessesary ranks

%use tau as tolerance for both rank decisions
%rank1 = rank([U z b]);
%rank2 = rank([U z]);
%[~,rank4] = size(U);
%rank3 = rank4+1;

% check if the norm of the componet of z orthogonal to U is not less than
% tau.  "Not less than" is used to allow tau = nan to effectively skip this
% check.
val = ~(norm(z-U*(U'*z)) < tau);

[~,rank_U] = size(U);
if val
    % if [U z b] is full rank, then we solve a LS problem
    LS = rank([U z b]) == rank_U+2;
else
    % if z is not sufficiently linearly independent of U we cannot compute
    % transfer function value.
    return
end


% if ~(rank3 == rank2)
%     %if M0 is not unique, we can do nothing and must refuse this window
%     return
% elseif rank1 == rank2
%     %if b is in Range([U,z]), then we can find M0 and do not need a LS
%     %solution
%     val = 1;
%     LS = 0;
% else
%     %if b is not in Range([U,z]), then we need a LS solution, where M0 is
%     %guarenteed to be unique
%     val = 1;
%     LS = 1;
% end

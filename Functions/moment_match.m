
function [Mj,cond_num,res] = moment_match(s,n,W,k,tau1,tau2)
%Function that performs moment matching.  Original method from
%Burohman, Besselink, Scherpen 2020.  Modified for numerical
%implementation.

%%%%%% INPUTS %%%%%
%s in C is the point we wish to learn H(s)
%W is an orthogonal subspace for the current window of data used to
%calculate the moments.
%n is (approximate) order of the system
%k number of desired moments to match

%%%%% OUTPUTS %%%%%
%Mj are the moments from 0 to k
%cond_num is the condition number of the linear system
%res is the relative residual of the solution to the linear system


%% Check interpolation conditions
% Check if it is possible to find a unique solution for M_0 at s for this
% window
[unique_cond, exist_cond] = check_interp(s,W,n,tau1,tau2);
%if we cant interpolate, quit
if ~unique_cond || ~exist_cond
    Mj = NaN(k+1,1);
    cond_num = NaN;
    res = NaN;
    return
end

%% Calculate 0th moment of system.
%Calculate an estimate for M_0
gamma_sig = calc_gamma(s,n,0);
z = [zeros(n+1,1);gamma_sig];

A = [W -z];
b = [gamma_sig;zeros(n+1,1)];

[Q,R] = orthogonalize_Wz(W,-z);
x = R\(Q'*b);

M0 = x(end);

%Calculate the relative residual
res = norm(A*x - b)/norm(b);

cond_num = cond(A);

%% Calculate higher order moments
% NEED TO ADD CHECKING RANK CONDITIONS
Mj = NaN(k+1,1);
Mj(1) = M0;

for i = 1:k
    %make the right hand side of the linear system (equation 30, lemma 17)
    rhs = Make_RHS(i,Mj(1:i),n,s);
    %if Mk is not unique, we can't calculate it
    if isnan(Mj(i))
        %we cannot calculate M_{k+1} if we can't calculate M_k, so quit
        %execution here.
        return
    end
    x = R\(Q'*rhs);
    Mi = x(end);
    Mj(i+1) = Mi;
end
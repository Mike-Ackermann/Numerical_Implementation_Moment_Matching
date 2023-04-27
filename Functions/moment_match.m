
function [Mj,cond_num,res,LS] = moment_match(s,n,W,k)
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
[val, LS] = check_interp(s,W,n);
%if we cant interpolate, quit
if ~val
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
if ~LS
    [Q,R] = orthogonalize_Wz(W,-z);
    x = R\(Q'*b);
else
    x = A\b;
end

M0 = x(end);

%Calculate the relative residual
res = norm(A*x - b)/norm(b);
if ~LS %if we didnt use least squares, calculate regular condition number
    cond_num = cond(A);
else %if we did use least squares, calculate LS condition number
    %Demmel Thm 3.4
    bhat = A*x;
    arg = abs(bhat'*b)/(norm(bhat)*norm(b));
    if abs(arg) > 1
        if arg < 0
            arg = -1;
        else
            arg = 1;
        end
    end
    theta = acos(arg);

    K2A = cond(A);
    cond_num = (2*K2A/cos(theta))+(tan(theta)*K2A^2);
end
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
    x = A\rhs;
    Mi = x(end);
    Mj(i+1) = Mi;
end
function [Mj,cond_num,res,LS] = moment_match(s,n,W,k)
%Function that performs moment matching.  Original idea taken from
%Burohman, Besselink, Scherpen 2020.  Modified for numerical
%implementation.

%%%%%% INPUTS %%%%%
%s in C is the point we want to calculate moments at.  s should be very
%close to the unit circle, otherwise calculations in calc_gamma_UC will
%cause overflow or numerical instability.
%W is an orthogonal subspace for the current window of data used to
%calculate the moments.
%n is calculated order of the system
%k number of desired moments to match

%%%%% OUTPUTS %%%%%
%Mj are the moments from 0 to k
%cond_num is the condition number of the linear system
%res is the relative residual of the solution to the linear system
%LS is a vector of boolean values for if we had to use a least squares
%solution to calculate the moment.


%% Check interpolation conditions
% Check if it is possible to find a unique solution for M_0 at s for this
% window
[val, LS] = check_interp(s,W,n);
%val = 1;
if ~val %if we cant interpolate, quit
    Mj = NaN(k+1,1);
    cond_num = NaN;
    res = NaN;
    return
end

%% Calculate 0th moment of system.
%Calculate an estimate for M_0
gamma_sig = calc_gamma_UC(s,n,0);
z = [zeros(n+1,1);gamma_sig];

A_true = [W -z];
b = [gamma_sig;zeros(n+1,1)];
%%%%%% Hacky way of checking for rank deficiencies
lastwarn('', '');
%xi = A_true\b;


A = [W -z/norm(z)];
%%%%%%
%[U,S,V] = svd(A,0);
%U1 = U(:,1:14);
%S1 = S(1:14,1:14);
%V1 = V(:,1:14);
%xi = V1*(S1\(U1'*b));
xi = A\b;

xi(end) = xi(end)/norm(z);
%%%%%%%
%%%%%%%
%[Q,R] = qr(A_true);
%xi = R\(Q'*b);
%%%%%%%

%these two seem to have difference
%%%%%%%%%%%%%%%%%%%
%xi = A\(b/norm(b));
%xi(1:end-1) = xi(1:end-1)*norm(z);
%%%%%%%%%%%%%%%%%%%
[~, warnId] = lastwarn();


%if ~isempty(warnId)
    %using QR this also doesnt work
    %b2 = A'*b;
    %xi = (A'*A)\b2;
    
%    [~,q] = size(A);
%    [~,R] = qr([A b],0);
%    Qb = R(1:q,q+1);
%    R = R(1:q,1:q);
%    xi = R\Qb;

%end
%%%%%% 

M0 = xi(end);

%Calculate the residual
res = norm(A_true*xi - b)/norm(b);
threshhold = 10^(-12);
if ~LS %if we didnt use least squares, calculate regular condition number
    cond_num = cond(A);
else %if we did use least squares, calculate LS condition number
    %Demmel Thm 3.4
    bhat = A_true*xi;
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

if LS && res < threshhold
    %If residual is acceptabily small, count it as not LS.
    LS = 0;
end

%% Calculate higher order moments
Mj = NaN(k+1,1);
Mj(1) = M0;

for i = 1:k
    %make the right hand side of the linear system (equation 30, lemma 17)
    rhs = Mk_RHS(i,Mj(1:i),n,s);
    %if Mk is not unique, we can't calculate it
    if isnan(Mj(i))
        %we cannot calculate M_{k+1} if we can't calculate M_k, so quit
        %execution here.
        return
    end
    xi = A\rhs;
    xi(end) = xi(end)/norm(z);
    Mi = xi(end);
    Mj(i+1) = Mi;
end
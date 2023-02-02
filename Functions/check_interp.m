function [val,LS] = check_interp(s,U,n)
%Checks interpolation conditions (Modified conditions from Thm 13 in
%Burohman, Besselink, Scherpen 2020

%%%%% OUTPUTS %%%%%%
%val is a boolean for if we can interpolate (pass uniqueness criterion)
%LS is boolean for if we can only compute a Least Squares solution

%%%%% INPUTS %%%%%%
%U is the orthogonal subspace for teh current data window
%n is order of system
%s is the interpolation point

val = 0;
LS = 0;

gamma_sig = calc_gamma_UC(s,n,0);
gamma_sig = gamma_sig/norm(gamma_sig);
z = [zeros(n+1,1);gamma_sig];
b = [gamma_sig;zeros(n+1,1)];


%calculate nessesary ranks
rank1 = rank([U z b]);
rank2 = rank([U z]);
rank4 = rank(U);
rank3 = rank4+1;

if ~(rank3 == rank2)
    %if M0 is not unique, we can do nothing and must refuse this window
    return
elseif rank1 == rank2
    %if b is in Range([U,z]), then we can find M0 and do not need a LS
    %solution
    val = 1;
    LS = 0;
else
    %if b is not in Range([U,z]), then we need a LS solution, where M0 is
    %guarenteed to be unique
    val = 1;
    LS = 1;
end

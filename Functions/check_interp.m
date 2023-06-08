function [unique_cond,exist_cond] = check_interp(s,U,n,tau1,tau2)
% NEED TO CHANGE TO MATCH PAPER

%Checks interpolation conditions (Modified conditions from Thm 13 in
%Burohman, Besselink, Scherpen 2020

%%%%% OUTPUTS %%%%%%
%unique_cond is a boolean for if we can interpolate (pass uniqueness criterion)
%LS is boolean for if we can only compute a Least Squares solution

%%%%% INPUTS %%%%%%
%U is the orthogonal subspace for the current data window
%n is order of system
%s is the interpolation point

gamma_sig = calc_gamma(s,n,0);
z = [zeros(n+1,1);gamma_sig];
b = [gamma_sig;zeros(n+1,1)];

norm_z = norm(z);

v = z-(U*(U'*z));
b_perp = b-(U*(U'*b));

%%%%%%%% UNIQUENESS CONDITION %%%%%%%%%%%%%
unique_cond = norm(v) > (tau1*norm_z);

%%%%%%%% EXISTENCE CONDITION %%%%%%%%%%%%%
exist_cond = norm(b_perp - ((v*(v'*b_perp))/(norm(v)^2))) < (tau2*norm_z);


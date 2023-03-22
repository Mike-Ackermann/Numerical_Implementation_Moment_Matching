function RHS = Make_RHS(k,Mj,n,s)
%Makes the right hand side of the linear system given in equation 30 from 
%Burohman, Besselink, Scherpen 2020.  Used to calculate transfer function
%derivatives directly from data.

%%%%% INPUTS %%%%%%
%k is order of moment to match
%Mj is vector of moments 0 to k-1
%n is order of system
%w is frequency to match

%%%%%% OUTPUTS %%%%%%%
%RHS is a vector given in equation 30 from Burohman, Besselink, Scherpen
%2020


%calulate the top half of RHS
gamma_sig_1 = calc_gamma(s,n,k);

%calculate bottom half of RHS
gamma_sig_2 = zeros(n+1,1);
for j = 0:k-1
    gamma_sig_2 = gamma_sig_2...
        + nchoosek(k,j)*Mj(j+1)*calc_gamma(s,n,k-j);
end

RHS = [gamma_sig_1;gamma_sig_2];
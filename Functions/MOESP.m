function n = MOESP(U,Y)
%Function to find an approximation to the system order.
%%%%%% INPUTS %%%%%%
%U,Y are FULL data vectors from time domain simulation of system
%%%%%% OUTPUTS %%%%%
%n is (numerical) order of system. Typically gives close to number of
%non-zero Hankle singular values.

n_surr = floor(length(U)/3)-1; 

Hu = HankMat(U,n_surr);
Hy = HankMat(Y,n_surr);
[M,~] = size(Hu);

[~, R] = qr([Hu.' Hy.'],0);
R22 = R(end-M+1:end,end-M+1:end);
n = rank(R22);

%MOESP requirest persistently exciting inputs
persistent_exitate = rank(Hu) > n;
if ~persistent_exitate
    fprintf('Not Persistently Exitated\n')
    n = NaN;
end


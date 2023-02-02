
function n = Test_subspacef(U,Y)
%Function to find an approximation to the system order. Taken from section
%9.3.1 Theorem 9.5 of "Matrix pencils in time and frequency domain system
%identification", Ionita, Antous, 2012.

%%%%%% INPUTS %%%%%%
%U,Y are FULL data vectors from time domain simulation of system
%%%%%% OUTPUTS %%%%%
%n is (numerical) order of system. Typically gives close to number of
%non-zero Hankle singular values.

%could set this to be any value less than length(U)/2
n_surr = floor(length(U)/3)-1; 

Hu = HankMat(U,n_surr);
Hy = HankMat(Y,n_surr);
[M,~] = size(Hu);

[~, R] = qr([Hu.' Hy.'],0);
R22 = R(end-M+1:end,end-M+1:end);
%The rank of the system is the rank of R_{22} (Thm 9.5.1)
n = rank(R22);
%figure; surf(log10(abs(R22(1:100,1:100).'))); shading flat

%%%%%%%%% NOTE R22 %%%%%%%%%%%
%if A in R^{mxn),n = m + 1
% R22 = upper triangle with appended column
%last column = Q'*A(:,end)
%so Q is a full orthogonal basis for R^m, and Q*R(:,end) = A(:,end)


persistent_exitate = rank(Hu) > n;
%the precise statement of the theorem has rank(Hu) > n as an assumption to
%calculate n. Here we check it after calculating n.  In practice, this has
%not shown any issues, except where other numerical issues have already
%been known to occur.
if ~persistent_exitate
    fprintf('Not Persistently Exitated\n')
    n = NaN;
end


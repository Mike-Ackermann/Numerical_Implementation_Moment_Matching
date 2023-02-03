function [Ap,Bp,Cp,Ep] = Loewner(s,Hs,nmax,epsilon)
%MIGHT WANT TO ADD A TOL OPTION TO ADD INFORMATION ABOUT HOW GOOD OUR
%APPROXIMATIONS ARE OF THE TRANSFER FUNCTION (PAGE 67 MOR BOOK)
%computes Hermite Loewner
%s are sample points (needs to be even)
%Hps are Transfer function values at s
%removes redundencies by removing singular values less than 10^-10 (rel)
r = length(s)/2;

%make sure Hs is a column vector
[col,~] = size(Hs);
if col == 1
    Hs = Hs.';
end

s1 = s(1:r); s2 = s(r+1:end);
Hs1 = Hs(1:r); Hs2 = Hs(r+1:end);

Mi = zeros(r); Li = zeros(r);
Bi = Hs1;
Ci = Hs2.';

for i = 1:r
    for j = 1:r
        Li(i,j) = (Hs1(i)-Hs2(j))/(s1(i)-s2(j));
        Mi(i,j) = (s1(i)*Hs1(i)-s2(j)*Hs2(j))/(s1(i)-s2(j));
    end
end

Ei = -Li;
Ai = -Mi;

%% Keep real
%%%%%%%%%%%%%%%%%%%%%%%
%check for repeated shifts to not multiply by transformation matrix
repeated_shift = NaN(r/2,1);
count = 1;
for i = 1:2:r
   if s(i) == s(i+1)
       repeated_shift(count) = 1;
   else
       repeated_shift(count) = 0;
   end
   count = count + 1;
end
idx = 1:r/2;
idx = idx(logical(repeated_shift));

% % repeated_shift = imag(s) == 0;
% % repeated_idx = zeros(ceil(r/2));
% % isolated_idx = zeros(ceil(r/2));
% % for i = 1:2:r
% %     if repeated_shift(i) && repeated_shift(i+1)
% %         repeated_idx(floor(i/2)) = 1;
% %     elseif repeated_shift(i)
% %         isolated_idx(floor(i/2)) = 1;
% %     end
% % end

%put the real values in s first
% real_idx = imag(s) == 0;
% indicies = 1:r;
% real_idx = indicies(real_idx);
% s_imag = s(~real_idx);
% s_real = s(real_idx);
I = eye(2);
T1 = [1 1; 1i -1i];
T1c = repmat({T1},1,r/2);

T1inv = (1/2)*[1 -1i; 1 1i];
% count = 1;
% T = [];
% Tinv = [];
% while count <= r
%     if ismember(count,real_idx)
%         T =blkdiag(T,1);
%         Tinv = blkdiag(Tinv,1);
%         count = count + 1;
%     else
%         T = blkdiag(T,T1);
%         Tinv = blkdiag(Tinv,T1inv);
%         count = count + 2;
%     end
% end
T1invc = repmat({T1inv},1,r/2);
for i = 1:length(idx)
   T1c{idx(i)} = I;
   T1invc{idx(i)} = I;
end

Tinv = blkdiag(T1c{:});
T = blkdiag(T1invc{:});

B = Tinv*Bi;
C = Ci*T;
A = Tinv*Ai*T;
E = Tinv*Ei*T;
L = Tinv*Li*T;
M = Tinv*Mi*T;

E = real(E);
A = real(A);
B = real(B);
C = real(C);
L = real(L);
M = real(M);


%% remove redundencies
[Y,theta1,X1] = svd([L M]);
[Y2,theta2,X] = svd([L;M]);

sVal = diag(theta2);
sValScaled = sVal/sVal(1);
%semilogy(1:r,sValScaled);
%hold on
%semilogy(1:r,10^(-10)*ones(r,1));

%epsilon was 10^(-13)
pHat = find(sValScaled < epsilon,1);
p = min([pHat,nmax]);%truncate at max allowed degree or at tolerence level

%semilogy(p,sValScaled(p));

%only remove reduncencies if they are present
if all(size(p)>0)
    Ep = Y(:,1:p)'*E*X(:,1:p);
    Ap = Y(:,1:p)'*A*X(:,1:p);
    Bp = Y(:,1:p)'*B;
    Cp = C*X(:,1:p);
else
    Ep = E;
    Ap = A;
    Bp = B;
    Cp = C;
end

%Hr = @(z) Cp*((z*Ep-Ap)\Bp);
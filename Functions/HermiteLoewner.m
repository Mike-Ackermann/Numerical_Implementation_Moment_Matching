
function [Ap,Bp,Cp,Ep] = HermiteLoewner(s,Hs,Hps,nmax,epsilon)
%computes Hermite Loewner
%s are sample points
%Hs are Transfer function values at s
%Hps are Transfer function derivative values at s
%removes redundencies by removing singular values less than 10^-10 (rel)
r = length(s);
num_real = 0; %initialize num_real, might be overwritten
is_real = imag(s) == 0;
if any(is_real)
    indicies = 1:r;
    real_idx = indicies(is_real);
    num_real = length(real_idx);
    imag_idx = indicies(~is_real);
    s_imag = s(imag_idx);
    s_real = s(real_idx);
    H_imag = Hs(imag_idx);
    H_real = Hs(real_idx);
    Hp_imag = Hps(imag_idx);
    Hp_real = Hps(real_idx);
    
    s = [s_real;s_imag];
    Hs = [H_real;H_imag];
    Hps = [Hp_real;Hp_imag];
end


[col,~] = size(Hs);
if col == 1
    %fprinf('here')
    Hs = Hs.';
end

Mi = zeros(r); Li = zeros(r);
%B = Hs';
%C = Hs;
Bi = Hs;
Ci = Hs.';

for i = 1:r
    for j = 1:r
        if i ~= j
            if s(i) == s(j)
                Li(i,j) = Hps(i);
                Mi(i,j) = (Hs(i)+s(i)*Hps(i));
            else
                Li(i,j) = (Hs(i)-Hs(j))/(s(i)-s(j));
                Mi(i,j) = (s(i)*Hs(i)-s(j)*Hs(j))/(s(i)-s(j));
            end
        else
            Li(i,j) = Hps(i);
            Mi(i,j) = (Hs(i)+s(i)*Hps(i));
        end
    end
end

Ei = -Li;
Ai = -Mi;



%make everything real

%%%%%%%%%%%%%%%%%%%%%%%


num_imag = r-num_real;
T1inv = [1 1; 1i -1i];
T1 = (1/2)*[1 -1i; 1 1i];
if ~(num_real == 0)
    I = eye(num_real);
    T1invc = repmat({T1inv},1,num_imag/2);
    T1c = repmat({T1},1,num_imag/2);
    Tinv = blkdiag(I,T1invc{:});
    T = blkdiag(I,T1c{:});
else
    T1c = repmat({T1},1,r/2);
    T1invc = repmat({T1inv},1,r/2);
    T = blkdiag(T1c{:});
    Tinv = blkdiag(T1invc{:});
end

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
%A =Ai;
%B = Bi;
%C =Ci;
%E = Ei;
%L = Li;
%M = Mi;

%% remove redundencies
[Y,theta1,X1] = svd([L M]);
[Y2,theta2,X] = svd([L;M]);

sVal = diag(theta2);
sValScaled = sVal/sVal(1);
%semilogy(1:r,sValScaled);
%hold on
%semilogy(1:r,10^(-10)*ones(r,1));

%epsilon = 10^-5; %10^(-10)
pHat = find(sValScaled < epsilon,1);
p = min([pHat,nmax]);%truncate at max allowed degree or at tolerence level
%semilogy(p,sValScaled(p));

%ensure that p stays even
if mod(p,2) == 1
    p = p-1;
end


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
%Ep = real(Ep);
%Ap = real(Ap);
%Bp = real(Bp);
%Cp = real(Cp);

%Hr = @(z) Cp*((z*Ep-Ap)\Bp);
function [Q,R] = orthogonalize_Wz(W,z)
% W is tall and skinny orthogonal matrix
% z is any vector
% returns Q, R such that QR = [Wz]
[m,n] = size([W z]);

R = sparse(diag(ones(n,1)));
Q = [W zeros(m,1)];
v = z;
for i = 1:n-1
    w = W(:,i);
    rin = w'*v;
    R(i,n) = rin;
    v = v-rin*w;
end
R(n,n) = norm(v);
Q(:,end) = v/R(n,n);
%R = sparse(R);
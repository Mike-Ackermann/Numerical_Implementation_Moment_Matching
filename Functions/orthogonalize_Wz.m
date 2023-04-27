function [Q,R] = orthogonalize_Wz(W,z)
% W is tall and skinny orthogonal matrix
% z is any vector
% returns Q, R such that QR = [Wz]
[m,n] = size(W);

u = W'*z;
v = z - (W*u);
norm_v = norm(v);

R = sparse([diag(ones(n,1)), u;zeros(1,n), norm_v]);
Q = [W v/norm_v];






% R = sparse(diag(ones(n,1)));
% Q = [W zeros(m,1)];
% v = z;
% for i = 1:n-1
%     w = W(:,i);
%     rin = w'*v;
%     R(i,n) = rin;
%     v = v-rin*w;
% end
% R(n,n) = norm(v);
% Q(:,end) = v/R(n,n);
%R = sparse(R);
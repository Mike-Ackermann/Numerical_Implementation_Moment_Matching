function Y = runDTSys(A,B,C,D,U,T)
%Simulates the discrete time system given by 
%x' = Ax + Bu
%y  = Cx + Du
%for T time steps and input U.

n = length(A);
t_len = length(T);

x = zeros(n,t_len+1);
Y = zeros(t_len,1);
for t = 1:t_len
    x(:,t+1) = A*x(:,t) + B*U(t);
    Y(t) = C*x(:,t) + D*U(t);
end
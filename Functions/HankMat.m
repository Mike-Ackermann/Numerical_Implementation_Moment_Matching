function Hx = HankMat(x,n)
%constructs Hankle Matrix

%x is either input or output data (U or Y) for current window
%n is calculated dimention of system

T = length(x)-1;

Hx = NaN(n+1,T-n+1);
for j = 0:T-n
    for i = 0:n
        Hx(i+1,j+1) = x(i+j+1);
    end
end
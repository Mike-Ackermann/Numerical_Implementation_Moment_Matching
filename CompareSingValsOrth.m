
% For showing the condition number of orth column basis better than with
% Hankel matricies

load RandImagEx1.mat

fs = 1e3;
Ts = 1/fs;
t_end = 1;

t_eval = 0:Ts:t_end;
T = length(t_eval);
U = randn(T,1);
Y = runDTSys(A,B,C,D,U,t_eval);
w = 0.5;
s = exp(w*1i);

n = length(A);
%% Compute Orthogonal Subspace
t = 3*n;
t_end = 1+t;
U_i = U(1:t_end);
Y_i = Y(1:t_end);

Hu = HankMat(U_i,n);
Hy = HankMat(Y_i,n);
W = orth([Hu;Hy]);
A = [Hu;Hy];
%% Get and plot results

gamma_sig = calc_gamma(s,n,0);
z = [zeros(n+1,1);gamma_sig];

Wz = [W z];
Az = [A z];

%This is better I think
S1 = svd(Az);
S2 = svd(Wz);
no = length(S2);
na = length(S1);
%figure; semilogy(S1,'-*'); hold on; semilogy(S2,'-*')

%% PLot
ColorMat = [0 0.4470 0.7410;...
            0.8500 0.3250 0.0980;...
            0.9290 0.6940 0.1250;...
            0.4940 0.1840 0.5560;...
            0.4660 0.6740 0.1880;...
            0.3010 0.7450 0.9330;...
            0.6350 0.0780 0.1840];
f=figure;
f.Position = [476 445 700 280];
semilogy(S1,'Color',ColorMat(1,:),'LineWidth',2)
hold on
semilogy(S2,'Color',ColorMat(2,:),'LineWidth',2)
semilogy(1,S1(1),'o','Color',ColorMat(1,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(na,S1(na),'o','Color',ColorMat(1,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(1,S2(1),'o','Color',ColorMat(2,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(no,S2(no),'o','Color',ColorMat(2,:),'MarkerSize',10,'LineWidth',1.5)
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 16;
%specify tick location and labels
xlabel('$j$','interpreter','latex','fontsize',25)
%set limits of plot
%labels
ylabel('Singular Value','interpreter','latex','fontsize',20)
legend('$[\mathbf G_{\hat n,k}~\mathbf z]$','$[\mathbf U_k~\mathbf z]$','autoupdate','off','Interpreter','latex')
% xlabel('$r$','interpreter','latex','fontsize',30)
ylim([10^-(17),10^(5)])
xlim([-5,205])

%% Report singular values
fprintf('Max S val, Min S val\n')
fprintf('Orthogonalized max: %e, Orthogonalized min: %e\n',S2(1),S2(end))
fprintf('Not orthogonalized max: %e, Not orthogonalized min: %e\n',S1(1),S1(end))




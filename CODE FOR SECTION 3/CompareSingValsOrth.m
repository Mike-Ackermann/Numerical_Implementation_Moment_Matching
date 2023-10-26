
% For showing the condition number of orth column basis better than with
% Hankel matricies

load RandImagEx1.mat
n = length(A);

T = 3*n;
t_eval = 0:T;

%U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);
%Y = Yt.*(1+(1e-8)*randn(length(Yt),1));
%Y = Yt;
w = 0.5;
s = exp(w*1i);
%% Compute Orthogonal Subspace
Hu = HankMat(U,n);
Hy = HankMat(Y,n);
U_trunc = orth([Hu;Hy]);
[U_full, ~, ~] = svd([Hu;Hy],'econ');
G = [Hu;Hy];
%% Get and plot results

gamma_sig = calc_gamma(s,n,0);
z = [zeros(n+1,1);-gamma_sig];

U_trunc_z = [U_trunc z];
U_full_z = [U_full z];
Gz = [G z];

%This is better I think
S1 = svd(Gz);
S2 = svd(U_trunc_z);
S3 = svd(U_full_z);

n1 = length(S2);
n2 = length(S1);
n3 = length(S3);

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
semilogy(S3,'--','Color',ColorMat(3,:),'LineWidth',2)
semilogy(1,S1(1),'o','Color',ColorMat(1,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(n2,S1(n2),'o','Color',ColorMat(1,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(1,S2(1),'o','Color',ColorMat(2,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(n1,S2(n1),'o','Color',ColorMat(2,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(1,S3(1),'*','Color',ColorMat(3,:),'MarkerSize',10,'LineWidth',1.5)
semilogy(n3,S3(n3),'*','Color',ColorMat(3,:),'MarkerSize',10,'LineWidth',1.5)
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
legend('$[\mathbf G_{n}~\mathbf z(\sigma)]$',...
    '$[\mathbf U_c~\mathbf z(\sigma)]$',...
    '$[\mathbf U~\mathbf z(\sigma)]$',...
    'autoupdate','off','Interpreter','latex','Location','southwest')
% xlabel('$r$','interpreter','latex','fontsize',30)
ylim([10^-(17),10^(5)])
xlim([-5,205])

%% Report singular values
fprintf('-----Condition Numbers------\n')
fprintf('[G_n  z]:  %e\n',S1(1)/S1(end))
fprintf('[U_c  z]:  %e\n',S2(1)/S2(end))
fprintf('[U  z]:    %e\n',S3(1)/S3(end))





load RandImagEx1.mat
T = 1000;
t_eval = 0:T;

rng(12345);
U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);

w = 0.5;
z = exp(1i*w);

I = eye(length(A));
Hz_true = C*((z*I-A)\B);

clear opts
opts.der_order = 0;
opts.num_windows = 20;
opts.num_windows_keep = 10;
opts.tau2 = 10^-10;

n_max = 200;
n_skip = 1;
n_start = 20;
n_vec = n_start:n_skip:n_max;
num_n = length(n_vec);
err_vec = nan(num_n,1);
nstd_vec = nan(num_n,1);
for k = 1:num_n
    opts.n = n_vec(k);
    [Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
    err_vec(k) = abs(Hz-Hz_true)/abs(Hz_true);
    nstd_vec(k) = nstd_Hz;
    opts = rmfield(opts,'W');
end

%%
N = MOESP(U,Y);


f = figure;
f.Position = [476 445 700 280];
semilogy(n_vec,err_vec,'LineWidth',2)
hold on
semilogy(n_vec,nstd_vec,'LineWidth',2)
xline(N,'color','#D95319','LineWidth',2)
xline(100,'LineWidth',2)
legend('$\epsilon_{rel}$','$s_W$' ,'$N$','$n$','interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
%specify tick location and labels
%xticks([1e-4,1e-3,1e-2,1e-1,1])
%xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'})
xlim([20,200])
xlabel('$\tilde n$','interpreter','latex','fontsize',25)
%labels
% ylabel('$|H_1(\sigma)-M_{0}^{(\tilde n)}(\sigma)|/|H_1(\sigma)|$','interpreter','latex','fontsize',20)
ylabel('$\epsilon_{rel}(\tilde n)$ and $s_W(\tilde n)$','interpreter','latex','fontsize',20)
%lgd = legend();
%lgd.Location = 'best';

%% Plot hankel singular values
S = PlotHankSingVals(A,B,C,0);
S_MP = find(S/S(1)<10^(-13),1);
f=figure;
%f.Position = [476 445 700 280];
semilogy(S/S(1),'-*','LineWidth',2)
xline(N,'color','#D95319','LineWidth',2)
ax = gca;

Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
%labels
ylabel('Hankel Singular Values','fontsize',20)
legend('HSV','$\hat n$','interpreter','latex','fontsize',20,Location='southwest')
%lgd.Location = 'best';
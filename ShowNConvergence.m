
load RandImagEx1.mat

fs = 1e3;
Ts = 1/fs;
t_end = 1;

t_eval = 0:Ts:t_end;
T = length(t_eval);
w = .5;
U = multiCos(t_eval,fs,w);
%U = FiltNoise(fs,t_end);
Y = runDTSys(A,B,C,D,U,t_eval);

z = exp(1i*w);

I = eye(length(A));
Hz_true = C*((z*I-A)\B);

clear opts
opts.tol = 10^(-5);
opts.der_order = 0;
opts.num_est = 1;

n_max = 200;
n_skip = 1;
n_start = 10;
n_vec = n_start:n_skip:n_max;
num_n = length(n_vec);
err_vec = nan(num_n,1);
for k = 1:num_n
    opts.n = n_vec(k);
    [Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
    err_vec(k) = abs(Hz-Hz_true)/abs(Hz_true);
end

%%
nhat = MOESP(U,Y);
if isnan(nhat)
    nhat = 83;
end

figure
semilogy(n_vec,err_vec,'LineWidth',2)
hold on
xline(nhat,'color','#D95319','LineWidth',2)
xline(100,'LineWidth',2)
legend('$\epsilon_{rel}$','$\hat n$','$\tilde n = n$','interpreter','latex')

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
xlabel('$m$','interpreter','latex','fontsize',25)
%labels
ylabel('$|H_1(\sigma)-M_{0,m}(\sigma)|/|H_1(\sigma)|$','interpreter','latex','fontsize',20)
%lgd = legend();
%lgd.Location = 'best';
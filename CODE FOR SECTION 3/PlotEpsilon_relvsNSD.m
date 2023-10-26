% To show relationship between standard deviation and error

load RandEx1.mat
rng(239843)
n_true = length(A);
T = 1000;
t_eval = 0:T;

U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
log_min_freq = -3;
freqs = logspace(log_min_freq,log10(.99*pi),num);
z = exp(1i*freqs);

clear opts
opts.num_windows = 20;
opts.num_windows_keep = 10;
opts.n = n_true;

%%
tic
[Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc
%%
I = eye(length(A));
H = @(s) C*((s*I-A)\B);
H_true = zeros(num,1);
parfor i = 1:num
    H_true(i) = H(z(i));
end
%% Calculate and report error
err = abs(Hz(:,1)-H_true);
err2 = norm(err);
err2rel = norm(err)/norm(H_true);

fprintf('2-norm error in TF estimates         : %.5e\n',err2)
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
%plot relative error and cond nums
f=figure;
f.Position = [476 445 700 280];
relerr = abs(Hz(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,1),'LineWidth',2)
legend('$\epsilon_{rel}$','$s_W$','Interpreter','latex')%,'Condtion Numbers')
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlim([10^-3,pi])
xticks([10^-3, 10^-2, 10^-1, 10^0, pi])
xticklabels({'10^{-3}', '10^{-2}', '10^{-1}', '10^0', '\pi'})
xlabel('$\omega$','interpreter','latex','fontsize',25)
ylabel('$s_W$, $\epsilon_{rel}$','interpreter','latex','fontsize',25)
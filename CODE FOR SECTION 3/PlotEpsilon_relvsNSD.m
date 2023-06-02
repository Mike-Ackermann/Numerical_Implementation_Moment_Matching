% To show relationship between standard deviation and error

load RandEx1.mat

fs = 1e3;%1e3?
Ts = 1/fs;
t_end = 1;

[A, B, C, D] = ConvDiscSISO(A,B,C,D,Ts);
t_eval = 0:Ts:t_end;
T = length(t_eval);


U = randn(T,1);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
log_min_freq = -3;
freqs = logspace(log_min_freq,log10(.99*pi),num);
r = 1; % radius of points
z = r*exp(1i*freqs);

n_true = length(A);
clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 10;
opts.n = 100;
%%
tic
[Hz,nstd_Hz,cond_nums,residuals,LS_vec,opts] = CalculateTFVals(U,Y,z,opts);
toc
%%
num = length(z);
sysd = ss(A,B,C,D,Ts);

I = eye(length(A));
H = @(s) C*((s*I-A)\B);
H_true = zeros(num,1);
parfor i = 1:num
    H_true(i) = H(z(i));
end
%if close to eps, just set them equal
H_true(abs(H_true) < 1e-15) = Hz(abs(H_true) < 1e-15);
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
legend('$\epsilon_{rel}$','NSD','Interpreter','latex')%,'Condtion Numbers')
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
ylabel('NSD, $\epsilon_{rel}$','interpreter','latex','fontsize',25)
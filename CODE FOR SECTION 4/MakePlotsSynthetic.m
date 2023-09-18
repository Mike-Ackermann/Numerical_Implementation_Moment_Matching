% Script to make plots for the synthetic random example

load Rand1000.mat
n_true = length(A);

rng(287346238);
T = 1000;
t_eval = 0:T;
U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);
% load Reproduce_MakePlotsSynthetic.mat;

num = 400;
log_min_freq = -2; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
z = exp(1i*freqs);


clear opts
opts.der_order = 1;
opts.num_windows = 20;
opts.num_windows_keep = 10;
opts.tau1 = 10^-10;
opts.tau2 = 10^-10;
opts.skip_condition_numbers = true;

%%
tic
[Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc
%%
num = length(z);

I = eye(n_true);
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
parfor i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end


%% Calculate and report error
err = abs(Hz(:,1)-H_true);
err2 = norm(err);
err2rel = norm(err)/norm(H_true);

%fprintf('2-norm error in TF estimates         : %.5e\n',err2)
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
if opts.der_order == 1
    err_der = abs(Hz(:,2)-Hp_true);
    err2D = norm(err_der);
    err2relD = norm(err_der)/norm(Hp_true);
    fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
end

%plot value estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz(:,1)),'--','LineWidth',2)
legend('$|H(e^{\mathbf i \omega})|$',...
    '$|M_0(e^{\mathbf i \omega})|$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-2),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)
ylabel('Magnitude','Interpreter','latex','FontSize',20)

%plot derivative estimates on top of true
figure;
loglog(freqs,abs(Hp_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz(:,2)),'--','LineWidth',2)
legend('$|H''(e^{\mathbf i \omega})|$',...
    '$|M_1(e^{\mathbf i \omega})|$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-2),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)
ylabel('Magnitude','Interpreter','latex','FontSize',20)

% Plot error vs standard deviation
figure;
relerr = abs(Hz(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,1),'LineWidth',2)
legend('$\epsilon_{rel}$','$s_W$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-2),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)

% Plot error vs standard deviation for derivative
figure;
relerrp = abs(Hz(:,2)-Hp_true)./abs(Hp_true);
loglog(freqs,relerrp,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,2),'LineWidth',2)
legend('$\epsilon_{rel}$','$s_W$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-2),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)



% % Plot the relative derivative
% figure;
% relder = abs(Hp_true)./abs(H_true);
% loglog(freqs,relder,'LineWidth',2)
% %legend('$\frac{H''(e^{\mathbf i \omega})}{H(e^{\mathbf i \omega})}$','Interpreter',...
% %    'latex','Location','northwest')
% xlim([10^(-2),pi])
% ax = gca;
% Default_TW = ax.TickLength;
% Default_LW = ax.LineWidth;
% ax.TickLength = Default_TW * 2;
% ax.LineWidth = Default_LW * 2;
% ax.FontSize = 16;
% xlabel('$\omega$','Interpreter','latex','FontSize',20)
% ylabel('$\frac{H''(e^{\mathbf i \omega})}{H(e^{\mathbf i \omega})}$','Interpreter',...
%     'latex')


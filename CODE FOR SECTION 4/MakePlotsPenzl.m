% Script for generating plots for the Penzl example

load Penzl_disc.mat
n_true = length(A);
rng(9876576)
T = 10000;
t_eval = 0:T;
U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);

 
% load Reproduce_Penzl.mat

num = 140;
log_min_freq = -5; %lowest frequency wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
z = exp(1i*freqs);

% calculate true values for error calculations
I = eye(n_true);
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
parfor i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end
%% Plot for n = nhat
clear opts
opts.der_order = 1;
opts.num_windows = 20;
opts.num_windows_keep = 10;
opts.tau1 = 10^-10;
opts.tau2 = 10^-10;
opts.skip_condition_numbers = true;
fprintf('Plot for nhat\n')

tic
[Hz_nhat,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc

% Calculate and report error
err = abs(Hz_nhat(:,1)-H_true);
err2rel = norm(err)/norm(H_true);
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
err_der = abs(Hz_nhat(:,2)-Hp_true);
err2relD = norm(err_der)/norm(Hp_true);
fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
%% plot
%plot value estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz_nhat(:,1)),'--','LineWidth',2)
legend('True $H(e^{\mathbf i \omega})$',...
    'Recovered $H(e^{\mathbf i \omega})$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(log_min_freq),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)

% Plot error vs standard deviation
figure;
relerr = abs(Hz_nhat(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,1),'LineWidth',2)
legend('$\epsilon_{rel}$','$s_W$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(log_min_freq),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)
%% Upping n to 900
fprintf('Plot for n = 300\n')
clear opts
opts.der_order = 1;
opts.num_windows = 40;
opts.num_windows_keep = 10;
opts.tau1 = 10^-10;
opts.tau2 = 10^-10;
opts.n = 900;
opts.skip_condition_numbers = true;

tic
[Hz_300,nstd_Hz_300,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc

% Calculate and report error
err = abs(Hz_300(:,1)-H_true);
err2rel = norm(err)/norm(H_true);
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)

err_der = abs(Hz_300(:,2)-Hp_true);
err2relD = norm(err_der)/norm(Hp_true);
fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
%% plot
%plot value estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz_300(:,1)),'--','LineWidth',2)
legend('$|H(e^{\mathbf i \omega})|$',...
    '$|M_0(e^{\mathbf i \omega})|$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(log_min_freq),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xticks([1e-5,1e-3,1e-1,pi])
xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
xlabel('$\omega$','Interpreter','latex','FontSize',20)
ylabel('Magnitude')
text(2e-5,2e-1,'(a)','FontSize',30)



%%
% Plot error vs standard deviation
f=figure;
f.Position = [476 445 700 280];
relerr_300 = abs(Hz_300(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr_300,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz_300(:,1),'LineWidth',2)
legend('$\epsilon_{rel}$','$s_W$','Interpreter',...
    'latex','Location','northeast')
xlim([10^(log_min_freq),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%ax.FontSize = 16;
ax.FontSize = 20;
xticks([1e-5,1e-3,1e-1,pi])
xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
yticks([1e-15,1e-2])
yticklabels({'10^{-15}','10^{-2}'})
xlabel('$\omega$','Interpreter','latex','FontSize',24)
text(2e-5,5e-13,'(b)','FontSize',30)



%% Plot the relative derivative
f=figure;
f.Position = [476 445 700 280];
loglog(freqs,abs(Hp_true./H_true),'LineWidth',2)
legend('$\frac{H''(e^{\mathbf i \omega})}{H(e^{\mathbf i \omega})}$',...
    'Interpreter','latex','Location','northeast','FontSize',20)
xlim([10^(log_min_freq),pi])

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 20;
xticks([1e-5,1e-3,1e-1,pi])
xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
yticks([1e0,1e3])
yticklabels({'10^{0}','10^{3}'})
xlabel('$\omega$','Interpreter','latex','FontSize',24)
text(2e-5,5e0,'(c)','FontSize',30)





% %% Up to n = 1006
% % Might take a very long time!
% % can uncomment: 
% % load Penzl1006_data.mat 
% % and comment out call to CalculateTFValues.m to plot quickly
% % for precomputed data
% 
% fprintf('Plot for n = 1006\n')
% clear opts
% opts.der_order = 1;
% opts.num_windows = 20;
% opts.num_windows_keep = 10;
% opts.tau1 = 10^-10;
% opts.tau2 = 10^-10;
% opts.n = 1006;
% 
% % tic
% % [Hz_1006,nstd_Hz_1006,cond_nums_1006,residuals_1006,opts] = CalculateTFVals(U,Y,z,opts);
% % toc
% load Penzl1006_data.mat
% 
% % Calculate and report error
% err = abs(Hz_1006(:,1)-H_true);
% err2rel = norm(err)/norm(H_true);
% fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
% 
% err_der = abs(Hz_1006(:,2)-Hp_true);
% err2relD = norm(err_der)/norm(Hp_true);
% fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
% 
% %% plot
% %plot value estimates on top of true
% figure;
% loglog(freqs,abs(H_true),'LineWidth',2)
% hold on
% loglog(freqs,abs(Hz_1006(:,1)),'--','LineWidth',2)
% legend('True $H(e^{\mathbf i \omega})$',...
%     'Recovered $H(e^{\mathbf i \omega})$','Interpreter',...
%     'latex','Location','northwest')
% xlim([10^(-4),pi])
% ax = gca;
% Default_TW = ax.TickLength;
% Default_LW = ax.LineWidth;
% ax.TickLength = Default_TW * 2;
% ax.LineWidth = Default_LW * 2;
% ax.FontSize = 16;
% xlabel('$\omega$','Interpreter','latex','FontSize',20)
% 
% % Plot error vs standard deviation
% figure;
% relerr_1006 = abs(Hz_1006(:,1)-H_true)./abs(H_true);
% loglog(freqs,relerr_1006,'LineWidth',2)
% hold on
% loglog(freqs,nstd_Hz_1006(:,1),'LineWidth',2)
% legend('$\epsilon_{rel}$','$s_W$','Interpreter',...
%     'latex','Location','northwest')
% xlim([10^(-4),pi])
% ax = gca;
% Default_TW = ax.TickLength;
% Default_LW = ax.LineWidth;
% ax.TickLength = Default_TW * 2;
% ax.LineWidth = Default_LW * 2;
% ax.FontSize = 16;
% xlabel('$\omega$','Interpreter','latex','FontSize',20)

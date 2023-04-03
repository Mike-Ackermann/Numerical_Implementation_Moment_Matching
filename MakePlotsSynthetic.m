% Script to make plots for the synthetic random example

load Rand1000.mat

fs = 1e3;%1e3?
Ts = 1/fs;
t_end = 2;

[A, B, C, D] = ConvDiscSISO(A,B,C,D,Ts);
t_eval = 0:Ts:t_end;
T = length(t_eval);
U = randn(T,1);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 400;
log_min_freq = -2; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
r = 1; % radius of points
z = r*exp(1i*freqs);


clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 1;
opts.num_est = 20;

%%
tic
[Hz,nstd_Hz,cond_nums,residuals,LS_vec,opts] = CalculateTFVals(U,Y,z,opts);
toc
%%
num = length(z);
sysd = ss(A,B,C,D,Ts);

I = eye(length(A));
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
parfor i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end
%if close to eps, just set them equal
H_true(abs(H_true) < 1e-15) = Hz(abs(H_true) < 1e-15);

%% Calculate and report error
err = abs(Hz(:,1)-H_true);
err2 = norm(err);
err2rel = norm(err)/norm(H_true);

fprintf('2-norm error in TF estimates         : %.5e\n',err2)
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
if opts.der_order == 1
    err_der = abs(Hz(:,2)-Hp_true);
    err2D = norm(err_der);
    err2relD = norm(err_der)/norm(Hp_true);
    fprintf('2-norm error in Derivative TF estimates         : %.5e\n',err2D)
    fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
end

%plot value estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz(:,1)),'--','LineWidth',2)
legend('True $|H(e^{\mathbf i \omega})|$',...
    'Recovered $|H(e^{\mathbf i \omega})|$','Interpreter',...
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
legend('True $|H''(e^{\mathbf i \omega})|$',...
    'Recovered $|H''(e^{\mathbf i \omega})|$','Interpreter',...
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
legend('$\epsilon_{rel}$','NSTD','Interpreter',...
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
legend('$\epsilon_{rel}$','NSTD','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-2),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)



%% I dont think this is really good to add to the paper, it doesnt really add anything
% Plot the relative derivative
figure;
relder = abs(Hp_true)./abs(H_true);
loglog(freqs,relder,'LineWidth',2)
%legend('$\frac{H''(e^{\mathbf i \omega})}{H(e^{\mathbf i \omega})}$','Interpreter',...
%    'latex','Location','northwest')
xlim([10^(-2),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)
ylabel('$\frac{H''(e^{\mathbf i \omega})}{H(e^{\mathbf i \omega})}$','Interpreter',...
    'latex')


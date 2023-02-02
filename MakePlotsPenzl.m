% Script for generating plots for the Penzl example

load Penzl_disc.mat

fs = 1e4;%1e3?
Ts = 1/fs;
t_end = 10;

t_eval = 0:Ts:t_end;
T = length(t_eval);
U = FiltNoise(fs,t_end);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
log_min_freq = -4; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
r = 1; % radius of points
z = r*exp(1i*freqs);

clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 10;

%% Plot for n = nhat
tic
[Hz_nhat,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc
num = length(z);
n_true = length(A);
I = eye(n_true);
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
for i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end
H_true(abs(H_true) < 1e-15) = Hz_nhat(abs(H_true) < 1e-15);

% Calculate and report error
err = abs(Hz_nhat(:,1)-H_true);
err2 = norm(err);
err2rel = norm(err)/norm(H_true);

fprintf('2-norm error in TF estimates         : %.5e\n',err2)
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
if opts.der_order == 1
    err_der = abs(Hz_nhat(:,2)-Hp_true);
    err2D = norm(err_der);
    err2relD = norm(err_der)/norm(Hp_true);
    fprintf('2-norm error in Derivative TF estimates         : %.5e\n',err2D)
    fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
end

%plot value estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz_nhat(:,1)),'--','LineWidth',2)
legend('True $H(e^{\mathbf i \omega})$',...
    'Learned $H(e^{\mathbf i \omega})$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-4),pi])
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
legend('$\epsilon_{rel}$','NSTD','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-4),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)
%% Upping n to 300
clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 10;
opts.n = 300;

tic
[Hz_300,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc

num = length(z);
n_true = length(A);
I = eye(n_true);
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
for i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end
%if close to eps, just set them equal
H_true(abs(H_true) < 1e-15) = Hz_300(abs(H_true) < 1e-15);

% Calculate and report error
err = abs(Hz_300(:,1)-H_true);
err2 = norm(err);
err2rel = norm(err)/norm(H_true);

fprintf('2-norm error in TF estimates         : %.5e\n',err2)
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
if opts.der_order == 1
    err_der = abs(Hz_300(:,2)-Hp_true);
    err2D = norm(err_der);
    err2relD = norm(err_der)/norm(Hp_true);
    fprintf('2-norm error in Derivative TF estimates         : %.5e\n',err2D)
    fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
end

%plot value estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz_300(:,1)),'--','LineWidth',2)
legend('True $H(e^{\mathbf i \omega})$',...
    'Learned $H(e^{\mathbf i \omega})$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-4),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)

% Plot error vs standard deviation
figure;
relerr = abs(Hz_300(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,1),'LineWidth',2)
legend('$\epsilon_{rel}$','NSTD','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-4),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)

%% Change input to cosine

freqs = [10^(-3.9),10^(-3.5)];
acurate_num = freqs;
U = multiCos(t_eval,fs,freqs);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
log_min_freq = -4; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
freqs = sort([freqs, acurate_num]);
r = 1; % radius of points
z = r*exp(1i*freqs);

clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 10;
opts.n = 300;
opts.t0 = 6e4;

tic
[Hz_300cos,nstd_Hz_300cos,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc
num = length(z);
n_true = length(A);
I = eye(n_true);
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
for i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end
%if close to eps, just set them equal
H_true(abs(H_true) < 1e-15) = Hz_300cos(abs(H_true) < 1e-15);

% Calculate and report error
err = abs(Hz_300cos(:,1)-H_true);
err2 = norm(err);
err2rel = norm(err)/norm(H_true);

fprintf('2-norm error in TF estimates         : %.5e\n',err2)
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
if opts.der_order == 1
    err_der = abs(Hz_300cos(:,2)-Hp_true);
    err2D = norm(err_der);
    err2relD = norm(err_der)/norm(Hp_true);
    fprintf('2-norm error in Derivative TF estimates         : %.5e\n',err2D)
    fprintf('Relative 2-norm error in Derivative TF estimates: %.5e\n',err2relD)
end

%plot value estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz_300cos(:,1)),'--','LineWidth',2)
legend('True $H(e^{\mathbf i \omega})$',...
    'Learned $H(e^{\mathbf i \omega})$','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-4),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)

% Plot error vs standard deviation
figure;
relerr = abs(Hz_300cos(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz_300cos(:,1),'LineWidth',2)
legend('$\epsilon_{rel}$','NSTD','Interpreter',...
    'latex','Location','northwest')
xlim([10^(-4),pi])
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)

%% Try Cosine convergence

freqs = 10^(-3.75);
U = multiCos(t_eval,fs,freqs);
Y = runDTSys(A,B,C,D,U,t_eval);

r = 1; % radius of points
z = r*exp(1i*freqs);

clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 10;
opts.t0 = 6e4;

n_vec = [100,200,300,400,500,600,700];
num_n = length(n_vec);
Hz_n = nan(num_n,1);
nstd_Hz_n = nan(num_n,1);
tic
for k = 1:num_n
    opts.n = n_vec(k);
    [Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
    Hz_n(k) = Hz;
    nstd_Hz_n(k) = nstd_Hz;
end
toc
num = length(z);
n_true = length(A);
I = eye(n_true);
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
for i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end
%if close to eps, just set them equal
H_true(abs(H_true) < 1e-15) = Hz_300cos(abs(H_true) < 1e-15);

%plot error decay vs n
relerr_n = abs(Hz_n-H_true)/abs(H_true);
figure;
semilogy(n_vec,abs(relerr_n),'-*','LineWidth',2)
hold on
semilogy(n_vec,abs(nstd_Hz_n),'-*','LineWidth',2)
legend('$\epsilon_{rel}$',...
    'NSTD','Interpreter',...
    'latex','Location','northwest')
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$n$','Interpreter','latex','FontSize',20)


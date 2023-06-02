%Just a quick script to plot the error
%test change

load RandEx1.mat

fs = 1e3;%1e3?
Ts = 1/fs;
t_end = 10;

[A, B, C, D] = ConvDiscSISO(A,B,C,D,Ts);
t_eval = 0:Ts:t_end;
T = length(t_eval);

acurate_num = .07;
U = randn(T,1);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
log_min_freq = -5; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
r = 1; % radius of points
z = r*exp(1i*freqs);

n_true = length(A);
clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 10;
opts.n = n_true;
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
for i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end
%if close to eps, just set them equal
H_true(abs(H_true) < 1e-15) = Hz(abs(H_true) < 1e-15);

if opts.der_order == 1
    err_der = abs(Hz(:,2)-Hp_true)./abs(Hp_true);
end
%% Calculate and report error
err = abs(Hz(:,1)-H_true);
err2 = norm(err);
err2rel = norm(err)/norm(H_true);

fprintf('2-norm error in TF estimates         : %.5e\n',err2)
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)

% Plot Condition Number and Derivative
relder = (abs(Hp_true)./abs(H_true));
f = figure;
%f.Position = [476 445 560 280];
f.Position = [476 445 700 280];
loglog(freqs,cond_nums/max(cond_nums),'LineWidth',2)
hold on
loglog(freqs,relder/max(relder),'LineWidth',2)
%need to get the prime
legend('$\kappa_2([\mathbf W\,\mathbf z])$','$\frac{H''_0(\sigma)}{H_0(\sigma)}$','Interpreter','latex')

ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
ax.FontSize = 16;
xlabel('$\omega$','Interpreter','latex','FontSize',20)
ylabel('Normalized Quantity','Interpreter','latex','FontSize',20)
xlim([10^-5,pi])
ylim([10^-4,5])
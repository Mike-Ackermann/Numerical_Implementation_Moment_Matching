%Just a quick script to plot the error
%test change

load RandEx1.mat
n_true = length(A);
T = 3*n_true;
t_eval = 0:T;

U = randn(T+1,1);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
log_min_freq = -5; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
z = exp(1i*freqs);


clear opts
opts.num_windows = 1;
opts.num_windows_keep = 1;
opts.n = n_true;
%%
tic
[Hz,nstd_Hz,cond_nums,residuals,opts] = CalculateTFVals(U,Y,z,opts);
toc
%%
I = eye(n_true);
H = @(s) C*((s*I-A)\B);
Hp = @(s) C*((s*I-A)\(-I*((s*I-A)\B)));
H_true = zeros(num,1);
Hp_true = zeros(num,1);
for i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end

%% Plot

% Plot Condition Number and Derivative
relder = (abs(Hp_true)./abs(H_true));
f = figure;
f.Position = [476 445 700 280];
loglog(freqs,cond_nums/max(cond_nums),'LineWidth',2)
hold on
loglog(freqs,relder/max(relder),'LineWidth',2)
legend('$\kappa_2([\mathbf U\,\mathbf z(\sigma)])$','$\frac{H''_0(\sigma)}{H_0(\sigma)}$','Interpreter','latex')

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
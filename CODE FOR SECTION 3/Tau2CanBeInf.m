load RandEx1.mat

fs = 1e3;%1e3?
Ts = 1/fs;
t_end = 1;

[A, B, C, D] = ConvDiscSISO(A,B,C,D,Ts);
t_eval = 0:Ts:t_end;
T = length(t_eval);

w0 = 1e-2;
U = cos(w0*fs*t_eval);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
max_dist = .0036;
%log_min_freq = -4; %lowest frequency/Ts wanted in frequency range
%freqs = logspace(log_min_freq,log10(.99*pi),num);
%freqs = sort([freqs, w]);
freqs = union(linspace(w0-max_dist,w0,num/2), linspace(w0,w0+max_dist,num/2));

z = exp(1i*freqs);

n_true = length(A);
clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 20;
opts.n = n_true;
opts.tau1 = 1e-10;
opts.tau2 = 1e-10;
%%
[Hz1,nstd_Hz1] = CalculateTFVals(U,Y,z,opts);
%%
opts.tau2 = Inf;
[Hz2,nstd_Hz2] = CalculateTFVals(U,Y,z,opts);
%%
num = length(z);
I = eye(n_true);
H = @(s) C*((s*I-A)\B);
H_true = zeros(num,1);
for i = 1:num
    H_true(i) = H(z(i));
end
relerr1 = abs(H_true-Hz1)./abs(H_true);
relerr2 = abs(H_true-Hz2)./abs(H_true);
f = figure;
f.Position = [476 445 700 280];
semilogy(freqs, relerr2)
hold on
semilogy(freqs, nstd_Hz2,'--')
semilogy(freqs, relerr1,'*')
semilogy(freqs, nstd_Hz1,'*')
xlabel('$\omega$','Interpreter','latex')
legend('$\epsilon_{rel}; \tau_2 = \infty$',...
    '$\overline s_x; \tau_2 = \infty$',...
    '$\epsilon_{rel}; \tau_2 = 10^{-10}$',...
    '$\overline s_x; \tau_2 = 10^{-10}$',...
    'interpreter','latex', 'Location','southwest')
xticks([min(freqs), w0, max(freqs)]);
xticklabels({'$\omega_0 - 3.6\times 10^{-3}$','$\omega_0$','$\omega_0 + 3.6\times 10^{-3}$'})


ax = gca;
% make all lines thicker
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
end
% Make all ticks and lines thicker
ax.TickLength = [.02,.05];
ax.LineWidth = 1;
% Make font size larger
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
%DELETE THIS FILE

%Assumes system matrices A,B,C,D
fs = 1e3;%1e3?
Ts = 1/fs;
t_end = 3;

[A, B, C, D] = ConvDiscSISO(A,B,C,D,Ts);
t_eval = 0:Ts:t_end;
T = length(t_eval);

U = randn(T,1);
Y = runDTSys(A,B,C,D,U,t_eval);

num = 100;
log_min_freq = -4; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
%freqs = sort([freqs, acurate_num]);
r = 1; % radius of points
z = r*exp(1i*freqs);
%w = .1;
%z = exp(1i*w);
%flip the interpolation points
%z = -real(z) + 1i*imag(z);
%freqs = freqs;

n_true = length(A);
clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 0;
opts.num_est = 150;
opts.n = 100;
%opts.t0 = 6e4;
%opts.n = n_true;
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
%plot relative error and cond nums
figure;
relerr = abs(Hz(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'*','LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,1),'*','LineWidth',2)
loglog(freqs(LS_vec), abs(relerr(LS_vec,1)),'*','LineWidth',2)
%loglog(freqs,cond_nums,'LineWidth',2)
legend('Relative Error','Normalized STD','AutoUpdate','off')%,'Condtion Numbers')
%plot where we used LS
%loglog(freqs'.*LS_used, relerr.*LS_used,'*','LineWidth',2)

% ax = gca;
% Default_TW = ax.TickLength;
% Default_LW = ax.LineWidth;
% ax.TickLength = Default_TW * 2;
% ax.LineWidth = Default_LW * 2;
% ax.FontSize = 16;
% xlim([-inf,pi/Ts])
% xticks([10^-2, 10^0, 10^2,pi/Ts])
% xticklabels({'10^{-5}','10^{-3}','10^{-1}','\pi'})
% xlabel('$\theta$','interpreter','latex','fontsize',25)
% ylabel('NSD, $\epsilon_{rel}$','interpreter','latex','fontsize',25)

figure
subplot(2,2,1)
loglog(freqs, abs(H_true),'LineWidth',2)
title('Frequency Response')
subplot(2,2,2)
loglog(freqs, cond_nums,'LineWidth',2)
title('Condition Numbers')
subplot(2,2,3)
loglog(freqs, abs(Hp_true)./abs(H_true),'LineWidth',2)
title('Relative derivative (true)')
subplot(2,2,4)
loglog(freqs,relerr,'LineWidth',2)
title('Relative error')

figure
loglog(freqs, abs(Hp_true)/max(abs(Hp_true)),'LineWidth',2)
hold on
loglog(freqs, cond_nums/max(cond_nums),'LineWidth',2)


if opts.der_order == 1
    figure;
    loglog(freqs,err_der,'LineWidth',2)
    title('Derivative Error')
end


figure
loglog(freqs,cond_nums,'LineWidth',2)
xlabel('frequency')
ylabel('condition number')


%plot estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz(:,1)),'--','LineWidth',2)
legend('True H','Aprox H')

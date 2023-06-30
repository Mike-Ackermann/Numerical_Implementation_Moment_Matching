% Calculates values and derivatives of Heat model and reports error

%I DONT NEED THIS FILE WHEN I PUBLISH THE CODE, EVERYTHING I NEED IS
% IN THE OTHER HEAT MODEL FILE
load heat-disc.mat
A = full(E\A); B = full(E\B); C = full(C);
D = 0;
n_true = length(A);

% T = 1000;
% t_eval = 0:T;
% U = randn(T+1,1);
% Y = runDTSys(A,B,C,D,U,t_eval);
load Reproduce_HeatModel.mat

num = 500;
log_min_freq = -4; %lowest frequency/Ts wanted in frequency range
freqs = logspace(log_min_freq,log10(.99*pi),num);
%%%%%%%%%%%%%%%%%%%%
% %Inverse logspace
% delta = diff(freqs);
% delta = fliplr(delta);
% freqs2 = [10^log_min_freq, 10^(log_min_freq)+cumsum(delta)];
% freqs = freqs2;
%%%%%%%%%%%%%%%%%%%%
z = exp(1i*freqs);

clear opts
opts.der_order = 1;
opts.num_windows = 20;
opts.num_windows_keep = 10;
opts.tau1 = 10^-10;
opts.tau2 = 10^-10;
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
parfor i = 1:num
    H_true(i) = H(z(i));
    Hp_true(i) = Hp(z(i));
end

%% Calculate and report error
err = abs(Hz(:,1)-H_true);
err2rel = norm(err)/norm(H_true);
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)

err_der = abs(Hz(:,2)-Hp_true);
err2rel_der = norm(err_der)/norm(Hp_true);
fprintf('Relative 2-norm error in TF derivatives: %.5e\n',err2rel_der)



figure;
relerr = abs(Hz(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,1),'LineWidth',2)
legend('$\epsilon_{rel}$','$s_W$','AutoUpdate','off','Interpreter',...
    'latex')
title('Value Error')


figure
loglog(freqs,abs(Hp_true),'LineWidth',2)
hold on
loglog(freqs, abs(Hz(:,2)),'LineWidth',2)
title('Derivative Response')
legend('True','Recovered')

figure;
loglog(freqs,abs(err_der)./abs(Hp_true),'LineWidth',2)
title('Derivative Error')



%plot estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz(:,1)),'--','LineWidth',2)
legend('True H','Aprox H')
title('Frequency Response')
% Calculates values and derivatives of Heat model and reports error

load heat-disc.mat
A = full(E\A); B = full(E\B); C = full(C);
D = 0;

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
%%%%%%%%%%%%%%%%%%%%
% %Inverse logspace
% delta = diff(freqs);
% delta = fliplr(delta);
% freqs2 = [10^log_min_freq, 10^(log_min_freq)+cumsum(delta)];
% freqs = freqs2;
%%%%%%%%%%%%%%%%%%%%
r = 1; % radius of points
z = r*exp(1i*freqs);

n_true = length(A);
clear opts
opts.tol = 10^(-1);
opts.noise = false;
opts.der_order = 1;
opts.num_est = 10;
%opts.n = 300;
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
err2rel = norm(err)/norm(H_true);
fprintf('Relative 2-norm error in TF estimates: %.5e\n',err2rel)
if opts.der_order == 1
    err_der = abs(Hz(:,2)-Hp_true);
    err2rel_der = norm(err_der)/norm(Hp_true);
    fprintf('Relative 2-norm error in TF derivatives: %.5e\n',err2rel_der)
end


figure;
relerr = abs(Hz(:,1)-H_true)./abs(H_true);
loglog(freqs,relerr,'LineWidth',2)
hold on
loglog(freqs,nstd_Hz(:,1),'LineWidth',2)
loglog(freqs(LS_vec), abs(relerr(LS_vec,1)),'*','LineWidth',2)
%loglog(freqs,cond_nums,'LineWidth',2)
legend('Relative Error','Normalized STD','AutoUpdate','off')%,'Condtion Numbers')
title('Value Error')
%plot where we used LS
%loglog(freqs'.*LS_used, relerr.*LS_used,'*','LineWidth',2)

if opts.der_order == 1
    figure
    loglog(freqs,abs(Hp_true),'LineWidth',2)
    hold on
    loglog(freqs, abs(Hz(:,2)),'LineWidth',2)
    title('Derivative Response')
    legend('True','Recovered')

    figure;
    loglog(freqs,abs(err_der)./abs(Hp_true),'LineWidth',2)
    title('Derivative Error')
end


%plot estimates on top of true
figure;
loglog(freqs,abs(H_true),'LineWidth',2)
hold on
loglog(freqs,abs(Hz(:,1)),'--','LineWidth',2)
legend('True H','Aprox H')
title('Frequency Response')
% show that we have similar quality estimates from only a very low number
% of windows

%% construct data to use to find estimates of transfer function values
load RandEx1.mat

fs = 1e4;
Ts = 1/fs;
t_end = 1;
t_eval = 0:Ts:t_end;
T = length(t_eval);
%U = randn(T,1);
%Y = runDTSys(A,B,C,D,U,t_eval);
load Data_Files/ReproduceLowWindowOK.mat

%% set up
num_est_vec = 25:50:200;
num_est_vec = [10,20,30,70,150];
num_trial = length(num_est_vec);
H_true = @(s) C*((s*eye(length(A))-A)\B);

w = .01;
z = exp(1i*w);
n = 100;
norm_std_vec = NaN(num_trial,1);
Hs_vec = NaN(num_trial,1);

clear opts
opts.tol = 10^(-1);
opts.der_order = 0;
opts.n = 100;
for i = 1:num_trial
    opts.num_est = num_est_vec(i);
    [Hz,nstd_Hz,cond_nums,residuals] = CalculateTFVals(U,Y,z,opts);
    norm_std_vec(i) = abs(nstd_Hz*Hz); %unnormalize
    Hs_vec(i) = Hz;
end

%% plot
load ColorMat.mat

param = exp(1i*2*pi*(0:100)/100);

figure;
axis equal
abs_err_vec = NaN(num_trial,1);
leg = cell(num_trial+1,1);
Hs_true = H_true(z);
plot(Hs_true,'k.','MarkerSize',30)
leg{1} = '$H_0(\sigma)$';
hold on
plot(norm_std_vec(1)*param + Hs_vec(1),'Color',ColorMat(1,:),'LineWidth',2)
abs_err_vec(1) = abs(Hs_true-Hs_vec(1));
leg{2} = strcat('$n_w = $',num2str(num_est_vec(1)));
for i = 2:num_trial
    plot(norm_std_vec(i)*param + Hs_vec(i),'Color',ColorMat(i,:),'LineWidth',2)
    abs_err_vec(i) = abs(Hs_true-Hs_vec(i));
    leg{i+1} = strcat('$n_w = $',num2str(num_est_vec(i)));
end
for i = 1:num_trial
    plot(Hs_vec(i),'*','Color',ColorMat(i,:),'MarkerSize',10)
end


legend(leg,'Interpreter','latex')

max_abs_err = .5*10^(floor(log10(max(abs_err_vec)))+1);
max_abs_err_str = num2str(log10(max_abs_err));
ax = gca;
Default_TW = ax.TickLength;
Default_LW = ax.LineWidth;
%double tick width and length
ax.TickLength = Default_TW * 2;
ax.LineWidth = Default_LW * 2;
%change font size
ax.FontSize = 14;
%specify tick location and labels
xlim([real(Hs_true) - max_abs_err,real(Hs_true) + max_abs_err])
xticks([real(Hs_true) - max_abs_err,real(Hs_true),real(Hs_true) + max_abs_err])
xticklabels({strcat('$-10^{',max_abs_err_str,'}$'),'$0$',strcat('$10^{',max_abs_err_str,'}$')}')
ax.TickLabelInterpreter='latex';
xlabel('$Re(M_0)-Re(H_0(\sigma))$','interpreter','latex','fontsize',20)
ylabel('$Im(M_0)-Im(H_0(\sigma))$','interpreter','latex','fontsize',20)
%set limits of plot
ylim([imag(Hs_true) - max_abs_err,imag(Hs_true) + max_abs_err])
yticks([imag(Hs_true) - max_abs_err,imag(Hs_true),imag(Hs_true) + max_abs_err])
yticklabels({strcat('$-10^{',max_abs_err_str,'}$'),'$0$',strcat('$10^{',max_abs_err_str,'}$')}')


xlim([real(Hs_true) - max_abs_err,real(Hs_true) + max_abs_err])
xticks([real(Hs_true) - max_abs_err,real(Hs_true),real(Hs_true) + max_abs_err])
xticklabels({strcat('$-10^{',max_abs_err_str,'}$'),'$0$',strcat('$10^{',max_abs_err_str,'}$')}')
%labels
lgd = legend();
lgd.Location = 'EastOutside';
%axis equal





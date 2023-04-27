
function [Hs,n_std_Hs,cond_nums,residuals,LS_vec,opts] = ...
    CalculateTFVals(U,Y,z_vec,opts)
% Finds values (and derivatives) of underlying dynamical system that
% produces output Y given input U.

%%%%%%%%% INPUTS %%%%%%%%%
% U is vector of time domain input
% Y is vector of time domain output corresponding to U
% Z_VEC is vector of complex number where we would like to learn frequency
%      informaion
% OPTS.TOL is maximum value for normalized standard deviation of estimates
% OPTS.DER_ORDER is number of derivatives to match
% OPTS.NUM_EST is number of estimates of frequency information to calculate
% OPTS.N is order of underlying dynamical system
% OPTS.T0 is first component in U and Y to use to calculate moments.
%         Increase to reduce influence of transient response

%%%%%%%%% OUTPUTS %%%%%%%%%
% HS are discovered frequency information
% N_STD_HS is normalized standard deviation in the estimates for components
%          of HS
% COND_NUMS are condition numbers of the matrix used to calculate values of
%          transfer function of underlying dynamical system
% RESIDUALS is vector of average residuals from linear system used to
%           calculate values of transfer function of underlying dynamical
%           system

def.tol = 10^(-3);
def.der_order = 0;
def.num_est = 10;
def.n = nan;
def.t0 = 1;
def.W = nan;
if nargin<4
    opts=def;
else
    %Merge default values into opts
    A=fieldnames(def);
    for m=1:length(A)
        if ~isfield(opts,A(m))
            dum=char(A(m));
            dum2=getfield(def,dum);
            opts=setfield(opts,dum,dum2);
        end
    end
end
% if n is not supplied, discover nhat using MOESP
if isnan(opts.n)
    if length(U) > 1000
        n = MOESP(U(1:1000),Y(1:1000));
    else
        n = MOESP(U,Y);
    end
else
    n = opts.n;
end
fprintf('Using n = %d\n',n)
num_est = opts.num_est;
t = 3*n;
T = length(U);


%% Precompute Orthogonal Subspaces
if iscell(opts.W)
    W_cell = opts.W;
else
    num = length(z_vec);
    t0 = opts.t0;
    num_skip = floor((T-t-t0)/num_est);
    %cell array for storing the orthogonal bases for each window
    W_cell = cell(num_est,1);
    count = 1;
    for t_start = t0:num_skip:T-t
        t_end = t_start + t;
        U_i = U(t_start:t_end);
        Y_i = Y(t_start:t_end);
    
        Hu = HankMat(U_i,n);
        Hy = HankMat(Y_i,n);
        W = orth([Hu;Hy]);
        W_cell{count} = W;
        count = count + 1;
    end
end

%% Actually calculate the transfer function values and derivatives
der_order = opts.der_order;
Hs = zeros(num,1 + der_order); %stores accepted TF vals and derivatives
n_std_Hs = zeros(num,1 + der_order); %stores average normalized standard devs
cond_nums = zeros(num,1); %stores average condition numbers of LS problem (val only)
residuals = zeros(num,1); %stores the average residual of the LS problem
LS_vec = false(num,1); %stores if we could solve linear system "exactly"
for k = 1:num
    z = z_vec(k);
    Mj_est_vec = zeros(num_est,1+der_order); %stores all estimates of val and ders
    cond_vec = zeros(num_est,1); %stores all condition numbers
    res_vec = zeros(num_est,1); %Least Squares residual
    LS_vec_temp = zeros(num_est,1);
    for current_SS = 1:num_est
        W = W_cell{current_SS}; %get precomputed subspace
        %calculate and store estimates, cond nums, and LS used data
        [Mj, cond_num,res,LS] = moment_match(z,n,W,der_order);
        Mj_est_vec(current_SS,:) = Mj.';
        cond_vec(current_SS) = cond_num;
        res_vec(current_SS) = res;
        LS_vec_temp(current_SS) = LS;
    end

    for i = 0:der_order
        %take out all instances where we failed to calculate a moment
        idx = ~isnan(Mj_est_vec(:,i+1));
        thin_data = Mj_est_vec(idx,i+1);
        thin_cond = cond_vec(idx);
        thin_res = res_vec(idx);
        %thin_LS = LS_vec_temp(idx);
        if isempty(thin_data)
            thin_data = NaN;
            thin_cond = NaN;
        else
            min_accuracy = opts.tol;
            num_to_accept = min(10,num_est);
            [sorted_thin_res,accepted_res] = sort(thin_res);
            %only lowest 10 res accepted
            accepted_res = accepted_res(1:num_to_accept);
            %reject = 0;
            if sorted_thin_res(num_to_accept) > min_accuracy
                accepted_res = false(length(thin_data),1);
            end
            
            thin_data = thin_data(accepted_res);
            thin_res = thin_res(accepted_res);
            thin_cond = thin_cond(accepted_res);
            if i == 0;
                LS_vec_temp = LS_vec_temp(accepted_res);
            end
        end
        if isempty(thin_data)
            avg_Mj = NaN;
            avg_cond = NaN;
        else
            %set the accepted M0 to be the mean of the calculated M0 values            
            avg_Mj = mean(thin_data);
            std_dev = std(thin_data);
            %normalize the standard deviation
            std_norm = std_dev/abs(avg_Mj);
            if abs(avg_Mj) < 1e-14 %if moment is too small standard deviation error bound doesnt work
                std_norm = 0;
            end
            avg_cond = mean(thin_cond);
            avg_res = mean(thin_res);
            cond_nums(k) = avg_cond;
            residuals(k) = avg_res;
            LS_vec(k) = any(LS_vec_temp); 
            Hs(k,i+1) = avg_Mj;
            n_std_Hs(k,i+1) = std_norm;
        end      
    end
end


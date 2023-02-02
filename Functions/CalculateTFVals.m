function [Hs,n_std_Hs,cond_nums,residuals,opts] = ...
    CalculateTFVals(U,Y,z_vec,opts)
% Finds num transfer function values evenly spaced around the unit circle.
%TESTING FINDING DERIVATIVE INFORMATION

%sigmas are the interpolation points
%Hs are the transfer function values
%n_std_HS is normalized standard deviation of HS

%w0 is vector of frequencies to use for determining n hat
%U is input data
%Y is output data
%tol is an error tolerence for the standard deviation of Hs
%num is number of transfer function values to find
%Ts is sampling frequency of discrete time system
%der_order is number of derivatives to calculate

def.tol = 10^(-3);
def.noise = false;
def.der_order = 0;
def.num_est = 10;
def.n = nan;
def.t0 = 1;
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
        n = Test_subspacef(U(1:1000),Y(1:1000));
    else
        n = Test_subspacef(U,Y);
    end
else
    n = opts.n;
end
fprintf('Using n = %d\n',n)
num_est = opts.num_est;
der_order = opts.der_order;
t = 3*n;
T = length(U);


%% Precompute Orthogonal Subspaces

num = length(z_vec);
t0 = opts.t0;
num_skip = floor((T-t-t0)/num_est);%only calculate num_est values of transfer function for each sigma
%cell array for storing the orthogonal bases for each window
W_cell = cell(num_est,1);
window_percent = NaN(num_est,1);
count = 1;
for t_start = t0:num_skip:T-t
    t_end = t_start + t;
    U_i = U(t_start:t_end);
    Y_i = Y(t_start:t_end);

    Hu = HankMat(U_i,n);
    Hy = HankMat(Y_i,n);
    W = orth([Hu;Hy]);
    %[W,~,~] = svd([Hu;Hy]);
    %W = W(:,1:end-1);
    %W = W(:,1:13);
    W_cell{count} = W;
    %calculate percentage of possible columns we keep
    [~,window_cols] = size(W);
    window_percent(count) = window_cols/(2*n+1);
    count = count + 1;
end
%% Caclulate containment gap between subspaces
containment_gaps = NaN(num_est);
num_cols = NaN(num_est,1);
for i = 1:num_est
    W_i = W_cell{i};
    [Wmi,Wni] = size(W_i);
    for j = 1:num_est
        W_j = W_cell{j};
        [~,Wnj] = size(W_j); 
        if Wnj == Wni
            containment_gaps(i,j) = norm((eye(Wmi)-W_i*W_i')*W_j*W_j');
        end
    end
    num_cols(i) = Wni;
end



%% Actually calculate the transfer function values and derivatives

Hs = zeros(num,1 + der_order); %stores accepted TF vals and derivatives
n_std_Hs = zeros(num,1 + der_order); %stores average normalized standard devs
cond_nums = zeros(num,1); %stores average condition numbers of LS problem (val only)
residuals = zeros(num,1); %stores the average residual of the LS problem
LS_used = zeros(num,1); %stores if we had to use LS
for k = 1:num
    z = z_vec(k);
    Mj_est_vec = zeros(num_est,1+der_order); %stores all estimates of val and ders
    cond_vec = zeros(num_est,1); %stores all condition numbers
    LS_vec = zeros(num_est,1); %stores which times LS was used
    res_vec = zeros(num_est,1); %Least Squares residual
    %parellelized for loop, change to regular "for" for debuging
    for current_SS = 1:num_est
        W = W_cell{current_SS}; %get precomputed subspace
        %calculate and store estimates, cond nums, and LS used data
        [Mj, cond_num,res,LS] = moment_match(z,n,W,der_order);
        Mj_est_vec(current_SS,:) = Mj.';

        cond_vec(current_SS) = cond_num;
        
        LS_vec(current_SS) = LS;

        res_vec(current_SS) = res;
    end

    for i = 0:der_order
        %take out all instances where we failed to calculate a moment
        idx = ~isnan(Mj_est_vec(:,i+1));
        thin_data = Mj_est_vec(idx,i+1);
        thin_LS = LS_vec(idx);
        thin_res = res_vec(idx);
        %if thin_data is empty set std_dev to unacceptable value, set Mj
        %and condition number to NaN
        if isempty(thin_data)
            std_norm = tol+10;
            avg_Mj = NaN;
            thin_cond = NaN;
        %need more than 1 data entry to compute std
        %elseif length(thin_data) == 1
        %    std_norm = tol + 20;
        %    avg_Mj = NaN;
        %    thin_cond = NaN;
        %else, we can calculate a moment!
        else
            num_not_use_LS = 5;
            %set the accepted M0 to be the mean of the calculated M0 values
            if (length(thin_LS) - sum(thin_LS)) >= num_not_use_LS
                %if we have less than num_not_use_LS instances of not
                %needing to use LS, then acceptthe LS estimates.  If we
                %have >= num_not_use_LS "exact" M0 calculations, only use those
                idx = ~thin_LS;
                thin_data = thin_data(idx);
                LS_used(k) = 0;
                
                %have to thin the cond vec here to get them to be the same
                %length
                thin_cond = cond_vec(~isnan(cond_vec));
                thin_cond = thin_cond(idx);
            else
                %if we do not have num_not_use_LS estimates that didn't use
                %LS, then take the num_accept estimates with lowest LS
                %residual.  If the worst estimate has greater than
                %10^min_accuracy residual, then mark as unable to calculate
                %the moment.
                LS_used(k) = 1;
                min_accuracy = opts.tol;
                num_to_accept = min(10,num_est);
                [sorted_thin_res,accepted_res] = sort(thin_res);
                %only lowest 10 res accepted
                accepted_res = accepted_res(1:num_to_accept);
                %reject = 0;
                if sorted_thin_res(num_to_accept) > min_accuracy
                    accepted_res = false(length(thin_LS),1);
                end
                
                thin_data = thin_data(accepted_res);%,i+1); %take out all 0s
                thin_LS = thin_LS(accepted_res);
                thin_res = thin_res(accepted_res);
                thin_cond = cond_vec(accepted_res);
                
            end
            
            avg_Mj = mean(thin_data);
            std_dev = std(thin_data);
            %normalize the standard deviation
            std_norm = std_dev/abs(avg_Mj);
            if abs(avg_Mj) < 1e-14 %if moment is too small standard deviation error bound doesnt work
                std_norm = 0;
            end
        end
        %only calculating the condition number for the interpolation step
        %so just one column
        
        avg_cond = mean(thin_cond);
        avg_res = mean(thin_res);
        cond_nums(k) = avg_cond;
        residuals(k) = avg_res;
        Hs(k,i+1) = avg_Mj;
        n_std_Hs(k,i+1) = std_norm;
    end

end
fprintf('\n')


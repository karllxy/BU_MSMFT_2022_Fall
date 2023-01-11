% MF731 HW4
% Edited by Xuyang Liu

% I refer to the matlab code using on the lecture and do some modifications.

% Global Variables
alpha = .95;
M = 100;
lambda = .97;
theta = .97;
K = 10;

data_file = 'MSFT_AAPL_Log_Returns.csv';
return_data = csvread(data_file);
N = length(return_data(:,1));
log_ret_data = return_data(:, 2:3);
mkt_caps = [448.77 575.11];
port_value = double(1000000);
dollar_pos = port_value * mkt_caps / sum(mkt_caps);
Sigma_hat = cov(log_ret_data(1:M,:));
mu_hat = mean(log_ret_data(1:M,:));

%Run the EWMA to get mean, covariance as of t=9/1/16. 
for i = M+1:N
    %obtain new sigma estimate
    Sigma_hat = theta*Sigma_hat + (1-theta)*transpose(log_ret_data(i,:)-mu_hat)...
        *(log_ret_data(i,:)-mu_hat);
    %obtain new mu estimate
    mu_hat = lambda*mu_hat+(1-lambda)*log_ret_data(i,:);
end

%Estimate the linearized loss portfolio VaR over [t,t+\Delta], given by
%   VaR(\alpha)(L') = -\theta'\mu + sqrt(\theta'\Sigma\theta)N^{-1}(\alpha)
th_port_VaR = -dot(dollar_pos,mu_hat) + (dollar_pos*Sigma_hat*transpose(dollar_pos))^(.5)...
    *norminv(alpha,0,1);

%Use the square root of time rule to project out K days.
th_port_VaR_Kday = (K)^(.5)*th_port_VaR;

%Use the square root of time rule to project out K days for the regulatory
%capital requirement
th_port_VaR_reg = 3*(K)^(.5)*th_port_VaR;


%Number of simlations.
num_sim = 100000;

shock_Kday_port_loss = zeros(num_sim,1);
shock_lastdat_port_loss = zeros(num_sim,1);
shock_sum_stock_return = zeros(num_sim,2);
shock_Kday_VaR_exceed = zeros(num_sim,1);
shock_Kday_VaR_reg_exceed = zeros(num_sim,1);
shock_lastday_VaR_exceed = zeros(num_sim,1);
%--------------------
Kday_port_loss = zeros(num_sim,1);
lastdat_port_loss = zeros(num_sim,1);
sum_stock_return = zeros(num_sim,2);
Kday_VaR_exceed = zeros(num_sim,1);
Kday_VaR_reg_exceed = zeros(num_sim,1);
lastday_VaR_exceed = zeros(num_sim,1);


for n = 1:num_sim
    shock_sim_log_ret = zeros(K,2);
    
    shock_sim_mu_hat = mu_hat;
    shock_sim_Sigma_hat = Sigma_hat;
    
    %assume a -5 sigma shock to the second log return
    shock_sim_log_ret(1,2) = shock_sim_mu_hat(1,2) ...
        - 5*(shock_sim_Sigma_hat(2,2))^(.5);
    %obtain first return using the conditional normal result
    shock_sim_correlation = shock_sim_Sigma_hat(1,2) / (shock_sim_Sigma_hat(1,1) ...
        *shock_sim_Sigma_hat(2,2))^(.5);
    shock_sim_log_ret(1,1) = shock_sim_mu_hat(1,1) + shock_sim_correlation * ...
        (shock_sim_Sigma_hat(1,1)/shock_sim_Sigma_hat(2,2))^(.5) * ...
        (shock_sim_log_ret(1,2)-shock_sim_mu_hat(1,2)) - 5 * (shock_sim_Sigma_hat(1,1)...
        *(1-shock_sim_correlation^2))^(.5); % the difference from lecture
    
    %get returns over days 2...K.
    for k = 2:K
        %update mu,Sigma by EWMA
        shock_sim_Sigma_hat = theta*shock_sim_Sigma_hat + ...
            (1-theta)*transpose(shock_sim_log_ret(k-1,:)-shock_sim_mu_hat)...
            *(shock_sim_log_ret(k-1,:)-shock_sim_mu_hat);
        shock_sim_mu_hat = lambda*shock_sim_mu_hat + (1-lambda)*shock_sim_log_ret(k-1,:);
        
        
        %sample return off normal distribution with new mean and variance.
        shock_sim_log_ret(k,:) = mvnrnd(shock_sim_mu_hat,shock_sim_Sigma_hat,1);
       
    end
    
    %linearized loss = -theta^T \sum_{k=1}^K return_k
    shock_sum_stock_return(n,:) = sum(shock_sim_log_ret,1);
    shock_Kday_port_loss(n,1) = -dot(dollar_pos,shock_sum_stock_return(n,:));
   
    %identify separately loss on last day
    shock_lastday_port_loss(n,1) = -dot(dollar_pos, shock_sim_log_ret(K,:));
    shock_Kday_VaR_exceed(n,1) = (shock_Kday_port_loss(n,1)>th_port_VaR_Kday);
    shock_Kday_VaR_reg_exceed(n,1) = (shock_Kday_port_loss(n,1)>th_port_VaR_reg);
    shock_lastday_VaR_exceed(n,1) = (shock_lastday_port_loss(n,1) > th_port_VaR);
end

%sort Kday losses and last day losses for shocked and un-shocked cases.
shock_sort_Kday_port_loss = sort(shock_Kday_port_loss(:,1),'ascend');
shock_sort_lastday_port_loss = sort(shock_lastday_port_loss(:,1),'ascend');


shock_avg_port_loss = mean(shock_Kday_port_loss(:,1));
shock_port_VaR_Kday=shock_sort_Kday_port_loss(ceil(num_sim*alpha),1);
shock_port_VaR_lastday = shock_sort_lastday_port_loss(ceil(num_sim*alpha),1);
shock_num_Kday_VaR_exceed = sum(shock_Kday_VaR_exceed(:,1));
shock_num_Kday_VaR_reg_exceed = sum(shock_Kday_VaR_reg_exceed(:,1));
shock_num_lastday_VaR_exceed = sum(shock_lastday_VaR_exceed(:,1));

%print the results

fprintf('Confidence: %0.2f\n', alpha);
fprintf('Number of days: %i\n', K);
fprintf('Initial one day VaR: %10.2f\n', th_port_VaR);
fprintf('Initial %i day VaR: %10.2f\n', K, th_port_VaR_Kday);
fprintf('3x Initial %i day VaR: %10.2f\n', K, th_port_VaR_reg);
fprintf('One day VaR after %i days (shock): %10.2f\n',...
    K,shock_port_VaR_lastday);
fprintf('%i day VaR (shock): %10.2f\n', ...
    K, shock_port_VaR_Kday);
fprintf('Average %i day loss: (shock) %10.2f\n',...
    K, shock_avg_port_loss);
fprintf('Pct exceedances for %i day VaR over inital %i day VaR (shock): %10.2f\n',...
    K,K,100*double(shock_num_Kday_VaR_exceed)/num_sim);
fprintf('Pct exceedances for %i day VaR over 3x inital %i day VaR (shock): %10.3f\n',...
    K,K,100*double(shock_num_Kday_VaR_reg_exceed)/num_sim);
fprintf('\n');


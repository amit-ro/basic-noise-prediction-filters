close all
clear all

% noise process and ideal estimator decleration and calculation
fs = 48000;
alpha = 0.9;
sigma_n_sq = 0.5;
prediction_order = 4;
N = sqrt(sigma_n_sq)*randn(fs*10,1);
gauss = randn(fs*10,1);
x = filter(1,[1,-alpha],gauss);
z = x + N;
mu = [0.001 0.01 0.1 0.2];
w_n = zeros(4,4);
c_n = zeros(4,101);
norm_error = zeros(4,101);
for i = 1:prediction_order
    p_vec(i) = (alpha^(i))/(1-alpha^2);
    if i == 1
        r_vec(i) = 1/(1-alpha^2) + sigma_n_sq;
    else
        r_vec(i) = (alpha^(i-1))/(1-alpha^2);
    end
end

R = toeplitz(r_vec);
w_star = inv(R)*p_vec';
norm_w_star_sq = norm(w_star)^2;


% steepest descent calculation and comparrison vs the ideal estimator
for i = 1:4
    for j =1:101
        c_n(i,j) = norm(w_n(:,i)-w_star)^2;
        norm_error(i,j) = 10*log10(c_n(i,j)/norm_w_star_sq);
        w_n(:,i) = w_n(:,i) + mu(i)*(p_vec'-R*w_n(:,i));
    end
end
iter = 0:1:100;
for k = 1:4
    plot(iter,norm_error(k,:));
    hold on
end
ylim([-100 20]);
legend ("\mu = 0.001", "\mu = 0.01", "\mu = 0.1", "\mu = 0.2");
title ("Normalized weight error")
ylabel("weight error [dB]")
xlabel("iteration")
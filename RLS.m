close all
clear all

% noise process and ideal estimator decleration and calculation
fs = 48000;
alpha = 0.9;
sigma_n_sq = 0.5;
N = sqrt(sigma_n_sq)*randn(fs*10,1);
gauss = randn(fs*10,1);
x = filter(1,[1,-alpha],gauss);
z = x + N;
L = 2;
lambda = 0.99;
delta = [10000000 100000 100 0.0001];
w_n = zeros(L, 4);
error = zeros(480000, 4);
c_n = zeros(480000,4);
w_star = [0.6593; 0.1978];
norm_w_star_sq = norm(w_star)^2; 

% RLS calculation and comparrison vs the ideal estimator
for j = 1:4
    p = 1/delta(j)*eye(L);
    for i = L:(480000-1)
        error(i+1, j) = z(i+1) - w_n(:,j)'*flip(z(i-1:i));
        k = ((1/lambda)*p*flip(z(i-1:i))/(1+(1/lambda)*flip(z(i-1:i))'*p*flip(z(i-1:i))));
        w_n(:,j) = w_n(:,j) + k*error(i+1, j);
        p = (1/lambda)*p - (1/lambda)*k*flip(z(i-1:i))'*p;
        c_n(i+1, j) = (norm(w_n(:,j) - w_star))^2;
    end
end

for i = 1:4
    NR(i) = 10*log10(var(z)/var(error(:,i)));
end


iter = 1:1:480000;
plot(iter, 10*log10(c_n(:,1)/norm_w_star_sq))
hold on
plot(iter, 10*log10(c_n(:,2)/norm_w_star_sq))
plot(iter, 10*log10(c_n(:,3)/norm_w_star_sq))
plot(iter, 10*log10(c_n(:,4)/norm_w_star_sq))
xlabel("iteration")
ylabel("filter coefficient error [dB]")
title("filter coefficient error as a fucntion of iteration L = 2")
legend("\delta = 10000000", "\delta = 100000", "\delta = 100", "\delta = 0.0001")

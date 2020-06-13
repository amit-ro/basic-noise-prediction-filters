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
mu = [0.01; 0.001; 0.0001];
w_n1 = zeros(1,3);
w_n2 = zeros(2,3);
w_n4 = zeros(4,3);
c_n1 = zeros(480000,3);
c_n2 = zeros(480000,3);
c_n4 = zeros(480000,3);
error1 = zeros(480000,3);
error2 = zeros(480000,3);
error4 = zeros(480000,3);

for i = 1:4
    p_vec(i) = (alpha^(i))/(1-alpha^2);
    if i == 1
        r_vec(i) = 1/(1-alpha^2) + sigma_n_sq;
    else
        r_vec(i) = (alpha^(i-1))/(1-alpha^2);
    end
end

R1 = toeplitz(r_vec(1));
R2 = toeplitz(r_vec(1:2));
R4 = toeplitz(r_vec(1:4));

w_star1 = inv(R1)*p_vec(1)';
w_star2 = inv(R2)*p_vec(1:2)';
w_star4 = inv(R4)*p_vec(1:4)';

norm_w_star1_sq = norm(w_star1)^2;
norm_w_star2_sq = norm(w_star2)^2;
norm_w_star4_sq = norm(w_star4)^2;

% LMS calculation and comparrison vs the ideal estimator
for i = 1:3
    error1(1, i) = z(1);
    for j = 1:480000
        if j < 2
            error1(j, i) = z(1);
        else
            error1(j, i) = z(j) - w_n1(:,i)'*z(j-1);
            w_n1(:,i) = w_n1(:,i) + mu(i)*error1(j, i)*z(j-1);
        end
        c_n1(j, i) = (norm(w_n1(:,i) - w_star1))^2;
        if j < 3
           error2(j, i) = z(j);
        else
           error2(j, i) = z(j) - w_n2(:,i)'*flip(z(j-2:j-1));
           w_n2(:,i) = w_n2(:,i) + mu(i)*error2(j, i)*flip(z(j-2:j-1));
        end
        c_n2(j, i) = (norm(w_n2(:,i) - w_star2))^2;
        if j < 5
           error4(j, i) = z(j);
        else
           error4(j, i) = z(j) - w_n4(:,i)'*flip(z(j-4:j-1));
           w_n4(:,i) = w_n4(:,i) + mu(i)*error4(j, i)*flip(z(j-4:j-1));
        end
        c_n4(j, i) = (norm(w_n4(:,i) - w_star4))^2;
    end
end

for i = 1:3
    NR1(i) = 10*log10(var(z)/var(error1(:,i)));
    NR2(i) = 10*log10(var(z)/var(error2(:,i)));
    NR4(i) = 10*log10(var(z)/var(error4(:,i)));
end


iter = 1:1:480000;
figure (1)
plot(iter, 10*log10(c_n1(:,1)/norm_w_star1_sq))
hold on
plot(iter, 10*log10(c_n1(:,2)/norm_w_star1_sq))
plot(iter, 10*log10(c_n1(:,3)/norm_w_star1_sq))
xlabel("iteration")
ylabel("filter coefficient error [dB]")
title("filter coefficient error as a fucntion of iteration L = 1")
legend("\mu = 0.01", "\mu = 0.001", "\mu = 0.0001")

figure (2)
plot(iter, 10*log10(c_n2(:,1)/norm_w_star2_sq))
hold on
plot(iter, 10*log10(c_n2(:,2)/norm_w_star2_sq))
plot(iter, 10*log10(c_n2(:,3)/norm_w_star2_sq))
xlabel("iteration")
ylabel("filter coefficient error [dB]")
title("filter coefficient error as a fucntion of iteration L = 2")
legend("\mu = 0.01", "\mu = 0.001", "\mu = 0.0001")

figure (3)
plot(iter, 10*log10(c_n4(:,1)/norm_w_star4_sq))
hold on
plot(iter, 10*log10(c_n4(:,2)/norm_w_star4_sq))
plot(iter, 10*log10(c_n4(:,3)/norm_w_star4_sq))
xlabel("iteration")
ylabel("filter coefficient error [dB]")
title("filter coefficient error as a fucntion of iteration, L = 4")
legend("\mu = 0.01", "\mu = 0.001", "\mu = 0.0001")
close all
clear all

% noise process definiton
fs = 48000;
alpha = 0.9;
sigma_n_sq = 0.5;
N = sqrt(sigma_n_sq)*randn(480000,1);
gauss = randn(480000,1);
x = filter(1,[1,-alpha],gauss);
z = x + N;

% calculation of the cross correlatiom amd auto correlation
for i = 1:5
    p_vec(i) = (alpha^(i))/(1-alpha^2);
    if i == 1
        r_vec(i) = 1/(1-alpha^2) + sigma_n_sq;
    else
        r_vec(i) = (alpha^(i-1))/(1-alpha^2);
    end
end
R1 = toeplitz(r_vec(1));
R2 = toeplitz(r_vec(1:2));
R3 = toeplitz(r_vec(1:3));
R4 = toeplitz(r_vec(1:4));
R5 = toeplitz(r_vec);

% calculation of the ideal estimator
w1 = inv(R1)*p_vec(1)';
w2 = inv(R2)*p_vec(1:2)';
w3 = inv(R3)*p_vec(1:3)';
w4 = inv(R4)*p_vec(1:4)';
w5 = inv(R5)*p_vec(1:5)';

% ideal estimator is applied and its performance is measured
z_p1 = filter([0; w1], 1, z);
z_p2 = filter([0; w2], 1, z);
z_p3 = filter([0; w3], 1, z);
z_p4 = filter([0; w4], 1, z);
z_p5 = filter([0; w5], 1, z);

e1 = z-z_p1;
e2 = z-z_p2;
e3 = z-z_p3;
e4 = z-z_p4;
e5 = z-z_p5;

NR1 = 10*log10(var(z)/var(e1));
NR2 = 10*log10(var(z)/var(e2));
NR3 = 10*log10(var(z)/var(e3));
NR4 = 10*log10(var(z)/var(e4));
NR5 = 10*log10(var(z)/var(e5));
  
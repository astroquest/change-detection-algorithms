close all
clear all
%% DISTRIBUTIONS

% healthy Gaussian
mu_0 = 0;
sigma = 0.5;

% mean change for faulty system
mu_1 = 0.5;


%% SIMULATION PARAMETERS

N = 1000;        % nr of discrete steps
k0 = N - 0.5*N;  % discrete instance at which fault occurs


%% DETECTOR SETTINGS

%%% detector = 0 (no detection), 
% = 1 (deterministic -> limit check), 
% = 2 (simple prob -> one-sample t-test),
% = 22 (simple prob -> moving window t-test),
% = 3 (adv. prob -> CUSUM)

detector = 22;


%% LIMIT CHECK

x_upper = 1.5;
x_lower = -x_upper;


%% SIMPLE PROBABILISTIC

alpha = 0.05;  % false alarm rate
N0 = 100;      % size of first sliding window
N1 = 100;      % size of second sliding window


%% ADVANCED PROBABILISTIC

h = 15;  % threshold
mu_1_t = 1;


%% SPACE ALLOCATION

x = zeros(N, 1);
mu_est = zeros(N,1);
var_est = zeros(N,1);
sigma_est = zeros(N,1);
t = zeros(N,1);
p_upper = zeros(N,1);
p_lower = zeros(N,1);

mu_0_est = zeros(N0,1);
mu_1_est = zeros(N0,1);
var_0_est = zeros(N1,1);
var_1_est = zeros(N1,1);

s = zeros(N,1);
S = zeros(N,1);
m = zeros(N,1);
g = zeros(N,1);

s_ideal = zeros(N,1);
S_ideal = zeros(N,1);
m_ideal = zeros(N,1);
g_ideal = zeros(N,1);


%% SIMULATION

rng('default')
rng(2)

for k=1:N
    % generate Gaussian samples
    if k < k0
        x(k) = normrnd(mu_0, sigma);
    elseif k >= k0
        x(k) = normrnd(mu_1, sigma);
    end

    if detector == 2
        % adjust estimation of mean and variance at every step
        if k == 1
            mu_est(k) = x(k);
            var_est(k) = x(k)^2;
        elseif k > 1
            mu_est(k) = mu_est(k-1) + 1/k*( x(k) - mu_est(k-1) );
            var_est(k) = ( k - 2 )/( k - 1 )*var_est(k-1) + 1/k*( x(k) - mu_est(k-1) )^2;
        end
        sigma_est(k) = sqrt(var_est(k));

        if k < k0
            t(k) = ( mu_est(k) - mu_0 )/( sigma_est(k) / sqrt(N) );
        elseif k >= k0
            t(k) = ( mu_est(k) - mu_1 )/( sigma_est(k) / sqrt(N) );
        end
     
        p_lower(k) = tcdf(t(k), k-1);
        p_upper(k) = tcdf(t(k), k-1, "upper");
        % if minimum p-value <= alpha/2 -> change detected
    end

    if detector == 22
        % start sliding window when enough samples are available
        if k >= N0 + N1
    
            for n0=1:N0
                m = k - N0 - N1 + n0;
                if n0 == 1
                    mu_0_est(n0) = x(m);
                    var_0_est(n0) = x(m)^2;
                elseif n0 > 1       
                    mu_0_est(n0) = mu_0_est(n0-1) + 1/n0*( x(m) - mu_0_est(n0-1) );
                    var_0_est(n0) = ( n0 - 2 )/( n0 - 1 )*var_0_est(n0-1) + 1/n0*( x(n0) - mu_0_est(n0-1) )^2;
                end
            end
    
            for n1=1:N1
                p = k - N1 + n1;
                if n1 == 1
                    mu_1_est(n1) = x(p);
                    var_1_est(n1) = x(p)^2;
                elseif n1 > 1
                    mu_1_est(n1) = mu_1_est(n1-1) + 1/n1*( x(p) - mu_1_est(n1-1) );
                    var_1_est(n1) = ( n1 - 2 )/( n1 - 1 )*var_1_est(n1-1) + 1/n1*( x(n1) - mu_1_est(n1-1) )^2;
                end
            end
    
            t(k) = ( mu_0_est(n0) - mu_1_est(n1) )/sqrt( (N0 - 1)*var_0_est(n0) + (N1 - 1)*var_1_est(n1) ) * sqrt( ( N0*N1*( N0 + N1 - 2 ) )/( N0 + N1 ) );
            p_lower(k) = tcdf(t(k), N0+N1-1);
            p_upper(k) = tcdf(t(k), N0+N1-1, "upper");

        end

    end

    if detector == 3
        % assume mu_0, sigma are known from healthy Gaussian distribution
        s(k) = ( mu_1_t - mu_0 )/sigma^2 * ( x(k) - ( mu_0 + mu_1_t )/2 );
        s_ideal(k) = ( mu_1 - mu_0 )/sigma^2 * ( x(k) - ( mu_0 + mu_1 )/2 );

        if k == 1
            S(k) = s(k);
            S_ideal(k) = s_ideal(k);
        elseif k > 1
            S(k) = S(k-1) + s(k);
            S_ideal(k) = S_ideal(k-1) + s_ideal(k);
        end
        m(k) = min(S);
        m_ideal(k) = min(S_ideal);
        g(k) = S(k) - m(k);
        g_ideal(k) = S_ideal(k) - m_ideal(k);
    end
end


%% PLOTS

K = 1:N;

lw = 1.2;

figure
hold on
plot(K, x, "LineWidth", lw)
title("Samples from changing Gaussian distribution k0=50")
if detector == 1
    yline(x_upper, 'k', "LineWidth", lw-0.2)
    yline(x_lower, 'k', "LineWidth", lw-0.2)
    title("Samples from changing Gaussian distribution with limit check")
    legend("samples", "upper limit", "lower limit", "location", "northwest")
end
grid

if detector == 2
    plot(K, mu_est, "LineWidth", lw)
    plot(K, sigma_est, "LineWidth", lw)
    legend("samples", "\mu est", "\sigma est", "location", "northwest")

    figure
    hold on
    plot(K, p_upper, "LineWidth", lw)
    plot(K, p_lower, "LineWidth", lw)
    yline(alpha/2, 'k', "LineWidth", lw-0.2)
    legend("upper p-value", "lower p-value", "\alpha/2")
    title("p-value fault detection")
    grid
end

if detector == 22
    figure
    hold on
    plot(K, p_upper, "LineWidth", lw)
    plot(K, p_lower, "LineWidth", lw)
    yline(alpha/2, 'k', "LineWidth", lw-0.2)
    legend("upper p-value", "lower p-value", "\alpha/2")
    title("p-value fault detection")
    grid
end

if detector == 3
    figure
    hold on
    plot(K, g_ideal, "LineWidth", lw)
    plot(K, g, "LineWidth", lw)
    yline(h, 'k', "LineWidth", lw-0.2)
    legend("actual \mu_1", "tuned \mu_1", "fault threshold", "location", "northwest")
    title("CUSUM fault detection")
    grid
end

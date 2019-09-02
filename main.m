%% PRM Matlab Workshop 2: Discrete Channel Simulation & Repetition Coding
%% 1c) Simulation of one channel use
clc, clear
% Probability that a bit will be swapped given that it
% is 1 is 0.05
% Probability that a bit will be erased due to buffer
% overflow given that it is 1 is 0.2
% Same probability applies for the case that the bit is 0

xProb = rand;

if (xProb <= 0.5)
     
    x = 0
    yProb = rand;
    
    if (yProb <= 0.75)
        y = '0'
    elseif (yProb > 0.75)&&(yProb <= 0.8)
        y = '1'
    else
        y = 'E'
    end
    
else
    
    x = 1
    yProb = rand;
    
    if (yProb <= 0.75)
        y = '1'
    elseif (yProb > 0.75)&&(yProb <= 0.8)
        y = '0'
    else
        y = 'E'
    end
    
end
%% 1d) Distribution of Swaps & its Poisson Approximation
clc, clear

% Length of bit string
n = 99;
% Probability of swap
p = 0.05;
alpha = p*n;

% All number of swaps that can occur
s = [1:99]';

% Distribution of Swaps (Binomial)
for k = 1:n
    P_S(k, 1) = factorial(n)/(factorial(k)*factorial(n-k)) * p^k * (1 - p)^(n-k);
end
% for k = 1:n
%     p_swap(k, 1) = nchoosek(n, k) * p^k * (1 -
%     p)^(n-k);
% end

figure(1)
subplot(2, 1, 1)
bar(s, P_S)
grid on
axis([0 n 0 0.2])
xlabel('Number of bit swaps')
ylabel('Probability of bit swaps')
title('PMF of Bit Swaps (Binomial)')

% Distribution of Swaps (Poisson Approximation)
P_S_poisson = (alpha.^s)./factorial(s) * exp(-alpha);

figure(1)
subplot(2, 1, 2)
bar(s, P_S_poisson)
grid on
axis([0 n 0 0.2])
xlabel('Number of bit swaps')
ylabel('Probability of bit swaps')
title('PMF of Bit Swaps (Poisson Approximation)')
%% 1e) Simulation of bit swapped with 99 trials
clc, clear

% Probability of bit swap
p2 = 0.05;
% Length of bit string
n = 99;

% Input bit string
X = zeros(n, 1);
% Output bit string
Y = zeros(n, 1);

for k = 1:n
    [X(k), Y(k)] = disChn(p2);
end

% Compare input X and output Y vectors to get a vector
% s that denotes bits that were swapped
S = (X ~= Y);
% Number of swaps (bit errors)
S_num = sum(S)
% Probability of swap
P_s = mean(S)
%% 1e) Empirical PMF Plot of bit swap with 10e3 trials
clc, clear

% Probability of bit swap
p2 = 0.05;
% Length of bit string
n = 99;
% Input bit string
X = zeros(n, 1);
% Output bit string
Y = zeros(n, 1);

% Number of trials
M = 10e3;
% Number of swaps (bit errors) from bit string of
% length n
S_num = zeros(n, 1);
% Probability of bit swap for each number of bit swap
P_s = zeros(n, 1);

for m = 1:M
    for k = 1:n
        [X(k), Y(k)] = disChn(p2);
    end
    
    % Compare input X and output Y vectors to get a vector
    % S that denotes bits that were swapped
    S = (X ~= Y);
    % Number of bits that was swapped
    S_num(m) = sum(S);
end

% Find the probability of each number of swaps
for k = 1:n
    P_s(k) = numel(find(S_num == k))/M;
end

figure(2)
subplot(2, 1, 1)
bar([1:n], P_s)
grid on
axis([0 n 0 0.2])
xlabel('Number of bit swaps')
ylabel('Probability of bit swaps')
title('Empirical PMF of Bit Swaps')

% Exact PMF of Swaps
alpha = p2*n;

% All number of swaps that can occur
s = [1:99]';

% Distribution of Swaps (Binomial)
for k = 1:n
    P_S(k, 1) = factorial(n)/(factorial(k)*factorial(n-k)) * p2^k * (1 - p2)^(n-k);
end

subplot(2, 1, 2)
bar(s, P_S)
grid on
axis([0 n 0 0.2])
xlabel('Number of bit swaps')
ylabel('Probability of bit swaps')
title('Exact PMF of Bit Swaps')
%% 1f) Expected values with different number of trials
clc, clear

% Length of bit string
n = 99;
% Input bit string
X = zeros(n, 1);
% Output bit string
Y = zeros(n, 1);

% Probability of bit swap
p2 = 0.05:0.05:0.8;
% Aveage number (expected value) of swaps from M trials
S_avg = zeros(length(p2), 1);
% Probability of bit swap for each number of bit swap
P_s = zeros(n, 1);

k0 = 1;
for M = [2 5 10 100 1000 10e3]
    
    k1 = 1;
    for p = 0.05:0.05:0.8
        
        for m = 1:M
            
            for k = 1:n
                [X(k), Y(k)] = disChn(p);
            end
            
            % Compare input X and output Y vectors to get a vector
            % S that denotes bits that were swapped
            S = (X ~= Y);
            % Number of bits that was swapped
            S_num(m) = sum(S);
            % s_num vector contains number of swaps collected from
            % 10e3 trials
        end
        S_avg(k1,1) = mean(S_num);
        
        k1 = k1 + 1;
    end
    
    figure(3)
    subplot(3, 2, k0)
    bar(p2, S_avg)
    grid on
    axis([0 0.9 0 n])
    xlabel('Probability of bit swaps')
    ylabel('Average number of bit swaps')
    title({['Number of trials M = ', num2str(M)], 'Bit string length n = 99'})
    
    k0 = k0 + 1;
end
%% 1h) Average number of decoder error for bit string length n = 99
clc, clear

% Length of bit string
n = 99;
% Input bit string
X = zeros(n, 1);
% Output bit string
Y = zeros(n, 1);

% Probability of bit swap
p2 = 0.05:0.05:0.45;
% Number of decoder error
error = 0;
errorAvg99 = zeros(length(p2), 1);

% Number of trials
M = 10e3;

k1 = 1;
for p = 0.05:0.05:0.45
    for m = 1:M
        for k = 1:n
            [X(k), Y(k)] = disChn(p);
        end
        
        S = (X ~= Y);
        if (sum(S) >= n/2)
            error = error + 1;
        end
    end
    errorAvg99(k1) = error/M;
    error = 0;
    k1 = k1 + 1;
end

figure(4)
bar(p2, errorAvg99)
ylim([0, 0.2])
grid on
legend('Bit String Length n = 99')
xlabel('Probability of bit swaps')
ylabel('Number of decoder errors')
title('Number of trials M = 10000')
%% 1i) Average number of decoder error for bit string length n = 5
clc

% Length of bit string
n = 5;
% Input bit string
X = zeros(n, 1);
% Output bit string
Y = zeros(n, 1);

% Probability of bit swap
p2 = 0.05:0.05:0.45;
% Number of decoder error
error = 0;
errorAvg5 = zeros(length(p2), 1);

% Number of trials
M = 10e3;

k1 = 1;
for p = 0.05:0.05:0.45
    for m = 1:M
        for k = 1:n
            [X(k), Y(k)] = disChn(p);
        end
        
        S = (X ~= Y);
        if (sum(S) >= n/2)
            % If number of bit swaps is more than half the number
            % of elements in Y, an error has occurred
            error = error + 1;
        end
    end
    errorAvg5(k1) = error/M;
    error = 0;
    k1 = k1 + 1;
end

figure(5)
bar(p2, [errorAvg99 ,errorAvg5])
ylim([0, 0.5])
grid on
legend({'Bit String Length n = 99', 'Bit String Length n = 5'})
xlabel('Probability of bit swaps')
ylabel('Average number of errors')
title('Number of trials M = 10000')
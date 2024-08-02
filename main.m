clear all;
%
L = 1e6;
%
KdB = 5; % dB
K = db2pow(KdB);
%
d = sqrt(K/(K+1));
sigma = 1/sqrt(K+1);
G = (randn(L, 1) + 1i*randn(L, 1))/sqrt(2);
phi_0 = 2*pi*rand(L, 1);

% Shadowing Characterization
k = 2; % half d.o.f.
lambda = 1.5;
Omega = 1/(k+lambda); %
xi = sqrt(Omega/2 * ncx2rnd(2*k, 2*lambda, [L, 1]));

% Received Signal
S = d*xi.*exp(1i*phi_0) + sigma*G;

% Signal-to-Noise Ratio
avgsnr = 1;
snr = avgsnr * abs(S).^2;

% =========================================================================
% Eq. (12): Raw Moments
nn = 1:7;
for n = nn
    rawMoment_sim(n) = mean(snr.^n);
    
    rawMoment_ana(n) = 0;
    for ii = 0:n
        rawMoment_ana(n) = rawMoment_ana(n) + factorial(n).*(sigma^2*avgsnr).^n...
            .* nchoosek(n,ii).*(Omega*d^2/sigma^2).^ii...
            .* laguerreL(ii, k-1, -lambda);
    end
end

figure;
semilogy(nn, rawMoment_sim, '-xk', 'LineWidth', 1.5); hold on;
semilogy(nn, rawMoment_ana, '--or', 'LineWidth', 1.5); hold on;
legend('Simulation', 'Analytical Eq. (12)', 'Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$\mathbf{E}\{ \gamma^n \}$', 'Interpreter', 'latex');

% =========================================================================
% Eq. (14): Probability Density Function
MGF_snr = @(s) (1-sigma^2*avgsnr*s).^(k-1)./(1-(sigma^2+Omega*d^2)*avgsnr*s).^k...
    .* exp(lambda*Omega*d^2*avgsnr*s./(1-(sigma^2+Omega*d^2)*avgsnr*s));
L = max(-1./[sigma^2*avgsnr, (sigma^2+Omega*d^2)*avgsnr])+1/2;
T = 1e3;
f_snr = @(x) integral(@(s) exp(x*s).*MGF_snr(-s), L-1i*T, L+1i*T, 'ArrayValued', true)/(1i*2*pi);

figure;
[y, x] = ksdensity(snr, 'Numpoints', 1e2); x = x(x>0);
histogram(snr, 'NumBins', 100, 'Normalization', 'pdf'); hold on;
plot(x, f_snr(x), '--r', 'LineWidth', 1.5);

% =========================================================================
% Eq. (15): Probability Density Function
F = @(s, x) MGF_snr(-s+avgsnr^(-1)/(sigma^2+Omega*d^2)) .* exp(-x*avgsnr^(-1)/(sigma^2+Omega*d^2));
L = 1/2;
T = 1e3;
f_snr = @(x) integral(@(s) exp(x*s).*F(s, x), L-1i*T, L+1i*T, 'ArrayValued', true)/(1i*2*pi);

figure;
[y, x] = ksdensity(snr, 'Numpoints', 1e2); x = x(x>0);
histogram(snr, 'NumBins', 100, 'Normalization', 'pdf'); hold on;
plot(x, f_snr(x), '--r', 'LineWidth', 1.5);

% =========================================================================
% Eq. (13): Probability Density Function
f_snr = @(x) 1/(sigma^2*avgsnr)/(1+Omega*d^2/sigma^2)^k...
    .* exp(-x/(sigma^2*avgsnr)-lambda)...
    .* HypergeomePsi2(k, lambda/(1+Omega*d^2/sigma^2),...
                         Omega*d^2/(sigma^2+Omega*d^2)/(sigma^2*avgsnr), x);

figure;
[y, x] = ksdensity(snr, 'Numpoints', 1e2); x = x(x>0);
histogram(snr, 'NumBins', 100, 'Normalization', 'pdf'); hold on;
plot(x, f_snr(x), '--r', 'LineWidth', 1.5);

% =========================================================================
% Eq. (22): Probability Density Function
if (floor(k) == k) && (k>=1)
    f_snr = @(x) exp(-lambda*Omega*d^2/(sigma^2+Omega*d^2))/(avgsnr*sigma^2)/(1+Omega*d^2/sigma^2)^k...
        * exp(-x/avgsnr/(sigma^2+Omega*d^2));
    
    sum_j = @(x) 0;
    for j = 0:k-1
        sum_j = @(x) sum_j(x) + (-1)^j*pochhammer(1-k,j)/factorial(j)...
            .* (Omega*d^2/sigma^2/lambda*x/(avgsnr*sigma^2)).^(j/2)...
            .* besseli(j, 2/(sigma^2+Omega*d^2)*sqrt(lambda*Omega*d^2*x/avgsnr));
    end

    f_snr = @(x) f_snr(x) .* sum_j(x);
end

figure;
[y, x] = ksdensity(snr, 'Numpoints', 1e2); x = x(x>0);
histogram(snr, 'NumBins', 100, 'Normalization', 'pdf'); hold on;
plot(x, f_snr(x), '--r', 'LineWidth', 1.5);

% =========================================================================
% Eq. (18): Cumulative Distribution Function
A = (1+Omega*d^2/sigma^2)^(-1);
B = Omega*d^2/(sigma^2+Omega*d^2);

F_snr = @(x) x/(sigma^2*avgsnr)*A^k*exp(-B*lambda)...
    .* exp(-A/(sigma^2*avgsnr)*x)...
    .* HypergeomePhi3n([1-k, 1],2,[-B, A, A*B*lambda]/(sigma^2*avgsnr), x);

figure;
[y, x] = ecdf(snr); 
plot(x, y, 'LineWidth', 1.5); hold on;
xx = linspace(min(x), max(x), 1e3);
plot(xx, F_snr(xx), '--r', 'LineWidth', 1.5);
axis([0 12 min(y) max(y)])

% =========================================================================
% Eq. (24): Cumulative Distribution Function - Kappa-Mu equivalent
F_kmu = @(g, kappa, mu, x) 1 - marcumq(sqrt(2*kappa*mu), sqrt(2*(1+kappa)*mu*x/g), mu);

if (floor(k) == k) && (k>=1)
    
    F_snr = @(x) 0;
    for j = 0:(k-1)
        C_j = (-1)^j/factorial(j)*(Omega*d^2/sigma^2)^j/(Omega*d^2/sigma^2+1)^(k-1)*pochhammer(1-k, j);
        g = (j+1+lambda*Omega*d^2/(sigma^2+Omega*d^2))*(sigma^2+Omega*d^2)*avgsnr;
        mu = j+1;
        kappa = lambda*Omega*d^2/(sigma^2+Omega*d^2)/(j+1);
        %
        F_snr = @(x) F_snr(x) + C_j*F_kmu(g, kappa, mu, x);
        %
    end
end

figure;
[y, x] = ecdf(snr); 
plot(x, y, 'LineWidth', 1.5); hold on;
xx = linspace(min(x), max(x), 1e3);
plot(xx, F_snr(xx), '--r', 'LineWidth', 1.5);
axis([0 12 min(y) max(y)])

% =========================================================================
% Outage Probability

F_snr_asymp = @(x) x/(sigma^2*avgsnr)*A^k*exp(-B*lambda);

xdB = -5:2.5:40;
for ix = 1:length(xdB)
    %
    x = db2pow(xdB(ix));
    %
    Pout_sim(ix) = mean(x*snr < 1);
    Pout_ana(ix) = F_snr(x^(-1));
    Pout_asymp(ix) = F_snr_asymp(x^(-1));
end
%
figure;
semilogy(xdB, Pout_sim, 'ok'); hold on;
semilogy(xdB, Pout_ana, '-k', 'LineWidth', 1.5); hold on;
semilogy(xdB, Pout_asymp, '--r', 'LineWidth', 1.5); hold on;
grid on;
xlabel('$\bar{\gamma}/\gamma_{\rm th}$ [dB]', 'Interpreter', 'Latex');
ylabel('$P_{\rm out}$', 'Interpreter', 'Latex');
legend('Theoretical', 'Simulated', 'Asymptotic', 'Interpreter', 'Latex');
%
clear all;
warning off;
%
L = 5e6;
%
KdB = 13; % dB
K = db2pow(KdB);
%
d = sqrt(K/(K+1));
sigma = 1/sqrt(K+1);
G = (randn(L, 1) + 1i*randn(L, 1))/sqrt(2);
phi_0 = 2*pi*rand(L, 1);

% Shadowing Characterization
k = 2.5; % half d.o.f.
lambda = 0;
Omega = 1/(k+lambda); %
xi = sqrt(Omega/2 * ncx2rnd(2*k, 2*lambda, [L, 1]));

% Received Signal
S = d*xi.*exp(1i*phi_0) + sigma*G;

% Signal-to-Noise Ratio
avgsnr = 1;
snr = avgsnr * abs(S).^2;

% =========================================================================
% Eq. (18): Prony Method

x0_ = [  0, 1e2, 1e3, 1e4, 1e5];
xL_ = [1e2, 1e3, 1e4, 1e5, 1e6];
M_  = [  4,   4,    4,  4,   4];
L_ = 1e2*ones(size(x0_));
S = length(x0_);

f = @(x) log(1+x);
% MGF = @(s) (1-sigma^2*s).^(k-1)./(1-(sigma^2+Omega*d^2)*s).^k...
%     .* exp(lambda*Omega*d^2*s./(1-(sigma^2+Omega*d^2)*s));


MGFg = @(s) (1-sigma^2*s).^(k-1)./(1-(sigma^2+Omega*d^2)*s).^k;

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

f_approx = @(x) 0;
for s = 1:S
    x0 = x0_(s);
    xL = xL_(s);
    L = L_(s);
    M = M_(s);

    [Tj, cj] = PronyMethod(M, L, f, x0, xL);
    Tj_{s} = Tj;
    cj_{s} = cj;

    for j=1:M
        f_approx = @(x) f_approx(x) + cj(j)*exp(-Tj(j)*x).*(x<xL).*(x>=x0);
    end
end
% 
% xx = linspace(min(x0_), max(xL_), 1e3);
% figure;
% plot(xx, f(xx), 'LineWidth', 1.5); hold on;
% plot(xx, f_approx(xx), ':b', 'LineWidth', 1.5,...
%     'MarkerIndices', 1:1e2:L+1); hold on;
% xlabel('$x$', 'Interpreter', 'LaTex');
% ylabel('$f(x)$', 'Interpreter', 'LaTex');
% legend('Exact', 'Approximation');
% set(gca, 'XScale', 'Log')

MGFg = @(s) (1-sigma^2*s).^(k-1)./(1-(sigma^2+Omega*d^2)*s).^k...
    .* exp(lambda*Omega*d^2*s./(1-(sigma^2+Omega*d^2)*s));

ErgodicCap = @(x) 0;
for s = 1:S
    Tj = Tj_{s};
    cj = cj_{s};
    x0 = x0_(s);
    xL = xL_(s);
    M = M_(s);
    
    for j=1:M
        C = 1/2;
        T = 1e3;       

        if x0>0
            Fp = @(p,x,z1,z2) (MGFg(-Tj(j)*x-x*p/z2)-MGFg(-Tj(j)*x-x*p/z1))./p.*exp(p);
            G = @(x,z1,z2) real(integral(@(p) Fp(p,x,z1,z2), C-1i*T, C+1i*T,...
                'AbsTol', 1e-20, 'RelTol', 1e-20)/(1i*2*pi));

            ErgodicCap = @(x) ErgodicCap(x) + cj(j)*(G(x,x0,xL));
        else
            Fp = @(p,x,z) MGFg(-Tj(j)*x-x*p/z)./p.*exp(p);
            G = @(x,z) real(integral(@(p) Fp(p,x,z), C-1i*T, C+1i*T,...
                'AbsTol', 1e-20, 'RelTol', 1e-20)/(1i*2*pi));

            ErgodicCap = @(x) ErgodicCap(x) + cj(j)*G(x,xL);
        end
    end
end

ErgodicCap_Inf = @(x) 0;
for s = 1:S
    Tj = Tj_{s};
    cj = cj_{s};
    x0 = x0_(s);
    xL = xL_(s);
    M = M_(s);
    
    for j=1:M
        % ErgodicCap = @(x) ErgodicCap(x)...
        %     + cj(j)*(exp(-Tj(j)*x*snr).*(snr<xL/x)...
        %             -exp(-Tj(j)*x*snr).*(snr<x0/x));
        C = 1/2;
        T = 1e3;       

        if x0>0
            Fp = @(p,x,z1,z2) (MGFg(-Tj(j)*x-x*p/z2)-MGFg(-Tj(j)*x-x*p/z1))./p.*exp(p);
            G = @(x,z1,z2) real(integral(@(p) Fp(p,x,z1,z2), C-1i*T, C+1i*T,...
                'AbsTol', 1e-20, 'RelTol', 1e-20)/(1i*2*pi));

            ErgodicCap_Inf = @(x) ErgodicCap_Inf(x) + cj(j)*exp(Tj(j))*G(x,x0,xL);
        else
            Fp = @(p,x,z) MGFg(-Tj(j)*x-x*p/z)./p.*exp(p);
            G = @(x,z) real(integral(@(p) Fp(p,x,z), C-1i*T, C+1i*T,...
                'AbsTol', 1e-20, 'RelTol', 1e-20)/(1i*2*pi));

            ErgodicCap_Inf = @(x) ErgodicCap_Inf(x) + cj(j)*exp(Tj(j))*G(x,xL);
        end
    end
end
% =========================================================================
% Ergodic Capacity
A = (1+Omega*d^2/sigma^2)^(-1);
B = Omega*d^2/(sigma^2+Omega*d^2);

F_snr_asymp = @(x) 1-1/(sigma^2*x)*A^k*exp(-B*lambda);

xdBsim = -10:5:40;
for ix = 1:length(xdBsim)
    %
    x = db2pow(xdBsim(ix));
    %
    EC_sim(ix) = mean(log(1+x*snr))/log(2);
end

xdBana = -10:2:40;
for ix = 1:length(xdBana)
    %
    x = db2pow(xdBana(ix));
    %
    % EC_ana(ix) = mean(f_approx(x*snr))/log(2);
    EC_ana(ix) = ErgodicCap(x)/log(2);
    EC_asymp(ix) = ErgodicCap_Inf(x)/log(2);
    EC_AWGN(ix) = log(1+x)/log(2);
end
%
figure(1);
plot(xdBsim, EC_sim, 'ok'); hold on;
plot(xdBana, EC_ana, '-k', 'LineWidth', 1.5); hold on;
plot(xdBana, EC_asymp, '--r', 'LineWidth', 1.5); hold on;
plot(xdBana, EC_AWGN, '--b', 'LineWidth', 1.5); hold on;

filename = sprintf('EC_sim_k%d.mat', 100*k);
save(filename, 'EC_sim');
filename = sprintf('EC_ana_k%d.mat', 100*k);
save(filename, 'EC_ana');
filename = sprintf('EC_asymp_k%d.mat', 100*k);
save(filename, 'EC_asymp');
save('EC_AWGN.mat', 'EC_AWGN');

xlabel('$\bar{\gamma}$ [dB]', 'Interpreter', 'Latex');
ylabel('$\bar{C}$ [bps/Hz]', 'Interpreter', 'Latex');

legend('Simulated', 'Approximated', 'Asymptotic', 'Interpreter', 'Latex');
set(gca, 'FontSize', 14);
% axis([min(xdBana) max(xdBana) 0 12])
grid on;
%
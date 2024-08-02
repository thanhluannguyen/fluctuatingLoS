function [Tj, cj] = PronyMethod(M, L, f, x0, xL)

    xx = linspace(x0, xL, L+1); % x0 = 1;
    h = xx(2)-xx(1);
    
    fvec = f(xx).';
    
    % Step 0: Construct Hf
    Hf = ones(L-M+1, M+1);
    for j = 1:(L-M+1)
        Hf(j, :) = fvec(j+(0:M)).';
    end
    
    % Step 1: Compute the Right Singular Vector - Do not use null(Hf)
    % p = fsolve(@(p) rightSingular(Hf, p), zeros(M, 1));
    % p = [p; 1];

    p = pinv(Hf(:,1:M))*Hf(:,end);
    p = [-p; 1];
    
    % Step 2: Compute the zeros zj, j = 1:M
    poly = poly2sym(flipud(p));
    r = solve(poly==0);

    zj = zeros(M, 1);
    for k = 1:M
        zj(k) = double(vpa(r(k)));
    end
    
    % Step 3: Compute the coefficients dj using Least-Square Method
    % Z = ones(L+1, M); 
    % for k = 2:(L+1)
    %     Z(k, :) = zj.^k;
    % end
    % dj = pinv(Z)*fvec;

    options = optimoptions('fsolve','Algorithm', 'levenberg-marquardt');
    x = fsolve(@(x) solveProny(x, zj, fvec), zeros(1, M), options);
    dj = x.';
    
    % Step 4: Compute Tj, cj
    Tj = -log(zj)/h;
    cj = dj.*exp(Tj*x0);
    
function F = solveProny(x, zj, fvec)
    M = length(zj);
    L = length(fvec)-1;
    %
    for k=0:L
        F(k+1) = 0;
        for j=1:M
            F(k+1) = F(k+1) + x(j)*zj(j)^k;
        end
        F(k+1) = F(k+1)-fvec(k+1);
    end
    %    
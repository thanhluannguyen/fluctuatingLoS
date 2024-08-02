function f = HypergeomePhi3n(B, c, D, t)
    n = length(D);
    %
    F = @(s) s.^(-c).*exp(D(n)./s);
    for i = 1:n-1
        F = @(s) F(s).*(1-D(i)./s).^(-B(i));
    end
    % Choose L slightly larger than the real part of all the singularities
    L = max([D(1:n-1), 0])+1/2;
    %
    T = 1e4;
    invLaplace = integral(@(s) exp(s*t).*F(s), L-1i*T, L+1i*T, 'ArrayValued', true,...
        'RelTol', 1e-6, 'AbsTol', 1e-10)/(2*pi*1i);
    %
    f = gamma(c)*t.^(1-c).*real(invLaplace);
end
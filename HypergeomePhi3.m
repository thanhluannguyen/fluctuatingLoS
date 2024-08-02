%
%   \Phi_3(b, c; x*t, y*t)
%
function f = HypergeomePhi3(b, c, x, y, t)
    F = @(s) s.^(b-c).*(s-x).^(-b).*exp(y./s);
    % Choose L slightly larger than the real part of all the singularities
    L = max(x, 0)+1/2;
    %
    T = 1e3;
    invLaplace = integral(@(s) exp(s*t).*F(s), L-1i*T, L+1i*T, 'ArrayValued', true,...
        'RelTol', 1e-6, 'AbsTol', 1e-10)/(2*pi*1i);
    %
    f = gamma(c)*t.^(1-c).*real(invLaplace);
end
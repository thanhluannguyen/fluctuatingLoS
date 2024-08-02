function f = HypergeomePsi2(a, w, z, t)
    %
    f = exp(z*t+w).*HypergeomePhi3n(1-a, 1, [-z, w*z], t);
    % f = exp(z*t+w).*HypergeomePhi3(1-a, 1, -z, w*z, t);
end
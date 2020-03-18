%% Setups
n = 400;
beta1 = 1;
beta2 = -0.5;
x1 = normrnd(0, 1, [400, 1]);
x2 = chi2rnd(1, [400, 1]);
u1 = evrnd(1, 0, [400, 1]);
u2 = evrnd(1, 0, [400, 1]);

%% (1) data generating

y = zeros(400,1);

for i= 1:400
    if x1(i)*beta1 + u1(i) > x2(i)*beta2 + u2(i)
        y(i) = 1;
    end
end


%% (2) log-likelihood function
F = inline('exp(-exp(x2*beta2 - x1*beta1))', 'x1', 'x2', 'beta1', 'beta2');



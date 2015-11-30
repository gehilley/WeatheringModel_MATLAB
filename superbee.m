function limiter = superbee(u)

n = length(u(1,:));
nx = length(u(:,1));

zer = zeros(nx,n);
one = ones(nx,n);
twos = 2.*ones(nx,n);

% Precondition u for nans and infs:

i = find(isnan(u));
u(i) = -2;

i = find(isinf(u) & (u > 0));
u(i) = 2;
i = find(isinf(u) & (u < 0));
u(i) = -2;

term1 = zer;
term2 = (2.*u > one).*one + (2.*u <= one).*u;
term3 = (u > twos).*twos + (u <= twos).*u;

limiter1 = (term1 > term2).*term1 + (term1 <= term2).*term2;
limiter = (limiter1 > term3).*limiter1 + (limiter1 <= term3).*term3;


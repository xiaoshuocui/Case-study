function [u,k] = conjugate(A,b,u0) %be careful about the dimension of A b and u0
p = b - A*u0;
k = 0;
r = b - A*u0;
u = u0;
while norm(r) > 1e-8
    m = A*p;
    a = (norm(r))^2/(p'*m);
    u = u + a*p;
    rpre = r;
    r = r - a*m;
    beta = (norm(r))^2/(norm(rpre))^2;
    p = r + beta*p;
    k = k+1;
end
disp(u)
disp(k)
end

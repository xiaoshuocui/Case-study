function [u,iter,soliter] = conjugate(A,b,u0,iter)
%Function to perform the conjugate gradient method. This method reuquires
%the matrix A to be real  and symmetric positive definite. f is the RHS of
%the equation Ax = b and u0 is the first guess at the solution.
r = b - A*u0; % This is the first rsidual, r0 in the notes.
p = r; % As part of the method we choose p0 = r0.
mag = norm(r); % This is the 2-norm/magnitude of r which we will also need later in the method. 
%This part is just initialising the variables we will require later.
u = u0;
soliter = [];
soliter(:,1) = u0;
while mag > 1e-8 
    N = A*p;% This is defined now, as we will require this product several times during the method.
    %alpha = mag/(p'*N); 
    alpha = p'*r/(p'*N);
    u = u + alpha*p;
    soliter(:,iter+2) = u;
    r = b - A*u;
    %Update the sol and residual (u and r).
    beta = -(p'*A*r)/(p'*N);
    %beta = norm(r)/mag; 
    p = r + beta*p;
    mag = norm(r);
    iter = iter + 1;
    
end
end


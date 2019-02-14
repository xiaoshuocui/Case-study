function [u,max_error,U] = multi(k,f,g,max_depth,smooth_it)
%Function designed to use a 2-grid iteration to solve the system on a
%coarser grid. We can use linear interpolation to restrict/prolong the
%resulting solution u.
shift = 2^k;
A = matrix(shift,shift); % Creates the (2^n - 1) square finite difference matrix which we can use with our multigrid method.

dx = 1/shift;
dy = 1/shift;
x = dx:dx:1-dx;
y = dy:dy:1-dy;
x1 = 0:dx:1;%Set up the mesh for x and y and pull together to solve the system.
y1 = 0:dy:1;%This is full mesh including boundary.
[X,Y] = meshgrid(x,y);
[X1,Y1] = meshgrid(x1,y1);
F = f(X,Y);
b = reshape(F,(shift-1)*(shift-1),1);


u = multigrid(A,b,max_depth,smooth_it);


n = shift; m = shift;
U = reshape(u,m-1,n-1);%matrix of approximate solution.
U1 = [zeros(1,n+1); [zeros(m-1,1), U, zeros(m-1,1)]; zeros(1,n+1)]; %Add back in boundary conditions.

u_exact = reshape(g(X,Y),(n-1)*(m-1),1);
U_exact = g(X1,Y1);
error1 = u - u_exact;
%Examine the errors, and check the maximum to observe how the errors change when we alter the mesh spacings in both the x and y directions.
max_error = norm(error1,Inf);


figure;
subplot(1,2,1)
title('approximate solution')
mesh(X1,Y1,U1)

subplot(1,2,2)
title('exact solution')
mesh(X1,Y1,U_exact)


end


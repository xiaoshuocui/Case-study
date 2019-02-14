function [x] = multigrid(A,b,max_depth,smooth_it)
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here
l = length(A);
if l == 1 || max_depth == 0
    x = A\b;
else
    shift = sqrt(length(A)) + 1; shift1 = shift/2;
    P1 = zeros(shift - 1,shift1 - 1);% making the 1D prolongation matrix.
    P1(1:shift + 1:end) = 1;
    P1(3:shift + 1:end) = 1;
    P1(2:shift + 1:end) = 2;
    P1 = 1/2*sparse(P1);
    
    R1 = 1/2*P1'; %Define the 1D restriction operator.
    
    R = kron(R1,R1); %Create the 2D restriction using the 'kron'.
    P = 4*R';
    Abar = R*A*P; %transform t
    h = @(l) 1e-8*(1000/24*(log(l)/log(2) + 1)- 976/24);
    smooth_method = @(A,b,count,x) smoother(tril(A),-triu(A,1),b,count,x);
    x = zeros(length(A),1);
    while norm(A*x - b) > h(l)
        x_s = smooth_method(A,b,smooth_it,x);
        
        r = b - A*x_s;
        r_bar = R*r;
        e_bar = multigrid(Abar,r_bar,max_depth-1,smooth_it);
        e = P*e_bar;
        x = x_s + e;
       
    end
end
end


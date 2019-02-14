function [u_s] = smoother(M,N,b,count,u0)

u_s = u0;
for j = 1:count
    u_s = M\(N*u_s + b);% update u_s
end
  
end




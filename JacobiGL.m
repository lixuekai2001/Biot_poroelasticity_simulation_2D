function [x w] = JacobiGL(alpha,beta,N);

% function [x] = JacobiGL(alpha,beta,N)
% Purpose: Compute the N'th order Gauss Lobatto quadrature 
%          points, x, associated with the Jacobi polynomial,
%          of type (alpha,beta) > -1 ( <> -0.5). 

x = zeros(N+1,1);
if (N==1) 
    x(1)=-1.0; x(2)=1.0; 
    w = [1; 1];
    return; 
end;

[xint,w] = JacobiGQ(alpha+1,beta+1,N-2);
x = [-1, xint', 1]';

if nargout == 2
    V = Vandermonde1D(N,x);
    M = inv(V*V');
    w = sum(M,2); % lumped GLL mass gives weights
end
return;

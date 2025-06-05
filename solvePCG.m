function [u,m] = solvePCG(A, f, u_s, tol, m_max)
% SOLVEPCG   Preconditioned Conjugate Gradients method for solving
%    a system of linear equations Au = f.
%
%    Input parameters:
%       A    : symmentric, positive definite NxN matrix
%       f    : right-hand side Nx1 column vector
%       u_s  : Nx1 start vector (initial guess)
%       C    : Preconditioner, split into C1 and C2. Splitting C 
%              in C1 und C2 offers computational advantages.
%	             Example: Symmetric Gauss-Seidel preconditioner:
%               C = (D+L)*inv(D)*(D+L')
%               D = diag(diag(A))
%              C1 = tril(A) 
%              C2 = inv(D)*triu(A) 
%       tol  : Break condition for the relative residuum
%       m_max: Maximum number of iterations to perform
%
%    Output parameters:
%	      m    : Number of actually performed iterations
%       u    : Solution vector

% Author : Andreas Klimke
% Date   : May 2003
% Version: 1.1
	
u = u_s;
r = f - A * u;
D = diag(diag(A));
C1 = tril(A) ;
C2 = D\triu(A) ; 
p = C2 \ (C1 \ r);
norm_f = norm(f);
  
m = 0;
while( (norm(r)/norm_f > tol) & (m < m_max))
  a = A * p;
  a_dot_p = a' * p;
  lambda = (r' * p) / a_dot_p;
  u = u + lambda * p;
  r = r - lambda * a;
  inv_C_times_r = C2 \ (C1 \ r);
  p = inv_C_times_r - ((inv_C_times_r' * a) / a_dot_p) * p;
  m=m+1;
end
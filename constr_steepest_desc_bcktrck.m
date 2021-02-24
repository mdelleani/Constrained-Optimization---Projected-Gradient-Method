function [xk, fk, gradfk_norm, deltaxk_norm, k] = ...
    constr_steepest_desc_bcktrck(x0, f, gradf, alpha0, ...
    kmax, tollgrad, c1, rho, btmax, gamma, tolx, Pi_X)
%
% function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq] = ...
%     constr_steepest_desc_bcktrck(x0, f, gradf, alpha0, ...
%     kmax, tollgrad, c1, rho, btmax, gamma, tolx, Pi_X)
%
% Projected gradient method (steepest descent) for constrained optimization.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% alpha0 = the initial factor that multiplies the descent direction at each
% iteration;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy.
% gamma = factor that multiplies the descent direction before the (possible) projection;
% tolx = a real scalar value characterizing the tolerance with respect to the norm of  |xk+1 - xk| in order to stop the method;
% Pi_X = projection function
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed


% Function handle for the armijo condition
farmijo = @(fk, alpha, xk, pk) ...
    fk + c1 * alpha * gradf(xk)' * pk;

% Initializations

xk = Pi_X(x0); % Project the starting point if outside the constraints
fk = f(xk);
k = 0;
gradfk_norm = norm(gradf(xk));
deltaxk_norm = tolx + 1;

while k < kmax && gradfk_norm >= tollgrad && deltaxk_norm >= tolx
    % Compute the descent direction
    pk = -gradf(xk);
    
    xbark = xk + gamma * pk;
    xhatk = Pi_X(xbark);    
    
    % Reset the value of alpha
    alpha = alpha0;
    
    % Compute the candidate new xk
    pik = xhatk - xk;
    xnew = xk + alpha * pik;
    
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo (w.r.t. pik) condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, xk, pik)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pik;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end
    
    % Update xk, fk, gradfk_norm, deltaxk_norm
    deltaxk_norm = norm(xnew - xk);
    xk = xnew;
    fk = fnew;
    gradfk_norm = norm(gradf(xk));
    
    % Increase the step by one
    k = k + 1;
   
end

end
function [gradfx] = findiff_grad(f, x, h, type)
%
% function [gradf] = findiff_grad(f, x, h, type)
%
% Function that approximate the gradient of f in x (column vector) with the
% finite difference (forward/centered) method.
%
% INPUTS:
% f = function handle that describes a function R^n->R;
% x = n-dimensional column vector;
% h = the h used for the finite difference computation of gradf
% type = 'fw' or 'c' for choosing the forward/centered finite difference
% computation of the gradient.
%
% OUTPUTS:
% gradfx = column vector (same size of x) corresponding to the approximation
% of the gradient of f in x.

gradfx = zeros(size(x));

 h = h * norm(x);

switch type
    case 'fw'
%       Without separability
%         for i=1:length(x)
%            xh = x;
%            xh(i) = xh(i) + h;
%            prova = (f(xh) - f(x))/ h;
%            prova
%            gradfx(i) = sum((f(xh) - f(x))/ h);
%  
%         end

%      EXPLOIT SEPARABILITY
        xh = x+h;
        gradfx = ((f(xh)-f(x))/h);
    case 'c'
%          Without separability
%         for i=1:length(x)
%             xh_plus = x;
%             xh_minus = x;
%             xh_plus(i) = xh_plus(i) + h;
%             xh_minus(i) = xh_minus(i) - h;
%             gradfx(i) = (f(xh_plus) - f(xh_minus))/(2 * h);
%         end
%      EXPLOIT SEPARABILITY
            xh_plus = x+h;
            xh_minus = x-h;
            gradfx = ((f(xh_plus)-f(xh_minus))/(2*h));
    otherwise % 
          % exact derivatives
         gradf = @(x) (x.^3 + x -1);
         gradfx = gradf(x);
end


end

function [y,grady] = quadobj(x,Q,f,c)

%y = -0.5*x'*Q*x+c;
y = -0.5*x'*Q*x +f'*x+c;
if nargout > 1
    grady = -Q*x +f;
    %grady=-Q*x;
end
end
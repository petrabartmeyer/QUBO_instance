function hess = quadhess(x,lambda,Q,H)
hess = Q;
jj = length(H); % jj is the number of equality constraints
for i = 1:jj
    hess = hess + lambda.eqnonlin(i)*H{i};
%    hess = hess + lambda.ineqnonlin(i)*H{i};
end

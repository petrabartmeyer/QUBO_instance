
function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
jj = length(H); % jj is the number of inequality constraints
y = zeros(1,jj);
for i = 1:jj
     yeq(i) = 1/2*x'*H{i}*x + k{i}'*x +d{i};
%    y(i) = 1/2*x'*H{i}*x +d{i};
end
y = [];
%yeq=[];
%gradyeq=[];
grady=[];
if nargout > 2
    gradyeq = zeros(length(x),jj);
%    grady=zeros(length(x),jj);
    for i = 1:jj
         gradyeq(:,i) = H{i}*x + k{i};
%          grady(:,i) = H{i}*x+k{i};
    end
end

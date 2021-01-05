 clc
clear all
 
print_solution = fopen('teste.txt','a+');
Q = textread("instance.txt");
(n,~) = size(Q);
f=ones(2*n,1);
Q=[Q, zeros(n);zeros(n),zeros(n)];
zeta=2*n;
c=0;

%%% H are the Hessian constraint
H{1}=[eye(n) -eye(n); -eye(n) eye(n)];
k{1}=zeros(zeta,1);
d{1}=-n;

options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H));


options.SubproblemAlgorithm='cg';
options.MaxIterations = 8000 ;
options.MaxFunctionEvaluations = 8000;
options.MaxPCGIter=8000;
options.TolProjCG=1e-8;
options.TolProjCGAbs=1e-8;
options.TolPCG=1e-12;
options.ConstraintTolerance=1e-8;
options.TolConSQP=1e-8;
options.FunctionTolerance=1e-6;
options.StepTolerance=1e-8;
options.OptimalityTolerance=1e-8;
options.TolCon=1e-8; 
fun = @(x)quadobj(x,Q,f,c);
nonlconstr = @(x)quadconstr(x,H,k,d);


x0 = ones(zeta,1); % column vector
tic
[x,fval,eflag,output,lambda] = fmincon(fun,x0,...
  [],[],[],[],zeros(zeta,1),ones(zeta,1),nonlconstr,options);
time=toc;

func_obj=x'*Q*x-f'*x;
feasibility=x'*H{1}*x;
fprintf(print_solution,'%d \t %1.3f \t %1.3f \t %1.3f \n',n,time,func_obj,feasibility);

 fclose(print_solution);

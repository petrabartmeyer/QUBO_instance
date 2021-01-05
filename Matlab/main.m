 clc
clear all
 
imprimir=fopen('teste.txt','a+');
 %problema=importdata('Instancias/bqp1000');
 %problema=mpsread('Instancias/bab1.mps')
 problema = textread('InstanceGintarras/in3000_4.txt');
 [nn,mm] = size(problema);
 n=sqrt(nn*mm);
 n= 3000;
 A = zeros(n,n);
 count = 1;
 for i=1:n
    for j=1:n
       A(i,j)= problema(1);
       count = count +1;
    end
 end
beta=diag(A);
Q=A;%+A' - diag(beta);
Qold=Q;
Qold=[Qold, zeros(n);zeros(n),zeros(n)];
%%% valores de f pioraram o tempo e a convergencia
f=ones(2*n,1);
Q=[Q, zeros(n);zeros(n),zeros(n)];
Q=Q+500000*diag(f);
zeta=2*n;
c=0;

%%% H são as restrições quadraticas
H{1}=[eye(n) -eye(n); -eye(n) eye(n)];
k{1}=zeros(zeta,1);
d{1}=-n;
%d{1}=n;
options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H));
%options.HessianApproximation ='bfgs';

options.SubproblemAlgorithm='cg';

%options.CheckGradients=true;
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
%[x,fval,eflag,output,lambda] = quadprog(Q,f,H,d);
[x,fval,eflag,output,lambda] = fmincon(fun,x0,...
  [],[],[],[],zeros(zeta,1),ones(zeta,1),nonlconstr,options);
tempo=toc;
func_obj=x'*Q*x-f'*x;
func_obj_old=x'*Qold*x;
factibilidade=x'*H{1}*x;
fprintf(imprimir,'%d \t %1.3f \t %1.3f \t %1.3f \n',n,tempo,func_obj_old,factibilidade);
% end
 fclose(imprimir);

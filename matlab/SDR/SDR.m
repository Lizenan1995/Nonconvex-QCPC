addpath(genpath('SDPNAL+v1.0'));
%%
N=20; % number of samples
d=3; % feature dimension
X=rand(N,d);
w_h=[5,0,0];
y_h=X*w_h';
D_h = squareform(pdist(y_h));

w_f=[2.5,0,1];
y_f=X*w_f';
y = y_f;

%%
% set up model
model = ccp_model(['SDR','QCQP']);
W = var_sdp(d+1, d+1);
model.add_variable(W);
temp1 = [X'*X,-X'*y];
temp2 = [-y'*X,y'*y];
C = [temp1;temp2];
model.minimize(inprod(C, W));

%%
% SDP equality constraints
z = zeros((d+1),(d+1));
z((d+1),(d+1)) = 1;
A = z;
model.add_affine_constraint(inprod(A, W) == 1);

%%
% SDP inequality constraints
cnt=0;
H={};
linear_coefficient={};
constant={};
for i=1:N
    for j=i+1:N
        for k=j+1:N
            sign_of_P=sign(D_h(i,k)-D_h(i,j));
            x_i=X(i,:)';
            x_j=X(j,:)';
            x_k=X(k,:)';
            P_ijk=sign_of_P*(x_k-x_j)*((2*x_i-x_j-x_k)');
            if nnz(P_ijk)~=0
                cnt=cnt+1;
                temp = zeros(d+1,d+1);
                temp(1:d,1:d) = P_ijk*2;
                H{cnt}=P_ijk*2;
                model.add_affine_constraint(inprod(temp, W) <= 0)
                linear_coefficient{cnt}=zeros(d,1);
                constant{cnt}=0;
            end
        end
    end
end

%%
% solve model
model.setparameter('maxiter',5000,'printlevel',1,'stopoption',0);
model.solve;
Ans = W.model.info.opt_solution;
wwT = Ans{1,1};
wwT = wwT((1:d),(1:d));
[V,D] = eig(wwT);
w = sqrt(D(d,d))*V(:,d);
formatSpec = 'the first stage finished';
fprintf(formatSpec)
%%
%{
cnt = 0;
violate = 0;
for i=1:length(H)
    cnt = cnt+1;
    P = H{i};
    formatSpec = 'constraint %d, value is %.4f\n';
    fprintf(formatSpec, cnt, w'*P*w);
    if(w'*P*w >0)
         violate = violate + 1;
    end
end
%}
%%
Q=2*(X'*X);
f=-2*X'*y;
c=y'*y;
w0=w;
options = optimoptions(@fmincon,'Display','iter',...
'PlotFcn',{@optimplotx,...
@optimplotfval,@optimplotfirstorderopt},...
'Algorithm','interior-point',...
'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H));
options.StepTolerance=1e-10;
problem.options=options;
problem.solver='fmincon';
problem.objective=@(x)quadobj(x,Q,f,c);
problem.x0=w0;
problem.nonlcon = @(x)quadconstr(x,H,linear_coefficient,constant);
[w,fval]=fmincon(problem);
formatSpec = 'the second stage finished';
fprintf(formatSpec);
%%

function [objs,Y,W,V,Z]=optimization_Phi_Psi_not(X,Psi,Phi,opts)    
%% parameters
K=opts.K;
lambda_1=opts.lambda_1;
lambda_2=opts.lambda_2;


rho_1=opts.rho_1;
rho_2=opts.rho_2;
MAX_ITER=opts.iter;


%% initialization
[n,t]=size(X);
[p,t]=size(Phi);
I_2=eye(n,p);


%W=randi(10,K,p);
W=zeros(K,p);
%Y=randi(10,n,K);
%Y=zeros(n,K);
Y=eye(n,K);

%Y=inv(Psi);
%Y=zeros(n,n);
V=0;
Z=Y;

Gamma_1=Y;
Gamma_2=W;

% MAX_ITER = 100;
obj_old = 0;
objs =[];



[Q_1,Lam_1]=eig(Psi'*Psi);
[Q_4,Lam_4]=eig(Phi*Phi');

obj_new  = getObj(X,Phi,Psi,Y,W);
objs=[objs,obj_new];

%start of optimization
for k = 1:MAX_ITER

    
    %Y-update
    Y =  hard_update_Y(W,X,Z,Psi,Phi,Lam_1,Q_1,Gamma_1,rho_1);
   
    %Z update
    h= Y-Gamma_1/rho_1;
    Z = sign(h).*max(abs(h)-lambda_1/rho_1,0);

    % W-update
    W=hard_update_W(V,X,Y,Psi,Phi,Lam_4,Q_4,Gamma_2,rho_2);
    
   
    % V update
    h= W-Gamma_2/rho_2;
    V = sign(h).*max(abs(h)-lambda_2/rho_2,0);


    
    % update Gamma
    Gamma_1 = Gamma_1 + rho_1*(Z-Y);
    Gamma_2 = Gamma_2 + rho_2*(V-W);
    %Y=Y-diag(Y);
    
    rho_1=min(rho_1*1.1, 1e5);
    rho_2=min(rho_2*1.1, 1e5);
    
    % stop condition
    % ++++++++++++++++++++++++++++ stop condition check
    obj_new  = getObj(X,Phi,Psi,Y,W);
    
    objs=[objs,obj_new];
    residual = abs(obj_old-obj_new); %/obj_new;
    disp(['obj-',num2str(k),'=',num2str(obj_new),',','residual-',num2str(k),'=',num2str(residual)]);
    
    if residual <  1e-4
        break;
    else
        obj_old = obj_new;
    end
    
    % update
    
    
end
% save result

end


function obj_new = getObj(X,Phi,Psi,Y,W)


% X−YΨWΦ‖

%obj_new = norm(X-Y*Psi*X)
obj_new = norm(X-Psi*Y*W*Phi)
%obj_new = norm(X-W*Phi);


end
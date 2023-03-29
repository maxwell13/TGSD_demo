function [objs,Y,Sigma,W,V,Z]=optimization_Phi_Psi_not_masked(X,Psi,Phi,opts)    
%% parameters
K=opts.K;
lambda_1=opts.lambda_1;
lambda_2=opts.lambda_2;


rho_1=opts.rho_1;
rho_2=opts.rho_2;
MAX_ITER=opts.iter;


%% initialization
[n,n2]=size(Psi);
[p,t]=size(Phi);

%W=randi(10,K,p);
W=zeros(K,p);
W(1,1)=.000001;
%Y=randi(10,n,K);
%Y=zeros(n,K);
Y=eye(n2,K);

Sigma=eye(K,K);
%Y=inv(Psi);
%Y=zeros(n,n);
V=0;
Z=Y;

Gamma_1=Y;
Gamma_2=W;

% MAX_ITER = 100;
obj_old = 0;
objs =[];

I_w=eye(size((Psi*Y)'*(Psi*Y)));

[Q_1,Lam_1]=eig(Psi'*Psi);
[Q_4,Lam_4]=eig(Phi*Phi');

% obj_new  = getObj(X,Phi,Psi,Y,W);
% objs=[objs,obj_new];

%start of optimization
for k = 1:MAX_ITER

    
        
    P=Psi*Y*W*Phi;
    D=update_d(P,X,opts.mask,opts.lambda_3);

    
    %Y-update
    Y =  hard_update_Y(Sigma,W,D,Z,Psi,Phi,Lam_1,Q_1,Gamma_1,rho_1);
   
    
    %Z update
    h= Y-Gamma_1/rho_1;
    Z = sign(h).*max(abs(h)-lambda_1/rho_1,0);

    % W-update
%     A=Psi*Y*Sigma;
%     W=inv(2*(A')*A + I_w*rho_2)*(2*A'*D*Phi'+rho_2*V+ Gamma_2);
    W=hard_update_W(Sigma,V,D,Y,Psi,Phi,Lam_4,Q_4,Gamma_2,rho_2);
    
    
   
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
     if (0==mod(k,25))
        obj_new  = getObj(opts.mask,D,X,Phi,Psi,Y,Sigma,W,lambda_1,lambda_2,opts.lambda_3)


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
    
end
% save result

end

% save('proteinRanks.mat','Gamma')

% [largest ind]= maxk(Gamma,100)
% p=sort(Gamma)
% plot(p)

function obj_new = getObj(mask,D,X,Phi,Psi,Y,Sigma,W,lambda_1,lambda_2,lambda_3)


% X−YΨWΦ‖

%obj_new = norm(X-Y*Psi*X)
temp=X(:)-D(:);
term3=norm(temp(setdiff(1:end,mask)));
obj_new = norm(D-Psi*Y*Sigma*W*Phi)+lambda_1*norm(Y,1)+lambda_2*norm(W,1)+lambda_3*term3;
%obj_new = norm(X-W*Phi);


end

function  [objs,Y,Sigma,W,V,Z] = optimization_Psi_Phi_orth(X,Psi,Phi,opts);

%% parameters
K=opts.K;
lambda_1=opts.lambda_1;
lambda_2=opts.lambda_2;

rho_1=opts.rho_1;
rho_2=opts.rho_2;
MAX_ITER=opts.iter;


%% initialization
[n,t]=size(X);
[hold,Yl]=size(Psi);
[p,t]=size(Phi);
I_2=eye(n,p);


%W=randi(10,K,p);
W=zeros(K,p);
W(1,1)=.000001;
%W(2,1)=1;
%Y=randi(10,n,K);
%Y=Y./norm(Y)
Y=zeros(Yl,K);
%Y(1,1)=1;
%Y=eye(n,K);
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



I_y=eye(size((W*Phi)*(W*Phi)'));
I_w=eye(size((Psi*Y)'*(Psi*Y)));

%[Q_1,Lam_1]=eig(Psi'*Psi);
%[Q_4,Lam_4]=eig(Phi*Phi');


obj_new  = getObj(X,Phi,Psi,Y,Sigma,W,lambda_1,lambda_2);
objs=[objs,obj_new];

% SigmaL = eye(size(Sigma));
% SigmaR = eye(size(Sigma));

XPhiT=X*Phi';
PsiTX=Psi'*X;
%start of optimization
for k = 1:MAX_ITER

    %Y-updatecd ..
    B=Sigma*W*Phi;
    Y=(2*PsiTX*B'+rho_1*Z+Gamma_1)*inv(2*(B*B')+rho_1*I_y+exp(-15));
   

    %Z update
    h= Y-Gamma_1/rho_1;
    Z = sign(h).*max(abs(h)-lambda_1/rho_1,0);
    %Z=inv(L+L'-rho_1*eye(size(L)))*(rho_1*Y-Gamma_1);

    

    A=Psi*Y*Sigma;
    W=inv(2*(A')*A + I_w*rho_2)*(2*A'*XPhiT+rho_2*V+ Gamma_2);


 
    
    % V update
    h= W-Gamma_2/rho_2;
    V = sign(h).*max(abs(h)-lambda_2/rho_2,0);

%      Sigma=Sigma.*SigmaR;
     %verifying th result updates properly
     
  %   norm(tempY*tempW-Y*Sigma*W)
 
    %Sigma=update_Sigma(Sigma,W,X,Y,Psi,Phi);
    % update Gamma
    
    
    Gamma_1 = Gamma_1 + rho_1*(Z-Y);
    Gamma_2 = Gamma_2 + rho_2*(V-W);
   
    rho_1=min(rho_1*1.1, 1e5);
    rho_2=min(rho_2*1.1, 1e5);
    
    % stop condition
    % ++++++++++++++++++++++++++++ stop condition check
    
    if (0==mod(k,25))
        obj_new  = getObj(X,Phi,Psi,Y,Sigma,W,lambda_1,lambda_2);

        objs=[objs,obj_new];
        residual = abs(obj_old-obj_new); %/obj_new;
        disp(['obj-',num2str(k),'=',num2str(obj_new),',','residual-',num2str(k),'=',num2str(residual)]);

        if residual <  1e-6;
            break;
        else
            obj_old = obj_new;
        end
    end
    
    % update
    
    
end
% save result

end


function obj_new = getObj(X,Phi,Psi,Y,Sigma,W,lambda_1,lambda_2)


% X−YΨWΦ‖

%obj_new = norm(X-Y*Psi*X)
obj_new = norm(X-Psi*Y*Sigma*W*Phi)+lambda_1*norm(Y,1)+lambda_2*norm(W,1);
%obj_new = norm(X-W*Phi);


end

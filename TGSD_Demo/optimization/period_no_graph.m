function[Wt]=period_no_graph(X,Phi,lambdat,rhot,iter)

obj_old=0
objs=[]
Gammat=0
Vt=0

XPhiT=X*Phi';
PhiPhiT=Phi*Phi';
rhotIphi=rhot*eye(size(Phi*Phi'));

for i=0:iter
   
Wt=(XPhiT+rhot*Vt+Gammat)*(PhiPhiT+rhotIphi)^(-1);

h= Wt-Gammat/rhot;
Vt = sign(h).*max(abs(h)-lambdat/rhot,0);

Gammat = Gammat + rhot*(Vt-Wt);


    if (0==mod(i,25))
    obj_new  = getObj(X,Phi,Wt,lambdat);

        objs=[objs,obj_new];
        residual = abs(obj_old-obj_new); %/obj_new;
        disp(['obj-',num2str(i),'=',num2str(obj_new),',','residual-',num2str(i),'=',num2str(residual)]);

        if residual <  1e-4
            fprintf('BREAK')
            break;
        else
            obj_old = obj_new;
        end
    end
end
end


function obj_new = getObj(X,Phi,Wt,lambdat)


% X−YΨWΦ‖

%obj_new = norm(X-Y*Psi*X)
obj_new = norm(X-Wt*Phi)+lambdat*norm(Wt,1)
%obj_new = norm(X-W*Phi);


end
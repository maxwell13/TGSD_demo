function   [objs,Y,Sigma,W,V,Z] = optimization(X,Psi,Phi,opts,type)



if(isempty(opts.mask))
    nomask=1
    %checks Psi 

    if(opts.Psi_orth==1 && opts.Phi_orth==1)
        [objs,Y,Sigma,W,V,Z]=optimization_Psi_Phi_orth(X,Psi,Phi,opts);

    elseif(opts.Psi_orth==1)
          [objs,Y,Sigma,W,V,Z]=optimization_Psi_orth_Phi_not(X,Psi,Phi,opts);

    elseif(opts.Phi_orth==1)  
        [objs,Y,Sigma,W,V,Z]=optimization_Phi_orth_Psi_not(X,Psi,Phi,opts);

    else
        [objs,Y,Sigma,W,V,Z]=optimization_Phi_Psi_not(X,Psi,Phi,opts);

    end
else 
    
    if strcmp(type,'col') | strcmp(type,'pred')
    
    
    
    [n,m]=size(X);
    temp2=ones(size(X));
    temp2(:,opts.mask)=0
    opts.mask=find(temp2==0)

    end 


    if strcmp(type,'row')

    [n,m]=size(X);
    temp2=ones(size(X));
    temp2(opts.mask,:)=0
    opts.mask=find(temp2==0)

    end

    if(opts.Psi_orth==1 && opts.Phi_orth==1)
        [objs,Y,Sigma,W,V,Z]=optimization_Psi_Phi_orth_masked(X,Psi,Phi,opts);

    elseif(opts.Psi_orth==1)
          [objs,Y,Sigma,W,V,Z]=optimization_Psi_orth_Phi_not_masked(X,Psi,Phi,opts);

    elseif(opts.Phi_orth==1)  
        [objs,Y,Sigma,W,V,Z]=optimization_Phi_orth_Psi_not_masked(X,Psi,Phi,opts);

    else
        [objs,Y,Sigma,W,V,Z]=optimization_Phi_Psi_not_masked(X,Psi,Phi,opts);

    end
    

end




function W =  hard_update_W(Sigma,V,X,Y,Psi,Phi,Lam_4,Q_4,Gamma_2,rho_2)

    A=Psi*Y*Sigma;

   [Q_3,Lam_3]=eig(A'*A);
    Pi=2*A'*X*Phi'+rho_2*V+Gamma_2;

    QPiQ=Q_3'*Pi*Q_4;
    

    
    temp0 = 2*diag(Lam_3) * diag(Lam_4)'+ rho_2;

    E=QPiQ./temp0;

    W=Q_3*E*Q_4';
  
end     
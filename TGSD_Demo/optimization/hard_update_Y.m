function Y =  hard_update_Y(Sigma,W,X,Z,Psi,Phi,Lam_1,Q_1,Gamma_1,rho_1)

    B=Sigma*W*Phi;
   [Q_2,Lam_2]=eig(B*B');
    Pi=2*Psi'*X*B'+rho_1*Z+Gamma_1;

    
    QPiQ=Q_1'*Pi*Q_2;
    

    temp0 = 2*diag(Lam_1) * diag(Lam_2)'+ rho_1;
    E=QPiQ./temp0;
    
    Y=Q_1*E*Q_2';
  
end
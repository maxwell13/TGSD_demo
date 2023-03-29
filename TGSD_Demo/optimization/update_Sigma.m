function Sigma =  update_Sigma(Sigma,W,X,Y,Psi,Phi)
D=Psi*Y;
Omega=D'*D;
G=W*Phi;
H=G*G';
L=D'*X*G';

Hpos=H.*(H>0);
Hneg= -1.*H.*(H<0);

Omegapos=Omega.*(Omega>0);
Omeganeg=-1.*Omega.*(Omega<0);

Lpos=L.*(L>0);
Lneg=-1.*L.*(L<0);

Sigma=Sigma.*(((Hneg*Sigma*Omegapos + Hpos*Sigma*Omeganeg+Lpos)./(Hpos*Sigma*Omegapos + Hneg*Sigma*Omeganeg+Lneg)).^(1/2));

Sigma=diag(diag(Sigma));
end 
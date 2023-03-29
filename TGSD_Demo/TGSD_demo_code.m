

%% version and libs 
% Signal Processing Toolbox  should be only neccsary for new DFT otherwise
% others should not be need

% MATLAB                                                Version 9.8         (R2020a)
% Simulink                                              Version 10.1        (R2020a)
% Bioinformatics Toolbox                                Version 4.14        (R2020a)
% Communications Toolbox                                Version 7.3         (R2020a)
% DSP System Toolbox                                    Version 9.10        (R2020a)
% Deep Learning Toolbox                                 Version 14.0        (R2020a)
% MET: Memory Efficient Tucker                          Version 1.0                 
% Mapping Toolbox                                       Version 4.10        (R2020a)
% Signal Processing Toolbox                             Version 8.4         (R2020a)
% Statistics and Machine Learning Toolbox               Version 11.7        (R2020a)
% Symbolic Math Toolbox                                 Version 8.5         (R2020a)
% Tensor Toolbox (Sandia National Labs)                 Version 2.6 



% add optimization to path


load('demo_data.mat')


% to generate a new GFT
%L=diag(sum(adj))-adj;
%[PsiGFT,e]= eigs(L,175);


% to generate a new DFT
% [~,t]=size(x);
% PhiDFT=(1/(t^(1/2)))*dftmtx(t);


parm.iter=500;
parm.K=7;
parm.lambda_1=.1;
parm.lambda_2=.1;
parm.lambda_3 = 1; 

parm.rho_1=parm.lambda_1/10;
parm.rho_2=parm.lambda_2/10;

%mask for missing values (mask is vecotorized location in "rand" case)
parm.mask = mask;

%change for differnt dictionaries
parm.Psi_orth=1;
parm.Phi_orth=1;


% type of mask see optimization for detials
type='rand';

%masked signal
X_masked=X;

%Sigma is unused
[objs,Y,Sigma,W,V,Z]=optimization(X_masked,PsiGFT,PhiDFT,parm,type);

%pred matrix with values reconstructed
pred_matrix = PsiGFT*Y*W*PhiDFT;











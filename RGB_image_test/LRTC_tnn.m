function [X, Out,psnrarr,ssimarr] = LRTC_tnn(B, Omega, opts)
%%% Solve the Low-Rank Tensor Completion (LRTC) based on Tensor Nuclear Norm (TNN) problem by M-ADMM
%
% min_X ||X||_*, s.t. P_Omega(X) = P_Omega(B)
%
% ---------------------------------------------
% Input:
%       M       -    d1*d2*d3 tensor
%       Omega   -    index of the observed entries
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.maxit      -   maximum number of iterations
%           opts.beta       -   stepsize for dual variable updating in ADMM
%           opts.max_beta   -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase beta
%           opts.DEBUG      -   0 or 1
%
% Output:
%       X       -    d1*d2*d3 tensor
%       Out
%         Out.PSNR
%         Out.Res
%         Out.ResT


% by Zemin Zhang 
% @INPROCEEDINGS{ZhangCVPR,
%   author={Z. Zhang and G. Ely and S. Aeron and N. Hao and M. Kilmer},
%   booktitle={the IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
%   title={Novel Methods for Multilinear Data Completion and De-noising Based on Tensor-{SVD}},
%   year={2014},
%   volume={},
%   number={},
%   pages={3842-3849},
% }
% @article{zhang2017exact,
%   title={Exact tensor completion using t-{SVD}},
%   author={Zhang, Zemin and Aeron, Shuchin},
%   journal={IEEE Transactions on Signal Processing},
%   volume={65},
%   number={6},
%   pages={1511--1526},
%   year={2017},
% }
%%

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'maxit');    maxit = opts.maxit;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'beta');        beta = opts.beta;                end
if isfield(opts, 'max_beta');      max_beta = opts.max_beta;        end

Nway = size(B);
X = B;
Y = zeros(Nway); %% ��������
M = Y; %% �������ճ���
Out.PSNR=[];Out.SSIM=[]; Out.iter=[];
psnrarr=zeros(1,maxit);
ssimarr=zeros(1,maxit);
for k = 1:maxit
    Xold = X;
    %% solve Y-subproblem
    [Y,~,~] = prox_tnn(X - M/beta,1/beta);
    %% solve X-subproblem
    X = Y + M/beta;
    X(Omega) = B(Omega);
    
    %% check the convergence

   
    res=norm(X(:)-Xold(:))/norm(X(:));
    
%     if k==1 || mod(k, 10) == 0
%         if isfield(opts, 'Xtrue')
%             fprintf('TNN: iter = %d   PSNR=%f   res=%f   real-res=%f\n', k, psnr, res, resT);
%         else
%             fprintf('TNN: iter = %d   res=%f   \n', k, res);
%         end
%     end
    if res < tol 
        break;
    end
    
    %% update Lagrange multiplier
    M = M + beta * (Y - X);
    beta = min(rho * beta, max_beta);
%     psnrarr(k)=TensorPSNR(X,opts.Xtrue);
%     ssimarr(k)=TensorSSIM(X,opts.Xtrue);
end
  Out.iter=k;
  XT=opts.Xtrue;
   Out.PSNR=TensorPSNR(X,XT);
   Out.SSIM=TensorSSIM(X,XT);
 

end
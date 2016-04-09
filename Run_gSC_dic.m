% Author:
% - Mehrtash Harandi (mehrtash.harandi at gmail dot com)
%
% This file is provided without any warranty of
% fitness for any purpose. You can redistribute
% this file and/or modify it under the terms of
% the GNU General Public License (GPL) as published
% by the Free Software Foundation, either version 3
% of the License or (at your option) any later version.

clear;
clc;
dbstop error;
Solver_Flag = 1;  %1: SPAMS, 2: CVX
%SPAMS toolbox is available from http://spams-devel.gforge.inria.fr/
%CVX is available from http://cvxr.com/cvx/


SR_lambda = 1e-1;    %sparse representation parameter
nAtoms = 128;        %size of the dictionary
dict_options.L = 20; %number of non-zero elements in OMP for dictionary learning

load('Cambridge');
trnIndex = find(grassSet.idx == 5);
tstIndex = find(grassSet.idx < 5);
trn.X = double(grassSet.X(:,:,trnIndex));
trn.y = grassSet.y(trnIndex);
tst.X = double(grassSet.X(:,:,tstIndex));
tst.y = grassSet.y(tstIndex);

addpath('iccv_ext_func')
fprintf('Learning the Grassmannian dictionary\n');
D = grassmann_dictionary_learning(trn.X,nAtoms,dict_options);
[gSC_alpha_trn,~,~] = gsc_func(trn.X,D,SR_lambda,Solver_Flag);
[gSC_alpha_tst,~,~] = gsc_func(tst.X,D,SR_lambda,Solver_Flag);

% %Classification
%Ridge regression
nClasses = max(trn.y);
nPoints = length(trn.y);
L = zeros(nClasses,nPoints);
L(sub2ind([nClasses,nPoints], trn.y, 1:nPoints)) = 1;
zeta = 1e-1;    %regularization parameter for ridge regression
big_alpha = gSC_alpha_trn*gSC_alpha_trn' + zeta*eye(nAtoms);
big_v = L*gSC_alpha_trn';
W = big_v/big_alpha;

[~,y_hat] = max(W*gSC_alpha_tst);
CRR = sum(y_hat == tst.y)/length(y_hat);
fprintf('Correct recognition accuracy with a dictionary of size %d : %.1f%%.\n',nAtoms,100*CRR);








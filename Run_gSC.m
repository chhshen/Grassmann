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


Solver_Flag = 1;  %1: SPAMS, 2: CVX
%SPAMS toolbox is available from http://spams-devel.gforge.inria.fr/
%CVX is available from http://cvxr.com/cvx/


SR_lambda = 1e-1;   %sparse representation parameter



load('Cambridge');

trnIndex = find(grassSet.idx == 5);
tstIndex = find(grassSet.idx < 5);
trn.X = double(grassSet.X(:,:,trnIndex));
trn.y = grassSet.y(trnIndex);
tst.X = double(grassSet.X(:,:,tstIndex));
tst.y = grassSet.y(tstIndex);

addpath('iccv_ext_func');
[gSC_alpha,gSC_qX,gSC_D] = gsc_func(tst.X,trn.X,SR_lambda,Solver_Flag);
%Classification-SRC
y_hat = Classify_SRC(gSC_D,trn.y,gSC_alpha,gSC_qX);
CRR = sum(double(y_hat == tst.y))/length(tst.y);
fprintf('Correct recognition accuracy with a labeled dictionary : %.1f%%.\n',100*CRR);







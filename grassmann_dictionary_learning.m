% Author:
% - Mehrtash Harandi (mehrtash.harandi at gmail dot com)
%
% This file is provided without any warranty of
% fitness for any purpose. You can redistribute
% this file and/or modify it under the terms of
% the GNU General Public License (GPL) as published
% by the Free Software Foundation, either version 3
% of the License or (at your option) any later version.

function D = grassmann_dictionary_learning(X,nAtoms,dict_options)

if (~isfield(dict_options,'L'))
    dict_options.L = 10;
end
if (~isfield(dict_options,'nIter'))
    dict_options.nIter = 5;
end
%initializing dictionary
% use the following two lines if you want to initialize the
% dictionary randomly.
% p = randperm(size(X,3));
% D = X(:,:,p(1:nAtoms));
fprintf('Initializing the dictionary using kmeans algorithm.\n');
D = kmeans_projection(X,nAtoms,5,false);

for tmpIter = 1:dict_options.nIter
    %sparse coding
    alpha = local_sparse_coding(X,D,dict_options.L);
    Dn = update_dict(X,D,alpha);
    if(tmpIter == 1)
        preCost = compute_dic_cost(X,D,alpha);
        fprintf('Initial cost -->%.3f\n',preCost);
    end
    postCost = compute_dic_cost(X,Dn,alpha);
    fprintf('Iter#%d: cost -->%.3f\n',tmpIter,postCost);
    D = Dn;
end
fprintf('-------\n');

end
%--------------------------------------------------------------------------
function Dn = update_dict(X,D,alpha)
sym_mat = @(X) real(0.5*(X+X'));
nAtoms = size(D,3);
[n,p,~] = size(X);
Dn = zeros(size(D));
for r = 1:nAtoms
    S = zeros(n);
    idx_alpha = find(alpha(r,:));
    if (isempty(idx_alpha))
        fprintf('a useless atom identified!\n');
        continue;
    end
    for tmpC1 = 1:length(idx_alpha)
        S = S + alpha(r,idx_alpha(tmpC1))*X(:,:,idx_alpha(tmpC1))*X(:,:,idx_alpha(tmpC1))';
        for j = 1:nAtoms
            if (j == r) || (alpha(j,idx_alpha(tmpC1)) == 0)
                continue;
            end
            S = S - alpha(j,idx_alpha(tmpC1))*alpha(r,idx_alpha(tmpC1))*D(:,:,j)*D(:,:,j)';
        end
    end
    S = sym_mat(S);
    [Dn(:,:,r),~] = eigs(S,p,'LA');
end


end

%--------------------------------------------------------------------------
function cost = compute_dic_cost(X,D,alpha)
[~,p,nPoints] = size(X);
k_DD = p - 0.5*grassmann_proj_dist(D);
k_DX = p - 0.5*grassmann_proj_dist(X,D);
cost = p*nPoints -2*trace(alpha'*k_DX) + trace(alpha'*k_DD*alpha);
end

%--------------------------------------------------------------------------
function centers = kmeans_projection(X,k,nIter,verbatim_flag)
% Initializations
MinCostVariation = 1e-3;

nPoints = size(X,3);

randVal = randperm(nPoints);
centers = X(:,:,randVal(1:k));
for iter = 1:nIter
    fprintf('.');
    %assign points and compute the cost
    [currCost,minIdx] = kmeans_cost(X,centers);
    for tmpC1 = 1:k
        idx = find(minIdx == tmpC1);
        if (isempty(idx))
            %zombie centers
            randVal = randperm(nPoints);
            centers(:,:,tmpC1) = X(:,:,randVal(1));
        else
            centers(:,:,tmpC1) = grassmann_mean_proj(X(:,:,idx));
        end
    end
    if (iter == 1)
        if (verbatim_flag)
            fprintf('kmeans: initial cost is %6.1f.\n',currCost);
        end
    else
        cost_diff = norm(preCost - currCost) ;
        if (cost_diff < MinCostVariation)
            if (verbatim_flag)
                fprintf('kmeans: done after %d iterations due to small relative variations in cost.\n',iter);
            end
            break ;
        elseif (verbatim_flag)
            fprintf('kmeans: Iter#%d, cost is %6.1f.\n',iter,currCost);
        end
    end
    preCost = currCost;
end
fprintf('\n');

end



%--------------------------------------------------------------------------
function [outCost,minIdx] = kmeans_cost(X,centers)
l_dist = grassmann_proj(X,centers);
[minDist,minIdx] = max(l_dist);
outCost = sum(minDist);
end

%--------------------------------------------------------------------------
function outMean = grassmann_mean_proj(X)
%this function computes the mean of a set of points
%{X_i}_{i=1}^{nPoints} on G(p,n).

[n,p,nPoints] = size(X);
tmpBig = zeros(n,n);
for tmpC1 = 1:nPoints
    tmpBig = tmpBig + X(:,:,tmpC1)*X(:,:,tmpC1)';
end
[outMean,~] = eigs(tmpBig,p);
end

%--------------------------------------------------------------------------

function dist_p = grassmann_proj(SY1,SY2)

MIN_THRESH = 1e-6;

same_flag = false;
if (nargin < 2)
    SY2 = SY1;
    same_flag = true;
end
p = size(SY1,2);


[~,~,number_sets1] = size(SY1);
[~,~,number_sets2] = size(SY2);

dist_p = zeros(number_sets2,number_sets1);

if (same_flag)
    %SY1 = SY2
    for tmpC1 = 1:number_sets1
        Y1 = SY1(:,:,tmpC1);
        for tmpC2 = tmpC1:number_sets2
            tmpMatrix = Y1'* SY2(:,:,tmpC2);
            tmpProjection_Kernel_Val = sum(sum(tmpMatrix.^2));
            if (tmpProjection_Kernel_Val < MIN_THRESH)
                tmpProjection_Kernel_Val = 0;
            elseif (tmpProjection_Kernel_Val > p)
                tmpProjection_Kernel_Val = p;
            end
            dist_p(tmpC2,tmpC1) = tmpProjection_Kernel_Val;
            dist_p(tmpC1,tmpC2) = dist_p(tmpC2,tmpC1);
        end
    end
else
    for tmpC1 = 1:number_sets1
        Y1 = SY1(:,:,tmpC1);
        for tmpC2 = 1:number_sets2
            tmpMatrix = Y1'* SY2(:,:,tmpC2);
            tmpProjection_Kernel_Val = sum(sum(tmpMatrix.^2));
            if (tmpProjection_Kernel_Val < MIN_THRESH)
                tmpProjection_Kernel_Val = 0;
            elseif (tmpProjection_Kernel_Val > p)
                tmpProjection_Kernel_Val = p;
            end
            dist_p(tmpC2,tmpC1) = tmpProjection_Kernel_Val;
        end
    end
end


end

%--------------------------------------------------------------------------
function alpha = local_sparse_coding(X,dicX,L)
[~,p,~] = size(dicX);
% nPoints = size(X,3);

%Creating kernel matrices
K_D = grassmann_proj(dicX);
K_XD = grassmann_proj(X,dicX);


%preparing for vector calculations
[KD_U,KD_D,~] = svd(K_D);
D = diag(sqrt(diag(KD_D)))*KD_U'/sqrt(p);
D_Inv = KD_U*diag(1./sqrt(diag(KD_D)));
qX = D_Inv'*K_XD/sqrt(p);
alpha = full(mexOMP(qX,D,struct('L',L)));

end



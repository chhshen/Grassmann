% Author:
% - Mehrtash Harandi (mehrtash.harandi at gmail dot com)
%
% This file is provided without any warranty of
% fitness for any purpose. You can redistribute
% this file and/or modify it under the terms of
% the GNU General Public License (GPL) as published
% by the Free Software Foundation, either version 3
% of the License or (at your option) any later version.

function [alpha,qX,D] = gsc_func(X,dicX,SR_lambda,Solver_Flag)

[~,p,nAtoms] = size(dicX);
nPoints = size(X,3);

%Creating kernel matrices
% K_D = p*ones(nAtoms);
% K_XD = zeros(nAtoms,nPoints);
% for tmpC1 = 1:nAtoms
%     for tmpC2 = tmpC1+1:nAtoms
%         K_D(tmpC2,tmpC1) = norm(dicX(:,:,tmpC1)'*dicX(:,:,tmpC2), 'fro')^2;
%         K_D(tmpC1,tmpC2) = K_D(tmpC2,tmpC1);
%     end
% end
% for tmpC1 = 1:nPoints
%     for tmpC2 = 1:nAtoms
%         K_XD(tmpC2,tmpC1) = norm(X(:,:,tmpC1)'*dicX(:,:,tmpC2), 'fro')^2;
%     end
% end
K_D = grassmann_proj(dicX);
K_XD = grassmann_proj(X,dicX);

%preparing for vector calculations
[KD_U,KD_D,~] = svd(K_D);
D = diag(sqrt(diag(KD_D)))*KD_U';
D_Inv = KD_U*diag(1./sqrt(diag(KD_D)));
qX = D_Inv'*K_XD;

switch Solver_Flag
    case 1
        alpha = full(mexLasso(qX,D,struct('mode',2,'lambda',SR_lambda,'lambda2',0)));
    otherwise 
        alpha = zeros(nAtoms,nPoints);       
        for tmpC1 = 1:nPoints 
            cvx_begin quiet;
            variable alpha_cvx(nAtoms,1);
            minimize( norm(qX(:,tmpC1) - D * alpha_cvx) +  SR_lambda*norm(alpha_cvx,1));
            cvx_end;
            alpha(:,tmpC1) = double(alpha_cvx);
        end
         
    
end

end

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
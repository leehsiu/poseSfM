function [outputArg1,outputArg2] = EDMNRSfM(C,m,c1mask,c1val,lamba1,lambda2)
% Author: Xiu Li, Tsinghua University
%
% Input:  IDX -- K-NNG neighboring point index
%         C   -- matrix of cosine of viewing angles
%         m   -- normalized image measurement
% Output: z   -- recovered depth
%         lgs -- recovered length of legs (see our paper)
%       edges -- recovered edges (see our paper)
%     maxEdgs -- recovered edges in the internal model (maximum distance model)
if(nargin<5)
    lambda2 = 20;
end
if(nargin<4)
    lambda1 = 1;
end


F = length(C); % number of frames
[N, k]= size(IDX); % number of points, number of neighbors

tmpI = speye(N);
A = struct([]);
tmp = sparse(N*k,N*N);

%for frame 1->F
%for point 1->N
%for Neighbor 1->K

for f = 1:F
    for i = 1:N
        for j = 1:k
            e = diag(tmpI(:,i)+tmpI(:,IDX(i,j)));
            eCe = e*C(f).c'*e;
            eCe = eCe(:);
            tmp((i-1)*k+j,:) = eCe;
        end
    end
    A(f).vecA = tmp;
end


% Input:  IDX -- K-NNG neighboring point index
%         C   -- matrix of cosine of viewing angles
%         m   -- normalized image measurement
% Output: z   -- recovered depth
%         lgs -- recovered length of legs (see our paper)
%       edges -- recovered edges (see our paper)   
%     maxEdgs -- recovered edges in the internal model (maximum distance model)


%  edges    == l
%  maxEdges == d 
%  lgs      == x (N points F frames).
%



%Constraint SfM 
% d.* mask = const 

cvx_begin
    cvx_precision low 
    variables x(N*F) d(N*k) l(N*k,F) Y(N,N,F);
    tmp = trace(Y(:,:,1));
    for i = 2:F
        tmp = tmp+trace(Y(:,:,i));
    end    
    minimize( tmp-lambda1*sum(x)-lambda2*sum(vec(l)) )        
    subject to
    sum(d) == 1; % maxEdges Normalization to 1
    x >= 0;      % x are lgs
    l >= 0;      % l are edges
    l <= d*ones(1,F);   
    for f = 1:F    
        idx1 = N*(f-1)+1;
        idx2 = N*f;
        [1, x(idx1:idx2)';x(idx1:idx2), Y(:,:,f)] == semidefinite(N+1);          
        (A(f).vecA)*vec(Y(:,:,f)) == l(:,f);  %Cosine * Lgs == Edges
    end
cvx_end


lgs = x;
edgs = l;
maxEdgs = d;

lgs = reshape(lgs, N, F);
lgs = lgs';
z = zeros(F,N);

for i = 1:F
    qi = [m(i).m; ones(1,N)];
    z(i,:) = lgs(i,:)./sqrt(sum(qi.^2));
end

end

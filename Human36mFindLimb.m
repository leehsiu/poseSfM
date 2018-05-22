clear,
close all;
load h36m;
[F,P] = size(SS);
F = F/3;

allEDM = zeros(P,P,F);

for i=1:F
    allEDM(:,:,i) = squareform(pdist(SS(i*3-2:i*3,:)'));
end

varEDM = var(allEDM,0,3);

constMask = varEDM;
constMask(varEDM<1e-20) = 1;
constMask(varEDM>1e-20) = 0;


constVal = allEDM(:,:,17).*constMask;

%constant constraints mask for constEDMSfM 
save('h36const','constMask','constVal');

allPose = dir(fullfile('../EDMPose/D3_Positions_mono/','*.cdf'));
pNum = length(allPose);

%Data preparation
dN = randi(pNum);
dN = 105;

allPoseCdf = cdfread(fullfile(allPose(dN).folder,allPose(dN).name));
allPoseMat = cell2mat(allPoseCdf);
seqLen = length(allPoseMat); 

%Intended Length
sampleLen = 400;
sampleStart = randi(seqLen-sampleLen);
W = allPoseMat(sampleStart:1:sampleStart+sampleLen-1,:);

SS = zeros(sampleLen*3,32);
for k = 1:sampleLen
    Pgt = reshape(W(k,:),[3 32]);
    cPgt = Pgt - mean(Pgt,2);
    SS(k*3-2:k*3,:) = cPgt;
    p2d = cPgt;
    p2d = p2d./(repmat(p2d(3,:),[3 1]));
    
%Visualization,unnecessary
%     plot3(-Pgt(1,:), Pgt(3,:), -Pgt(2,:),'bo');
%     axis equal;    
%     figure(2);
%     plot(p2d(1,:), p2d(2,:),'bo');
%     axis image;
%     axis([-0.5 0.5 -0.5 0.5]);
%     pause(0.01);
end


sCano = SS(1:3,:);
pL = pdist(sCano');
pLm = max(pL);
SS = SS/pLm*6;

%Change the x,y
St = SS;
SS(1:3:end,:) = St(2:3:end,:);
SS(2:3:end,:) = St(1:3:end,:);
SS(3:3:end,:) = St(3:3:end,:);



save('h36m.mat','SS');
%Create Data


for k = 1:sampleLen
    Pgt = SS(k*3-2:k*3,:);
    plot3(Pgt(1,:),Pgt(2,:),Pgt(3,:),'bo');
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    pause(0.01);
    
end


% for dN = 1:pNum
%     allPos = cdfread(fullfile(allPose(dN).folder,allPose(dN).name));
%     allPosMat = cell2mat(allPos);
%     seqLen = length(allPosMat);
%     for i=1:seqLen
%         cPos = allPosMat(i,:);
%         cPos = reshape(cPos,[3,32]);
%         cEDM = squareform(pdist(cPos'));
%         allEDM{nEDM} = cEDM(:);
%         nEDM= nEDM+1;
%     end   
% end
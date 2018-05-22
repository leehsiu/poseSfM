clear,
close all

%load dance  h36m data
load h36m
load h36mconst


S = SS;


offset = min(vec(S(3:3:end,:)));

S(3:3:end,:) = S(3:3:end,:) - offset;


[D, N] = size(S);
totframes = D/6;



for k = 1:totframes    
    rng(100*k);
    %t = 2*rand(3,1);
    t = 4*rand(3,1);
    S((3*k-2):(3*k),:) = S((3*k-2):(3*k),:) + t*ones(1,N);
    Pgth(k).P = S((3*k-2):(3*k),:);
    m(k).m = Pgth(k).P(1:2,:)./(ones(2,1)*Pgth(k).P(3,:));  
end

xmin = min(vec(S(1:3:end,:)));
xmax = max(vec(S(1:3:end,:)))+2;
ymin = min(vec(S(2:3:end,:)));
ymax = max(vec(S(2:3:end,:)))+2;
zmin = min(vec(S(3:3:end,:)));
zmax = max(vec(S(3:3:end,:)))+2;

sv = 2; % subsample views
Pgt = Pgth(sv:sv:end);
m = m(sv:sv:end);


N = length(m(1).m);
M = length(m);

% visibility is true for the example:

lambda1 = 1;
lambda2 = 20;

%IDX = getNeighborsVis(m,Kneighbors,visibt);
%visbc = num2cell(visibt,1);

%IDX to EDM Mask to make symmetry

C = getAngleCosFull(m);
%Get Full Cosine Matrix.



%Visualization of the cosine matrix
%figure(1);
%for i=1:100
%    imagesc(C(i).c);   
%    pause(0.1);
%end

%%

mc  = squeeze(struct2cell(m));
disp('NRSfM function');
tic;

%Add limb-length constraints.


[mu,D] = ConstNRSfM(IDX, C, m, lambda1, lambda2);
ts = toc;

xmin = xmin+offset;
xmax = xmax+offset;
ymin = ymin+offset;
ymax = ymax+offset;
zmin = zmin+offset;
zmax = zmax+offset;

% third part: display results, 
%%
res.Q2 = cell(1,M);
res.Pg = cell(1,M);
res.err3d = zeros(1,M); % RMSE for each surface
res.err3dper = zeros(1,M); % RMSE for each surface
for k=1:M
    Q2k=double([mu(k,visibt(:,k));mu(k,visibt(:,k));mu(k,visibt(:,k))]).*[m(k).m(:,visibt(:,k));ones(1,length(m(k).m(:,visibt(:,k))))];    
    P2 = Pgt(k).P(:,visibt(:,k));
    % get valid indices: some groundtruth points are 0
    mugth = P2(3,:);
    l = mugth>0;
    % fix scale of reconstructed surface
    Q2k_n = RegisterToGTH(Q2k(:,l),P2(:,l));
    
    P2(:,l) = P2(:,l)+offset;
    Q2k_n = Q2k_n+offset;
    
    figure(1)
    clf;
    plot3(Q2k_n(1,:),Q2k_n(2,:),Q2k_n(3,:),'b*');
    hold on;
    plot3(P2(1,l),P2(2,l),P2(3,l),'go');
    hold off;
    axis equal;
    axis( [xmin xmax ymin ymax zmin zmax] );
    grid on;
    pause(0.1);
    res.Q2{k} = Q2k_n;
    res.Pg{k} = P2(:,l);
 
    scale = norm(P2(:,l),'fro');    
    res.err3dper(k) = norm(Q2k_n - P2(:,l),'fro')/scale*100;    
    res.err3d(k) = sqrt(mean(sum((P2(:,l)-Q2k_n).^2)));
    fprintf('3D rmse =%.2f mm\t',res.err3d(k));
    fprintf('relative 3D error =%.2f %% \n',res.err3dper(k));          
end
meandepth = mean(res.err3d)
meanper = mean(res.err3dper)

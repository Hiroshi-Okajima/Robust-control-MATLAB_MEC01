clc;
clear;
close all;

kMax = 100; % Number of updates
numParticles = 50; % Number of particles
pfc = true; % If use pfc, pfc = true
omega = 1; % Inartia cofficient
c1 = 0.9; % Weight of the tarm directed to parsonal best
c2 = 0.9; % Weight of the tarm directed to global best

% Load plant and model
plantData = load('./plants/stable_mp.mat'); % Minimum phase plant
% plantData = load('./plants/stable_nmp.mat'); % Non minimum phase plant

% Detarmine the structure of the differential compensator
% NaN：Design variable
% Number(e.g. 0,1...)：Parameter of the differential compensator as it is
z = NaN;
rAd = [0 1 0 0;
    0 z 0 z;
    0 0 0 1;
    0 z 0 z
    ];
rBd = [0 0;
    1 z;
    0 0;
    0 1
    ];
rCd = [z z z z;
    z z z z
    ];
rDd = [z z;z z];

% Number of parameters to design
numParams = getNumParams(rAd,rBd,rCd,rDd);

% Range of initial values
posRange = [-10 10];
veloRange = [-5 5];

if ~pfc
    A = plantData.A;
    B = plantData.B;
    C = plantData.C;
    Bw = plantData.Bw;
    Dw = plantData.Dw;
    Am = plantData.Am;
    Bm = plantData.Bm;
    Cm = plantData.Cm;
else
    sizeA = size(plantData.A);
    sizeB = size(plantData.B);
    sizeC = size(plantData.C);
    sizeAf = size(plantData.Af);
    sizeBf = size(plantData.Bf);
    sizeCf = size(plantData.Cf);
    A = zeros(sizeA(1)+sizeAf(1),sizeA(2)+sizeAf(2),sizeA(3));
    B = zeros(sizeB(1)+sizeBf(1),sizeB(2),sizeB(3));
    C = zeros(sizeC(1),sizeC(2)+sizeCf(2),sizeC(3));
    for i = 1:size(plantData.A,3)
        A(:,:,i) = [plantData.A(:,:,i) zeros(sizeA(1),sizeAf(2))
                    zeros(sizeAf(1),sizeA(2)) plantData.Af];
        B(:,:,i) = [plantData.B(:,:,i); plantData.Bf];
        C(:,:,i) = [plantData.C(:,:,i) plantData.Cf];
    end
    Bw = [plantData.Bw; zeros(sizeBf(1:2))];
    Dw = plantData.Dw;
    Am = [plantData.Am zeros(size(plantData.Am,1),sizeAf(2))
          zeros(sizeAf(1),size(plantData.Am,2)) plantData.Af];
    Bm = [plantData.Bm; plantData.Bf];
    Cm = [plantData.Cm plantData.Cf];
end

Pm = ss(Am,Bm,Cm,0);
    
globalBest.eval = 1e+9;
globalBest.pos = [];
globalBest.rAd = rAd;
globalBest.rBd = rBd;
globalBest.rCd = rCd;
globalBest.rDd = rDd;

for k=0:kMax
    i = 1;
    while i<=numParticles
        if k==0
            % Initialize
            p = struct('pos',[],'velo',[],'omega',omega,'c1',c1,'c2',c2,'rAd',rAd,'rBd',rBd,'rCd',rCd,'rDd',rDd,'eval',[],'pBest',[],'pBestEval',[]);
            p.pos = posRange(1)+(posRange(2)-posRange(1))*rand(1,numParams);
            p.velo = veloRange(1)+(veloRange(2)-veloRange(1))*rand(1,numParams);
            if ~pfc
                p = culcEval(p,A,B,Bw,C,Dw,Pm);
            else
                p = culcEval(p,A,B,Bw,C,Dw,Pm,pfc,size(plantData.A,1));
            end
            if isinf(p.eval)
                % Initialize again if generarized plant is unstable
                continue
            else
                p.pBest = p.pos;
                p.pBestEval = p.eval;
                particles(i) = p;
            end
        else
            % Evaluate
            if ~pfc
                particles(i) = culcEval(particles(i),A,B,Bw,C,Dw,Pm);
            else
                particles(i) = culcEval(particles(i),A,B,Bw,C,Dw,Pm,pfc,size(plantData.A,1));
            end
        end
        if particles(i).pBestEval < globalBest.eval && ~isinf(particles(i).eval)
            % Update global best
            globalBest.eval = particles(i).pBestEval;
            globalBest.pos = particles(i).pBest;
        end
        % Update parameters
        particles(i) = update(particles(i),globalBest);
        i=i+1;
    end
end

[D,~,~,~,~] = getD(globalBest)
evalGamma = globalBest.eval

if ~pfc
    save(['./results/pso_',char(datetime('now','Format','yyyy-MM-dd''T''HHmmss')),'.mat'])
else
    save(['./results/pso_pfc_',char(datetime('now','Format','yyyy-MM-dd''T''HHmmss')),'.mat'])
end

 function num = getNumParams(A,B,C,D)
    num = nnz(isnan(A));
    num = num + nnz(isnan(B));
    num = num + nnz(isnan(C));
    num = num + nnz(isnan(D));
 end
 
 function obj = update(obj,globalBest)
% Update equations of position and velocity
obj.velo = obj.omega*obj.velo + ...
          obj.c1*rand*(obj.pBest-obj.pos) + ...
          obj.c2*rand*(globalBest.pos-obj.pos);
obj.pos = obj.pos + obj.velo;
end

 function obj =  culcEval(obj,A,B,Bw,C,Dw,Pm,pfc,numPlantState)
% Evaluation function

if ~exist('pfc','var')
    pfc = false;
    numPlantState = 0;
end

if pfc
    numPFCState = size(A,1) - numPlantState;
end
    
[~,Ad,Bd,Cd,Dd] = getD(obj);
[Am,Bm,Cm,~] = ssdata(Pm);

for v = 1:size(A,3)
    % generalized plant
    dA=A(:,:,v)-Am; dB=B(:,:,v)-Bm; dC=C(:,:,v)-Cm;
    barAt = [          A(:,:,v)-B(:,:,v)*Dd*C(:,:,v)                     -B(:,:,v)*Cd dA-B(:,:,v)*Dd*dC;
                                         Bd*C(:,:,v)                               Ad             Bd*dC;
              zeros(size(Am,1),size(A(:,:,v),2)) zeros(size(Am,1),size(Ad,2))               Am];
    barBt = [                              Bw                  -B(:,:,v)*Dd*Dw                                 dB;
                 zeros(size(Bd,1),size(Bw,2))                            Bd*Dw zeros(size(Bd,1),size(B(:,:,v),2));
             zeros(size(Bm,1),size(Bw,2)) zeros(size(Bm,1),size(Dw,2))                            Bm];
    if ~pfc
        barEt = [C(:,:,v) zeros(size(C(:,:,v),1),size(Ad,2)) dC];
    else
        barEt = [C(:,1:numPlantState),zeros(size(C,1),numPFCState),zeros(size(C,1),size(Ad,2)),dC(:,1:numPlantState),zeros(size(C,1),numPFCState)];
    end

    if v == 1
        barA = barAt; barB = barBt; barE = barEt;
    else
        barA = cat(3,barA,barAt);
        barB = cat(3,barB,barBt);
        barE = cat(3,barE,barEt);
    end
    if ~all(eig(barAt)<0)
        % If generalized plant is unstable
        obj.eval = Inf;
        return
    end
end

setlmis([])
[gamma2,~,~] = lmivar(1,[1 1]);

[X,~,~] = lmivar(1,[size(barA(:,:,v),1) 1]);
Sn = newlmi;
lmiterm([-Sn 1 1 X],1,1)

for v=1:size(A,3)
    lmi = newlmi;
    lmiterm([lmi 1 1 X],barA(:,:,v),1,'s')
    lmiterm([lmi 2 1 X],barE(:,:,v),1)
    lmiterm([lmi 1 3 0],barB(:,:,v))
    lmiterm([lmi 2 2 gamma2],-1,eye(size(barE(:,:,v),1),size(barE(:,:,v),1)))
    lmiterm([lmi 3 3 0],-eye(size(barB(:,:,v),2),size(barB(:,:,v),2)))
end

% Solve LMI
LMIs = getlmis;
ndec = decnbr(LMIs);
c = [1 zeros(1,ndec-1)];
[~,xopt] = mincx(LMIs,c,[1e-16,1000,10000,0,1]);

if ~isempty(xopt)
    gamma2opt = dec2mat(LMIs,xopt,gamma2);
    obj.eval = sqrt(gamma2opt);
    % Compare personal best
    if obj.pBestEval >= obj.eval
        obj.pBestEval = obj.eval;
        obj.pBest = obj.pos;
    end
else
    obj.eval = Inf;
end
end

function [D,Ad,Bd,Cd,Dd] = getD(obj)
% Transrate particle to state space differential compensator
x = obj.pos;
index = 1;
Ad = zeros(size(obj.rAd));Bd=zeros(size(obj.rBd));Cd=zeros(size(obj.rCd));Dd=zeros(size(obj.rDd));
for i=1:size(obj.rAd,1)
    for j=1:size(obj.rAd,2)
        if isnan(obj.rAd(i,j))
            Ad(i,j) = x(index);
            index=index+1;
        else
            Ad(i,j) = obj.rAd(i,j);
        end
    end
end
for i=1:size(obj.rBd,1)
    for j=1:size(obj.rBd,2)
        if isnan(obj.rBd(i,j))
            Bd(i,j) = x(index);
            index=index+1;
        else
            Bd(i,j) = obj.rBd(i,j);
        end
    end
end
for i=1:size(obj.rCd,1)
    for j=1:size(obj.rCd,2)
        if isnan(obj.rCd(i,j))
            Cd(i,j) = x(index);
            index=index+1;
        else
            Cd(i,j) = obj.rCd(i,j);
        end
    end
end
for i=1:size(obj.rDd,1)
    for j=1:size(obj.rDd,2)
        if isnan(obj.rDd(i,j))
            Dd(i,j) = x(index);
            index=index+1;
        else
            Dd(i,j) = obj.rDd(i,j);
        end
    end
end
D = ss(Ad,Bd,Cd,Dd);
end

clear;
close all;
clc;

numMaxItr = 10;
usePFC = true; % Available in the case that plant is 'nmp_plant.mat' and Ad=4x4, Bd=4x2, Cd=2x4, Dd=2x2.

% Option for mincx. See mincx documentation for details.
options = [1e-16,1000,10000,0,1]; 

% Give initial value of differential compensator
initialValueFilePath = './results/pso_pfc_2022-03-16T203233.mat';
load(initialValueFilePath,'D')
[Ad,Bd,Cd,Dd] = ssdata(D);

% Load plant and model
load('./plants/stable_mp.mat') % minimum phase plant
% load('./plants/stable_nmp.mat') % non minimum phase plant

if usePFC
    pfc = ss(Af,Bf,Cf,0);
else
    pfc = NaN;
end

ssVertices = cell(4,1);
for i =1:4
    ssVertices{i} = ss(A(:,:,i),B(:,:,i),C(:,:,i),0);
end
model = ss(Am,Bm,Cm,0);

prevGamma = 1e+10;
gamma = 0;
itr = 0;

search_xLog = cell(numMaxItr,3);
search_dLog = cell(numMaxItr,3);
ssD = ss(Ad,Bd,Cd,Dd);

xopt = [];
dopt = [];

while (itr < numMaxItr && checkContinueProcess(itr,search_dLog,gamma,1,1e-6))
itr = itr + 1;
prevGamma = gamma;

if usePFC
    [gamma,lyp,xopt] = search_x_pfc(ssVertices,Bw,Dw,model,ssD,pfc,xopt,options);
else
    [gamma,lyp,xopt] = search_x(ssVertices,Bw,Dw,model,ssD,xopt,options);
end

search_xLog{itr,1} = gamma;
search_xLog{itr,2} = lyp;
search_xLog{itr,3} = ssD;
if gamma == -1
    warning(['Process completed with infeasible LMI(search_x). '...
        'This may be due to errors in numerical calculations '...
        '(change options(3) to solve) or unstable generalized plant '...
        'with the initial value.'])
    break
end

if usePFC
    [gamma,ssD,dopt] = search_d_pfc(ssVertices,Bw,Dw,model,lyp,pfc,dopt,options);
else
    [gamma,ssD,dopt] = search_d(ssVertices,Bw,Dw,model,lyp,dopt,options);
end
search_dLog{itr,1} = gamma;
search_dLog{itr,2} = ssD;
search_dLog{itr,3} = lyp;
if gamma == -1
    warning(['Process completed with infeasible LMI(search_d). '...
        'This may be due to errors in numerical calculations '...
        '(change options(3) to solve) or unstable generalized plant '...
        'with the initial value.'])
    break
end
end

D = search_dLog{itr,2}
evalGamma = search_dLog{itr,1}

if ~usePFC
    save(['./results/itr_',char(datetime('now','Format','yyyy-MM-dd''T''HHmmss')),'.mat'])
else
    save(['./results/itr_pfc_',char(datetime('now','Format','yyyy-MM-dd''T''HHmmss')),'.mat'])
end

function cont = checkContinueProcess(itr,log,gamma,num_cont,diff)
    if itr == 0 || itr == 1
        cont = true;
        return
    end
    cont = false;
    for i = 1:num_cont
        if abs(log{max(itr-i,1),1}-gamma) > diff
            cont = true;
        end
    end
end
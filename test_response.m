clc;
clear;
close all;

usePFC = true;
matFile = './results/itr_pfc_2022-03-16T203855.mat'; % Load file that design results saved.
numPlants = 5; % Number of plants to simulate.

data = strcat(matFile);
data = load(data);

% load('./plants/stable_mp.mat')
load('./plants/stable_nmp.mat')

[Ad,Bd,Cd,Dd] = ssdata(data.D);
gamma = data.evalGamma;

rng(4,'twister');
for c = 1:numPlants
    lambda = rand(vertices,1);
    lambda = lambda/sum(lambda,1);
    
    Ap = zeros(size(Am));
    Bp = zeros(size(Bm));
    Cp = zeros(size(Cm));
    for i = 1: vertices
        Ap = Ap + lambda(i)*A(:,:,i);
        Bp = Bp + lambda(i)*B(:,:,i);
        Cp = Cp + lambda(i)*C(:,:,i);
    end
    Bw = Bw;Dw = Dw;
        
    in = size(Bp,2);
    out = size(Cp,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    sspan = [0 40];
    
    t = 0:0.01:40;
    nst = 2000;nen = 4001;
    ist = 1;ien = 4001;

    points = size(t,2);
    wu_size = size(Bw,2);
    wy_size = size(Dw,1);
    u_size = size(Bm,2);
    wu = zeros(points,wu_size);
    wy = zeros(points,wy_size);
    u = zeros(points,u_size);

    wu(nst:nen,:) = 0.5;
    wy(nst:nen,:) = 0.5;

    u(ist:ien,:) = 1;
    
    input_disturbance.time=t';
    input_disturbance.signals.values=wu;
    input_disturbance.signals.dimensions=wu_size;
    observation_noise.time=t';
    observation_noise.signals.values=wy;
    observation_noise.signals.dimenstions=wy_size;
    input_signal.time=t';
    input_signal.signals.values=u;
    input_signal.signals.dimensions=u_size;
    
    if ~usePFC
        out=sim('./simulink_models/MEC_R2016b','StopTime',num2str(40));
    else
        out=sim('./simulink_models/MEC_with_PFC_R2016b','StopTime',num2str(40));
    end
    
    figure(1)
    hold on
    pl(1,:)=plot(out.get('mec').Time,out.get('mec').Data,'linewidth',1,'Color','r');
    pl(2,:)=plot(out.get('no_mec').Time,out.get('no_mec').Data,'linewidth',1,'Color','b');
    pl(3,:)=plot(out.get('ideal').Time,out.get('ideal').Data,'--','linewidth',2,'Color','k');
    xlabel('Time [s]');ylabel('y');
    legend([pl(1),pl(2),pl(3)],{'w/ MEC','w/o MEC','Ideal'},'Location','best')
    hold off
end

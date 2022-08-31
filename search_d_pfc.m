function [gamma,d_ss,xopt] = search_d_pfc(ss_vertices,Bw,Dw,model_ss,lyapunov,pfc,xinit,options)

Am = model_ss.A;
Bm = model_ss.B;
Cm = model_ss.C;
Dm = model_ss.D;

Af = pfc.A;
Bf = pfc.B;
Cf = pfc.C;
Df = pfc.D;

canonical_form = true; % Available only if degree of D is 4.

if canonical_form
    Tal = [0 0;1 0;0 0;0 1];
    Tar = Tal';
    Tas = [0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0];
    Tbr = [0 0;0 1;0 0;0 0];
    Tbs = [0 0;1 0;0 0;0 1];
end

AFm = [Am zeros(size(Am,1),size(Af,2));zeros(size(Af,1),size(Am,2)) Af];
BFm = [Bm;Bf];
CFm = [Cm Cf];

p_state = size(Am,1);
p_in = size(Bm,2);
p_out = size(Cm,1);

pfc_state = size(Af,1);
pf_state = p_state + pfc_state;

d_state = size(lyapunov,1)-2*pf_state;
d_in = p_out;
d_out = p_in;

vertices = size(ss_vertices,1);

X = cell(3,3);
X{1,1} = lyapunov(1:pf_state,1:pf_state);
X{1,2} = lyapunov(1:pf_state,pf_state+1:pf_state+d_state);
X{1,3} = lyapunov(1:pf_state,pf_state+d_state+1:2*pf_state+d_state);
X{2,1} = lyapunov(pf_state+1:pf_state+d_state,1:pf_state);
X{2,2} = lyapunov(pf_state+1:pf_state+d_state,pf_state+1:pf_state+d_state);
X{2,3} = lyapunov(pf_state+1:pf_state+d_state,pf_state+d_state+1:2*pf_state+d_state);
X{3,1} = lyapunov(pf_state+d_state+1:2*pf_state+d_state,1:pf_state);
X{3,2} = lyapunov(pf_state+d_state+1:2*pf_state+d_state,pf_state+1:pf_state+d_state);
X{3,3} = lyapunov(pf_state+d_state+1:2*pf_state+d_state,pf_state+d_state+1:2*pf_state+d_state);

dA = zeros(size(Am,1),size(Am,2),vertices);
dB = zeros(size(Bm,1),size(Bm,2),vertices);
dC = zeros(size(Cm,1),size(Cm,2),vertices);
A = zeros(size(Am,1),size(Am,2),vertices);
B = zeros(size(Bm,1),size(Bm,2),vertices);
C = zeros(size(Cm,1),size(Cm,2),vertices);
AF = zeros(p_state+pfc_state,p_state+pfc_state,vertices);
BF = zeros(p_state+pfc_state,p_in,vertices);
CF = zeros(p_out,p_state+pfc_state,vertices);
dAF = zeros(p_state+pfc_state,p_state+pfc_state,vertices);
dBF = zeros(p_state+pfc_state,p_in,vertices);
dCF = zeros(p_out,p_state+pfc_state,vertices);

BFw = [Bw;zeros(size(Bf))];
DFw = Dw;

for i = 1:vertices
    A(:,:,i) = ss_vertices{i}.A;
    B(:,:,i) = ss_vertices{i}.B;
    C(:,:,i) = ss_vertices{i}.C;
    dA(:,:,i) = A(:,:,i) - Am;
    dB(:,:,i) = B(:,:,i) - Bm;
    dC(:,:,i) = C(:,:,i) - Cm;
    AF(:,:,i) = [A(:,:,i) zeros(size(A,1),size(Af,2)); zeros(size(Af,1),size(A,2)) Af];
    BF(:,:,i) = [B(:,:,i);Bf];
    CF(:,:,i) = [C(:,:,i) Cf];
    dAF(:,:,i) = AF(:,:,i) - AFm;
    dBF(:,:,i) = BF(:,:,i) - BFm;
    dCF(:,:,i) = CF(:,:,i) - CFm;
end

setlmis([])

[gamma2,~,~] = lmivar(1,[1 1]);
if canonical_form
    [Ad,~,sAd] = lmivar(2,[2 2]);
    [Bd,~,sBd] = lmivar(2,[1 1]);
else
    [Ad,~,sAd] = lmivar(2,[d_state d_state]);
    [Bd,~,sBd] = lmivar(2,[d_state d_in]);
end
[Cd,~,sCd] = lmivar(2,[d_out d_state]);
[Dd,n,sDd] = lmivar(2,[d_out d_in]);

S = zeros(vertices,1);

% In the case that design all parameters of D
if ~canonical_form
    for i = 1:vertices
        S(i) = newlmi;
        lmiterm([S(i) 1 1 0],AF(:,:,i)*X{1,1}+X{1,1}*AF(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-BF(:,:,i),CF(:,:,i)*X{1,1},'s')
        lmiterm([S(i) 1 1 Cd],-BF(:,:,i),X{2,1},'s')
        lmiterm([S(i) 1 1 0],dAF(:,:,i)*X{3,1}+X{1,3}*dAF(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-BF(:,:,i),dCF(:,:,i)*X{3,1},'s')

        lmiterm([S(i) 2 1 Bd],1,CF(:,:,i)*X{1,1})
        lmiterm([S(i) 2 1 Ad],1,X{2,1})
        lmiterm([S(i) 2 1 Bd],1,dCF(:,:,i)*X{3,1})
        lmiterm([S(i) 2 1 0],X{2,1}*AF(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,1}*CF(:,:,i)',BF(:,:,i)')
        lmiterm([S(i) 2 1 -Cd],-X{2,2},BF(:,:,i)')
        lmiterm([S(i) 2 1 0],X{2,3}*dAF(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,3}*dCF(:,:,i)',BF(:,:,i)')

        lmiterm([S(i) 2 2 Bd],1,CF(:,:,i)*X{1,2},'s')
        lmiterm([S(i) 2 2 Ad],1,X{2,2},'s')
        lmiterm([S(i) 2 2 Bd],1,dCF(:,:,i)*X{3,2},'s')

        lmiterm([S(i) 3 1 0],AFm*X{3,1})
        lmiterm([S(i) 3 1 0],X{3,1}*AF(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,1}*CF(:,:,i)',BF(:,:,i)')
        lmiterm([S(i) 3 1 -Cd],-X{3,2},BF(:,:,i)')
        lmiterm([S(i) 3 1 0],X{3,3}*dAF(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,3}*dCF(:,:,i)',BF(:,:,i)')

        lmiterm([S(i) 3 2 0],AFm*X{3,2})
        lmiterm([S(i) 3 2 -Bd],X{3,1}*CF(:,:,i)',1)
        lmiterm([S(i) 3 2 -Ad],X{3,2},1)
        lmiterm([S(i) 3 2 -Bd],X{3,3}*dCF(:,:,i)',1)

        lmiterm([S(i) 3 3 0],AFm*X{3,3}+X{3,3}*AFm')

    %     C_F*x_F - C_Fm*x_Fm
    %     lmiterm([S(i) 4 1 0],CF(:,:,i)*X{1,1}+dCF(:,:,i)*X{3,1})
    %     lmiterm([S(i) 4 2 0],CF(:,:,i)*X{1,2}+dCF(:,:,i)*X{3,2})
    %     lmiterm([S(i) 4 3 0],CF(:,:,i)*X{1,3}+dCF(:,:,i)*X{3,3})

    %     C*x - C_m*x_m
        lmiterm([S(i) 4 1 0],[C(:,:,i) zeros(p_out,pfc_state)]*X{1,1}+[dC(:,:,i) zeros(p_out,pfc_state)]*X{3,1})
        lmiterm([S(i) 4 2 0],[C(:,:,i) zeros(p_out,pfc_state)]*X{1,2}+[dC(:,:,i) zeros(p_out,pfc_state)]*X{3,2})
        lmiterm([S(i) 4 3 0],[C(:,:,i) zeros(p_out,pfc_state)]*X{1,3}+[dC(:,:,i) zeros(p_out,pfc_state)]*X{3,3})

        lmiterm([S(i) 5 1 0],BFw')
        lmiterm([S(i) 6 1 -Dd],-DFw',BF(:,:,i)')
        lmiterm([S(i) 6 2 -Bd],DFw',1)
        lmiterm([S(i) 7 1 0],dBF(:,:,i)')
        lmiterm([S(i) 7 3 0],BFm')

        lmiterm([S(i) 4 4 gamma2],-eye(2),1)

        lmiterm([S(i) 5 5 0],-eye(2))
        lmiterm([S(i) 6 6 0],-eye(2))
        lmiterm([S(i) 7 7 0],-eye(2))
    end
else
    % Design D as control canonical form
    for i = 1:vertices
        S(i) = newlmi;
        lmiterm([S(i) 1 1 0],AF(:,:,i)*X{1,1}+X{1,1}*AF(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-BF(:,:,i),CF(:,:,i)*X{1,1},'s')
        lmiterm([S(i) 1 1 Cd],-BF(:,:,i),X{2,1},'s')
        lmiterm([S(i) 1 1 0],dAF(:,:,i)*X{3,1}+X{1,3}*dAF(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-BF(:,:,i),dCF(:,:,i)*X{3,1},'s')

        lmiterm([S(i) 2 1 Bd],1,Tbr*CF(:,:,i)*X{1,1})
        lmiterm([S(i) 2 1 Ad],Tal,Tar*X{2,1})
        lmiterm([S(i) 2 1 Bd],1,Tbr*dCF(:,:,i)*X{3,1})
        lmiterm([S(i) 2 1 0],Tbs*CF(:,:,i)*X{1,1})
        lmiterm([S(i) 2 1 0],Tas*X{2,1})
        lmiterm([S(i) 2 1 0],Tbs*dCF(:,:,i)*X{3,1})
        lmiterm([S(i) 2 1 0],X{2,1}*AF(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,1}*CF(:,:,i)',BF(:,:,i)')
        lmiterm([S(i) 2 1 -Cd],-X{2,2},BF(:,:,i)')
        lmiterm([S(i) 2 1 0],X{2,3}*dAF(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,3}*dCF(:,:,i)',BF(:,:,i)')

        lmiterm([S(i) 2 2 Bd],1,Tbr*CF(:,:,i)*X{1,2},'s')
        lmiterm([S(i) 2 2 Ad],Tal,Tar*X{2,2},'s')
        lmiterm([S(i) 2 2 Bd],1,Tbr*dCF(:,:,i)*X{3,2},'s')
        lmiterm([S(i) 2 2 0],Tbs*CF(:,:,i)*X{1,2}+(Tbs*CF(:,:,i)*X{1,2})')
        lmiterm([S(i) 2 2 0],Tas*X{2,2}+(Tas*X{2,2})')
        lmiterm([S(i) 2 2 0],Tbs*dCF(:,:,i)*X{3,2}+(Tbs*dCF(:,:,i)*X{3,2})')

        lmiterm([S(i) 3 1 0],AFm*X{3,1})
        lmiterm([S(i) 3 1 0],X{3,1}*AF(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,1}*CF(:,:,i)',BF(:,:,i)')
        lmiterm([S(i) 3 1 -Cd],-X{3,2},BF(:,:,i)')
        lmiterm([S(i) 3 1 0],X{3,3}*dAF(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,3}*dCF(:,:,i)',BF(:,:,i)')

        lmiterm([S(i) 3 2 0],AFm*X{3,2})
        lmiterm([S(i) 3 2 -Bd],X{3,1}*CF(:,:,i)'*Tbr',1)
        lmiterm([S(i) 3 2 -Ad],X{3,2}*Tar',Tal')
        lmiterm([S(i) 3 2 -Bd],X{3,3}*dCF(:,:,i)'*Tbr',1)
        lmiterm([S(i) 3 2 0],X{3,1}*CF(:,:,i)'*Tbs')
        lmiterm([S(i) 3 2 0],X{3,2}*Tas')
        lmiterm([S(i) 3 2 0],X{3,3}*dCF(:,:,i)'*Tbs')

        lmiterm([S(i) 3 3 0],AFm*X{3,3}+X{3,3}*AFm')

    %     C_F*x_F - C_Fm*x_Fm
    %     lmiterm([S(i) 4 1 0],CF(:,:,i)*X{1,1}+dCF(:,:,i)*X{3,1})
    %     lmiterm([S(i) 4 2 0],CF(:,:,i)*X{1,2}+dCF(:,:,i)*X{3,2})
    %     lmiterm([S(i) 4 3 0],CF(:,:,i)*X{1,3}+dCF(:,:,i)*X{3,3})

    %     C*x - C_m*x_m
        lmiterm([S(i) 4 1 0],[C(:,:,i) zeros(p_out,pfc_state)]*X{1,1}+[dC(:,:,i) zeros(p_out,pfc_state)]*X{3,1})
        lmiterm([S(i) 4 2 0],[C(:,:,i) zeros(p_out,pfc_state)]*X{1,2}+[dC(:,:,i) zeros(p_out,pfc_state)]*X{3,2})
        lmiterm([S(i) 4 3 0],[C(:,:,i) zeros(p_out,pfc_state)]*X{1,3}+[dC(:,:,i) zeros(p_out,pfc_state)]*X{3,3})

        lmiterm([S(i) 5 1 0],BFw')
        lmiterm([S(i) 6 1 -Dd],-DFw',BF(:,:,i)')
        lmiterm([S(i) 6 2 -Bd],DFw'*Tbr',1)
        lmiterm([S(i) 6 2 0],DFw'*Tbs')
        lmiterm([S(i) 7 1 0],dBF(:,:,i)')
        lmiterm([S(i) 7 3 0],BFm')

        lmiterm([S(i) 4 4 gamma2],-eye(2),1)

        lmiterm([S(i) 5 5 0],-eye(2))
        lmiterm([S(i) 6 6 0],-eye(2))
        lmiterm([S(i) 7 7 0],-eye(2))
    end
end

LMIs = getlmis;
c = [1 zeros(1,n-1)];

if size(xinit,1) ~= n
    xinit = zeros(n,1);
end

if exist('xinit','var')
    [copt,xopt] = mincx(LMIs,c,options,xinit);
else
    [copt,xopt] = mincx(LMIs,c,options);
end

if size(xopt,1)~=0
    gamma2opt = dec2mat(LMIs,xopt,gamma2);
else
    d_ss = NaN;
    gamma = -1;
    xopt = NaN;
    return
end

gamma2opt = dec2mat(LMIs,xopt,gamma2);
Ad = dec2mat(LMIs,xopt,Ad);
Bd = dec2mat(LMIs,xopt,Bd);
Cd = dec2mat(LMIs,xopt,Cd);
Dd = dec2mat(LMIs,xopt,Dd);

if canonical_form
    Ad = Tal*Ad*Tar+Tas;
    Bd = Bd*Tbr+Tbs;
end

gamma = sqrt(gamma2opt);
d_ss = ss(Ad,Bd,Cd,Dd);

end
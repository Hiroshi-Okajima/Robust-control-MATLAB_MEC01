function [gamma,d_ss,d_opt] = search_d(ss_vertices,Bw,Dw,model_ss,lyapunov,xopt,options)

Am = model_ss.A;
Bm = model_ss.B;
Cm = model_ss.C;
Dm = model_ss.D;

p_state = size(Am,1);
p_in = size(Bm,2);
p_out = size(Cm,1);

canonical_form = true; % Available only if degree of D is 4.

if canonical_form
    Tal = [0 0;1 0;0 0;0 1];
    Tar = Tal';
    Tas = [0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0];
    Tbr = [0 0;0 1;0 0;0 0];
    Tbs = [0 0;1 0;0 0;0 1];
end

d_state = size(lyapunov,1)-2*p_state;
d_in = p_out;
d_out = p_in;

vertices = size(ss_vertices,1);

X = cell(3,3);
X{1,1} = lyapunov(1:p_state,1:p_state);
X{1,2} = lyapunov(1:p_state,p_state+1:p_state+d_state);
X{1,3} = lyapunov(1:p_state,p_state+d_state+1:2*p_state+d_state);
X{2,1} = lyapunov(p_state+1:p_state+d_state,1:p_state);
X{2,2} = lyapunov(p_state+1:p_state+d_state,p_state+1:p_state+d_state);
X{2,3} = lyapunov(p_state+1:p_state+d_state,p_state+d_state+1:2*p_state+d_state);
X{3,1} = lyapunov(p_state+d_state+1:2*p_state+d_state,1:p_state);
X{3,2} = lyapunov(p_state+d_state+1:2*p_state+d_state,p_state+1:p_state+d_state);
X{3,3} = lyapunov(p_state+d_state+1:2*p_state+d_state,p_state+d_state+1:2*p_state+d_state);

dA = zeros(size(Am,1),size(Am,2),vertices);
dB = zeros(size(Bm,1),size(Bm,2),vertices);
dC = zeros(size(Cm,1),size(Cm,2),vertices);
A = zeros(size(Am,1),size(Am,2),vertices);
B = zeros(size(Bm,1),size(Bm,2),vertices);
C = zeros(size(Cm,1),size(Cm,2),vertices);

for i = 1:vertices
    A(:,:,i) = ss_vertices{i}.A;
    B(:,:,i) = ss_vertices{i}.B;
    C(:,:,i) = ss_vertices{i}.C;
    dA(:,:,i) = ss_vertices{i}.A - Am;
    dB(:,:,i) = ss_vertices{i}.B - Bm;
    dC(:,:,i) = ss_vertices{i}.C - Cm;
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
        lmiterm([S(i) 1 1 0],A(:,:,i)*X{1,1}+X{1,1}*A(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-B(:,:,i),C(:,:,i)*X{1,1},'s')
        lmiterm([S(i) 1 1 Cd],-B(:,:,i),X{2,1},'s')
        lmiterm([S(i) 1 1 0],dA(:,:,i)*X{3,1}+X{1,3}*dA(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-B(:,:,i),dC(:,:,i)*X{3,1},'s')

        lmiterm([S(i) 2 1 Bd],1,C(:,:,i)*X{1,1})
        lmiterm([S(i) 2 1 Ad],1,X{2,1})
        lmiterm([S(i) 2 1 Bd],1,dC(:,:,i)*X{3,1})
        lmiterm([S(i) 2 1 0],X{2,1}*A(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,1}*C(:,:,i)',B(:,:,i)')
        lmiterm([S(i) 2 1 -Cd],-X{2,2},B(:,:,i)')
        lmiterm([S(i) 2 1 0],X{2,3}*dA(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,3}*dC(:,:,i)',B(:,:,i)')

        lmiterm([S(i) 2 2 Bd],1,C(:,:,i)*X{1,2},'s')
        lmiterm([S(i) 2 2 Ad],1,X{2,2},'s')
        lmiterm([S(i) 2 2 Bd],1,dC(:,:,i)*X{3,2},'s')

        lmiterm([S(i) 3 1 0],Am*X{3,1})
        lmiterm([S(i) 3 1 0],X{3,1}*A(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,1}*C(:,:,i)',B(:,:,i)')
        lmiterm([S(i) 3 1 -Cd],-X{3,2},B(:,:,i)')
        lmiterm([S(i) 3 1 0],X{3,3}*dA(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,3}*dC(:,:,i)',B(:,:,i)')

        lmiterm([S(i) 3 2 0],Am*X{3,2})
        lmiterm([S(i) 3 2 -Bd],X{3,1}*C(:,:,i)',1)
        lmiterm([S(i) 3 2 -Ad],X{3,2},1)
        lmiterm([S(i) 3 2 -Bd],X{3,3}*dC(:,:,i)',1)

        lmiterm([S(i) 3 3 0],Am*X{3,3}+X{3,3}*Am')

        lmiterm([S(i) 4 1 0],C(:,:,i)*X{1,1}+dC(:,:,i)*X{3,1})
        lmiterm([S(i) 4 2 0],C(:,:,i)*X{1,2}+dC(:,:,i)*X{3,2})
        lmiterm([S(i) 4 3 0],C(:,:,i)*X{1,3}+dC(:,:,i)*X{3,3})

        lmiterm([S(i) 5 1 0],Bw')
        lmiterm([S(i) 6 1 -Dd],-Dw',B(:,:,i)')
        lmiterm([S(i) 6 2 -Bd],Dw',1)
        lmiterm([S(i) 7 1 0],dB(:,:,i)')
        lmiterm([S(i) 7 3 0],Bm')

        lmiterm([S(i) 4 4 gamma2],-eye(2),1)

        lmiterm([S(i) 5 5 0],-eye(2))
        lmiterm([S(i) 6 6 0],-eye(2))
        lmiterm([S(i) 7 7 0],-eye(2))
    end
else
    for i = 1:vertices
        S(i) = newlmi;
        lmiterm([S(i) 1 1 0],A(:,:,i)*X{1,1}+X{1,1}*A(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-B(:,:,i),C(:,:,i)*X{1,1},'s')
        lmiterm([S(i) 1 1 Cd],-B(:,:,i),X{2,1},'s')
        lmiterm([S(i) 1 1 0],dA(:,:,i)*X{3,1}+X{1,3}*dA(:,:,i)')
        lmiterm([S(i) 1 1 Dd],-B(:,:,i),dC(:,:,i)*X{3,1},'s')

        lmiterm([S(i) 2 1 Bd],1,Tbr*C(:,:,i)*X{1,1})
        lmiterm([S(i) 2 1 Ad],Tal,Tar*X{2,1})
        lmiterm([S(i) 2 1 Bd],1,Tbr*dC(:,:,i)*X{3,1})
        lmiterm([S(i) 2 1 0],Tbs*C(:,:,i)*X{1,1})
        lmiterm([S(i) 2 1 0],Tas*X{2,1})
        lmiterm([S(i) 2 1 0],Tbs*dC(:,:,i)*X{3,1})
        lmiterm([S(i) 2 1 0],X{2,1}*A(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,1}*C(:,:,i)',B(:,:,i)')
        lmiterm([S(i) 2 1 -Cd],-X{2,2},B(:,:,i)')
        lmiterm([S(i) 2 1 0],X{2,3}*dA(:,:,i)')
        lmiterm([S(i) 2 1 -Dd],-X{2,3}*dC(:,:,i)',B(:,:,i)')

        lmiterm([S(i) 2 2 Bd],1,Tbr*C(:,:,i)*X{1,2},'s')
        lmiterm([S(i) 2 2 Ad],Tal,Tar*X{2,2},'s')
        lmiterm([S(i) 2 2 Bd],1,Tbr*dC(:,:,i)*X{3,2},'s')
        lmiterm([S(i) 2 2 0],Tbs*C(:,:,i)*X{1,2}+(Tbs*C(:,:,i)*X{1,2})')
        lmiterm([S(i) 2 2 0],Tas*X{2,2}+(Tas*X{2,2})')
        lmiterm([S(i) 2 2 0],Tbs*dC(:,:,i)*X{3,2}+(Tbs*dC(:,:,i)*X{3,2})')

        lmiterm([S(i) 3 1 0],Am*X{3,1})
        lmiterm([S(i) 3 1 0],X{3,1}*A(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,1}*C(:,:,i)',B(:,:,i)')
        lmiterm([S(i) 3 1 -Cd],-X{3,2},B(:,:,i)')
        lmiterm([S(i) 3 1 0],X{3,3}*dA(:,:,i)')
        lmiterm([S(i) 3 1 -Dd],-X{3,3}*dC(:,:,i)',B(:,:,i)')

        lmiterm([S(i) 3 2 0],Am*X{3,2})
        lmiterm([S(i) 3 2 -Bd],X{3,1}*C(:,:,i)'*Tbr',1)
        lmiterm([S(i) 3 2 -Ad],X{3,2}*Tar',Tal')
        lmiterm([S(i) 3 2 -Bd],X{3,3}*dC(:,:,i)'*Tbr',1)
        lmiterm([S(i) 3 2 0],X{3,1}*C(:,:,i)'*Tbs')
        lmiterm([S(i) 3 2 0],X{3,2}*Tas')
        lmiterm([S(i) 3 2 0],X{3,3}*dC(:,:,i)'*Tbs')

        lmiterm([S(i) 3 3 0],Am*X{3,3}+X{3,3}*Am')

        lmiterm([S(i) 4 1 0],C(:,:,i)*X{1,1}+dC(:,:,i)*X{3,1})
        lmiterm([S(i) 4 2 0],C(:,:,i)*X{1,2}+dC(:,:,i)*X{3,2})
        lmiterm([S(i) 4 3 0],C(:,:,i)*X{1,3}+dC(:,:,i)*X{3,3})

        lmiterm([S(i) 5 1 0],Bw')
        lmiterm([S(i) 6 1 -Dd],-Dw',B(:,:,i)')
        lmiterm([S(i) 6 2 -Bd],Dw'*Tbr',1)
        lmiterm([S(i) 6 2 0],Dw'*Tbs')
        lmiterm([S(i) 7 1 0],dB(:,:,i)')
        lmiterm([S(i) 7 3 0],Bm')

        lmiterm([S(i) 4 4 gamma2],-eye(2),1)

        lmiterm([S(i) 5 5 0],-eye(2))
        lmiterm([S(i) 6 6 0],-eye(2))
        lmiterm([S(i) 7 7 0],-eye(2))
    end
end

LMIs = getlmis;
c = [1 zeros(1,n-1)];

xinit = xopt;
if size(xopt,1) ~= n
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
    d_opt = NaN;
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

d_opt = xopt;
gamma = sqrt(gamma2opt);
d_ss = ss(Ad,Bd,Cd,Dd);

end
function [gamma,lyapunov,xopt] = search_x(ss_vertices,Bw,Dw,model_ss,d_ss,xinit,options)

Ad = d_ss.A;
Bd = d_ss.B;
Cd = d_ss.C;
Dd = d_ss.D;

Am = model_ss.A;
Bm = model_ss.B;
Cm = model_ss.C;
Dm = model_ss.D;

d_state = size(Ad,1);
d_in = size(Bd,2);
d_out = size(Cd,1);

p_state = size(Am,1);
p_in = size(Bm,2);
p_out = size(Cm,1);

vertices = size(ss_vertices,1);

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

sizes = [{[p_state,p_state]},{[p_state,d_state]},{[p_state,p_state]};
         {[d_state,p_state]},{[d_state,d_state]},{[d_state,p_state]};
         {[p_state,p_state]},{[p_state,d_state]},{[p_state,p_state]}];

X = zeros(3);
for i=1:3
    for j = 1:i
        if j==i
            [X(i,j),n,sX{i,j}] = lmivar(1,[sizes{i,j}(1),1]);
        else
            [X(i,j),n,sX{i,j}] = lmivar(2,[sizes{i,j}(1),sizes{i,j}(2)]);
        end
    end
end

S = zeros(vertices,1);
for i = 1:vertices
    S(i) = newlmi;
    % LMI(1,1)
    lmiterm([S(i) 1 1 X(1,1)],A(:,:,i)-B(:,:,i)*Dd*C(:,:,i),1,'s')
    lmiterm([S(i) 1 1 X(2,1)],-B(:,:,i)*Cd,1,'s')
    lmiterm([S(i) 1 1 X(3,1)],dA(:,:,i)-B(:,:,i)*Dd*dC(:,:,i),1,'s')
    lmiterm([S(i) 2 1 X(1,1)],Bd*C(:,:,i),1)
    lmiterm([S(i) 2 1 X(2,1)],Ad,1)
    lmiterm([S(i) 2 1 X(3,1)],Bd*dC(:,:,i),1)
    lmiterm([S(i) 2 1 X(2,1)],1,A(:,:,i)'-C(:,:,i)'*Dd'*B(:,:,i)')
    lmiterm([S(i) 2 1 X(2,2)],-1,Cd'*B(:,:,i)')
    lmiterm([S(i) 2 1 -X(3,2)],1,dA(:,:,i)'-dC(:,:,i)'*Dd'*B(:,:,i)')
    lmiterm([S(i) 3 1 X(3,1)],Am,1)
    lmiterm([S(i) 3 1 X(3,1)],1,A(:,:,i)'-C(:,:,i)'*Dd'*B(:,:,i)')
    lmiterm([S(i) 3 1 X(3,2)],-1,Cd'*B(:,:,i)')
    lmiterm([S(i) 3 1 X(3,3)],1,dA(:,:,i)'-dC(:,:,i)'*Dd'*B(:,:,i)')
    lmiterm([S(i) 2 2 -X(2,1)],Bd*C(:,:,i),1,'s')
    lmiterm([S(i) 2 2 X(2,2)],Ad,1,'s')
    lmiterm([S(i) 2 2 X(3,2)],Bd*dC(:,:,i),1,'s')
    lmiterm([S(i) 3 2 X(3,2)],Am,1)
    lmiterm([S(i) 3 2 X(3,1)],1,C(:,:,i)'*Bd')
    lmiterm([S(i) 3 2 X(3,2)],1,Ad')
    lmiterm([S(i) 3 2 X(3,3)],1,dC(:,:,i)'*Bd')
    lmiterm([S(i) 3 3 X(3,3)],Am,1,'s')

    % LMI(2,1)
    lmiterm([S(i) 4 1 X(1,1)],C(:,:,i),1)
    lmiterm([S(i) 4 1 X(3,1)],dC(:,:,i),1)
    lmiterm([S(i) 4 2 -X(2,1)],C(:,:,i),1)
    lmiterm([S(i) 4 2 X(3,2)],dC(:,:,i),1)
    lmiterm([S(i) 4 3 -X(3,1)],C(:,:,i),1)
    lmiterm([S(i) 4 3 X(3,3)],dC(:,:,i),1)

    % LMI(3,1)
    lmiterm([S(i) 5 1 0],Bw')
    lmiterm([S(i) 6 1 0],-Dw'*Dd'*B(:,:,i)')
    lmiterm([S(i) 7 1 0],dB(:,:,i)')
    lmiterm([S(i) 6 2 0],Dw'*Bd')
    lmiterm([S(i) 7 3 0],Bm')

    % LMI(2,2)
    lmiterm([S(i) 4 4 gamma2],-eye(2),1)

    % LMI(3,3)
    lmiterm([S(i) 5 5 0],-eye(2))
    lmiterm([S(i) 6 6 0],-eye(2))
    lmiterm([S(i) 7 7 0],-eye(2))
end

Xlmi = newlmi;
for i = 1:3
    for j = 1:i
        lmiterm([-Xlmi i j X(i,j)],1,1)
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
    lyapunov = [];
    gamma = -1;
    xopt = NaN;
    return
end

lyapunov_opt = cell(3,3);
for i = 1:3
    for j = 1:i
        lyapunov_opt{i,j} = dec2mat(LMIs,xopt,X(i,j));
    end
end

lyapunov = [lyapunov_opt{1,1} lyapunov_opt{2,1}' lyapunov_opt{3,1}';lyapunov_opt{2,1} lyapunov_opt{2,2} lyapunov_opt{3,2}';lyapunov_opt{3,1} lyapunov_opt{3,2} lyapunov_opt{3,3}];
gamma = sqrt(gamma2opt);

end
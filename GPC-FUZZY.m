clear, clc

Pot = 11; %Pot=3-11W;
Rco = 1.8*1.5e-1;
L = 8.64e-3; %Indcutor H
Co = 66.7e-6; %Capacitor F
Vo = 20; %Output voltage
Vg =15; %Input voltage Vg=12-15V;
Ro = (Vo^2)/Pot; %input load
Dcycle = 1 - Vg/Vo; %duty cicle
R = ((1 - Dcycle)^2)*Ro + Dcycle*(1 - Dcycle)*((Rco*Ro)/(Rco+Ro)); %
[A1, B1 , C1, D1] = BOOST(Vg, L, Dcycle, Rco, Ro, Co);

%% Discrete Model Transfer Function
Ts = 1e-3; %sample time
[Bc, Ac] = ss2tf(A1, B1, C1, D1);
Tz = c2d(tf(Bc, Ac), Ts) % Discretized
A = cell2mat(Tz.den); %   A(z-¹)
B = cell2mat(Tz.num); %   B(z-¹)
time=0:Ts:25e-3;
N = length(time)-1;     %   iteration
N2=8;                 %   prediction horizon
Nu=7;                 %   input horizon
lambda = 0.8;
%% constraints
Dumax=0.2;             %input rate limit
umax=0.75;   umin=0;    %max input
ymax=20;     ymin=-0.4;    %output limit
ref(1:N) =20;          %reference

%% Fuzzy Variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rage=[-Vo:1:Vo];
Ne = trapmf(rage, [-Vo -Vo -15 0]); %Negative set
Ze = trimf(rage, [-15 0 15]); %Zero set
Pe = trapmf(rage, [0 15 Vo Vo]); %Positive set
figure(4);
plot(rage, Ne); hold on; plot(rage,Ze); hold on; plot(rage,Pe);
% Ne =  gaussmf(rage,  [8.493 -20]);
% Ze= gaussmf(rage,   [8.493 -2.22e-16]);
% Pe= gaussmf(rage,   [8.493 20]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Find prediction matrices
% yfut = H *Dufut + P*Dupast + Q*ypast
D = [1,-1]; %delta
sizey = size(A,1);
if size(B,2)==sizey;B=[B,zeros(sizey,sizey)];end
AD = conv(A,D);
nA = length(AD);
nB = length(B);
%%%%  Initialise Recursion data
%%% nominal model    y =  Bo ut + B2 Dupast  + A2 ypast
A2 = -AD(1,1+1:nA);
B2 = B(1,1+1:nB);
Bo = B(1,1:1);
nB2 = nB-1;
nA2 = nA-1;
P1=Bo; P2=B2; P3=A2;
%%%%% Loop updating models using recursion
for i=2:N2;
    vecold = (i-2)*1+1:(i-1)*1;
    vecnew = (i-1)*1+1:i*1;
    Phi = P3(vecold,1:1);
    vecufut = 1:1*i;
    vecufut2 = 1:1*(i-1);
    P1(vecnew,vecufut) = [(P2(vecold,1:1)+Phi*Bo),P1(vecold,vecufut2)];
    vecupast = 1+1:nB2;
    vecypast = 1+1:nA2;
    temp = [P2(vecold,vecupast),zeros(1,1)] + Phi*B2;
    P2(vecnew,1:nB2) = temp;
    P3(vecnew,1:nA2) = [P3(vecold,vecypast),zeros(1,1)] + Phi*A2;
end
H=P1; P=P2; Q=P3;
%%cost function
%Dufut = Pr*rfut - Dk*Dupast - Nk*ypast
%J = Dufut'*S*Dufut + Dufut'*2X*[Dupast;ypast;rfut]
% Control horizon
P1 = H(:,1:1*Nu);
L = ones(N2, 1);
% Define performance index parameters
S = P1'*P1 + lambda*eye(Nu);
X = [P1'*P,P1'*Q,-P1'];
% Define the control law parameters
Nk = inv(S)*P1'*Q;
Dk = inv(S)*P1'*P;
% Remove advance knowledge on the set point
Pr = inv(S)*P1'*L;
X = [P1'*P,P1'*Q,-P1'*L];

%% Define constraint matrices
%  CC*Du(future) - dd - du*u(k-1))-ddu*Dupast-dy*ypast <= 0
d(1:2*Nu, 1) = ones(2*Nu,1)*Dumax;
C(1:Nu,1:Nu) = eye(Nu);
C(Nu+1:2*Nu,1:Nu) = -eye(Nu);
d(2*Nu+1:4*Nu, 1) = [ones(Nu,1)*umax;-ones(Nu,1)*umin];
d(4*Nu+1:4*Nu+2*N2,1) = [ones(N2,1)*ymax;-ones(N2,1)*ymin];
C(2*Nu+1:3*Nu,1:Nu) = tril(ones(Nu,Nu));
C(3*Nu+1:4*Nu,1:Nu) = -tril(ones(Nu,Nu));
d1(2*Nu+1:4*Nu,1) = [-ones(Nu,1);ones(Nu,1)];
% Combine loop constraints into single inequality condition
CC(1:4*Nu,1:1:Nu) = C;
dd(1:(4*Nu+2*N2),1) = d(1:4*Nu+2*N2,1);
du(1:4*Nu,1) = d1(1:4*Nu,1);
% Add output constraints details
CC=[CC;H(:,1:Nu*1);-H(:,1:Nu*1)];
ddu=[zeros(Nu*4*1,size(P,2));-P;P];
dy=[zeros(Nu*4*1,size(Q,2));-Q;Q];
du=[du;zeros(2*N2*1,1)];

%% Set up simulation parameters
nNk = size(Nk,2)/1;
nDk = size(Dk,2)/1;
init = max([nNk,nDk])+2;
y = zeros(1,init); %output
u = y; ufast=u; Du = u;
r = u; d=u;

opt = optimset('quadprog');
opt.Diagnostics='off';    %%%%% Switches of unwanted MATLAB displays
opt.LargeScale='off';     %%%%% However no warning of infeasibility
opt.Display='off';
opt.Algorithm='active-set';
u1=u; ue=u;
y1=y;
load fuzyda;
ee=0;
%% Closed-loop simulation
for i=init:N-1;  
    %% Update unconstrained control law
    d(1:1,i+1)= 0;
    % d(1:1,i+1) = awgn(sawtooth(y(:,i)),5,'measured')
    ypast = y(:, i:-1:i+1-nNk)
    Dupast = Du(:, i-1:-1:i-nDk) ;
    upast = u(:, i-1);
    rfut = ref(:,i+1);
    %% Unconstrained law - for overlay plots
    Dufast(:,i) = Pr*rfut - Nk*ypast(:) - Dk*Dupast(:);
    ufast(:,i)=u(:,i-1)+Dufast(1:1,i);
    %% Form constraint matrices and solve constrained optimisation
    %  CC*Du(future) - dd - du*u(k-1))-ddu*Dupast-dy*ypast <= 0
    %  J = Dufut'S*Dufut+2*X*[Dupast(:);ypast(:);rfut(:)]*Dufut
    dt = dd+du*upast(:)+ddu*Dupast(:)+dy*ypast(:);
    f=X*[Dupast(:);ypast(:);rfut(:)];
    [Duqp,fval,exitflag] = quadprog(S,f,CC,dt,[],[],[],[],[],opt); %Quadratic Program
    
    %% Fuzzy
    e(i) = ref(:,i)-y(:,i);
    de = ee -e(i);
    ee=e(i);
    DOF=interp1(rage, Pe, ee); %DOF of Rule 1
    DU1=min(Pe,DOF);
    DOF=interp1(rage,Ze,ee); %DOF of Rule 2
    DU2=min(Ze,DOF);
    DOF=interp1(rage,Ne,ee); %DOF of Rule 3
    DU3=min(Ne,DOF);
    DU=max(max(DU1,DU2), DU3);
    du1(i)=0.1*defuzz(rage, DU, 'bisector');
    if(du1(i)>Dumax) du1(i)=Dumax; end
    u1(:,i)=u1(:,i-1)+du1(i); % signal control   
     due(i) = (evalfis([de*4 e(i)*4], fuzyda))/4; 
      ue(:,i) = ue(:,i-1)+due(i);
  
    %% GPC
    Du(:,i) = Duqp(1:1);
       med(i) = (du1(i)+Du(:,i))/2;
    u(:,i) = u(:,i-1)+Du(:,i);
     
    % Ensure the constraints satisfied by proposed control law
    for j=1:1;
        if u(j,i)>u(j,i-1)+Dumax(j);u(j,i)=u(j,i-1)+Dumax(j);end
        if u(j,i)<u(j,i-1)-Dumax(j);u(j,i)=u(j,i-1)-Dumax(j);end
        if u(j,i)>umax(j); u(i)=umax(j);end
        if u(j,i)<umin(j); u(i)=umin(j);end
    end
    for j=1:1;
        if u1(j,i)>u1(j,i-1)+Dumax(j);u1(j,i)=u1(j,i-1)+Dumax(j);end
        if u1(j,i)<u1(j,i-1)-Dumax(j);u1(j,i)=u1(j,i-1)-Dumax(j);end
        if u1(j,i)>umax(j); u1(i)=umax(j);end
        if u1(j,i)<umin(j); u1(i)=umin(j);end
    end
   
    Du(:,i) = u(:,i)-u(:,i-1);
    % End of update to the control law

    % Simulate the process
    upast2 = u(:,i:-1:i-nDk);
    upast21 = u1(:,i:-1:i-nDk); %Fuzzy
    ypast2 = y(:, i:-1:i+2-nNk);
    ypast21 = y1(:, i:-1:i+2-nNk);
    y(:,i+1) = -A(:,1+1:nNk*1)*ypast2(:) + B*[upast2(:)] + d(:,i+1);
    y1(:,i+1) = -A(:,1+1:nNk*1)*ypast21(:) + B*[upast21(:)] + d(:,i+1); %Fuzzy 
    r(:,i+1) = ref(:,i+1);  
end
%% Plot
figure(1); clf reset
upper_u(1:N) =umax;
subplot(211); plot(time(1:end-5),y(5:end), 'g', 'LineWidth', 2); hold on; plot(time(1:end-5),ref(5:end), '--k'); hold on; plot(time(1:end-5),y1(5:end), 'r', 'LineWidth', 2); ylabel('Ouput')
%axis([1 20 0 ymax+5])
subplot(212); plot(u(5:end), 'LineWidth', 2); hold on; plot(upper_u(5:end), '--k'); hold on; plot(u1(5:end),'r' ,'LineWidth', 2); title('u(t)'); legend('GPC', 'umax' ,'GPC-Fuzzy')



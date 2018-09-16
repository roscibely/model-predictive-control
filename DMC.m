%Dynamic matrix control applied boost convert
clear, close, clc
R=8; %Predition Horizon
L=4; %Control Horizon
alpha =0.6; %Parameter determines how fast the trajectory reaches the setpoint
T=1e-3; %sample time
t=0:T:20e-3;
N = length(t); %Horizon model
umax =0.5; 
%First phase Vg=36V Pot=1000W 
A1=[-0.3004 -7.7425; 0.0616 -0.1297];
B1=[581.3690 66.4196]';
C1=[0.0198 0.9986]; D1=-0.7290;

%Open-loop
x=[0 0]'; u=1; y=0;
for k=1:N
    x(:,k+1)=A1*x(:,k)+B1*u;
    y(k)=C1*x(:,k)+D1*u;
end
figure(1); plot(t,y); title('Open-loop response'); 
figure(2); subplot(211); plot(x(1,:)); title('Current'); subplot(212); plot(x(2,:)); title('capacitor voltage'); 
%Dynamic matrix control 
a=y(2:end); %step
h(1)=a(1);
for j=2:length(a);
    h(j) = a(j) - a(j-1);
end
f=0;  
A = toeplitz(a(1:R), [a(1) zeros(1,L-1)]);
disp('Dynamic Matrix A'); disp(A);
Kc = inv(A'*A + f*eye(L))*A';
KT = Kc(1,:); %gain
disp('Matrix Kc '); disp(Kc);
yr(1)=0;
u=0; 
for k=1:1+40
    time(k) = (k-1)*T;
    r(k) =48; % set point
    for m=1:R 
        S(m) =0;
        for i=m+1:N-2 
            if k+m-i>0
                S(m)=S(m) +h(i)*deltau(k+m-i);
            end
        end
    end
    for i=1:R
        P(i)=0;
        for m=1:i
            P(i)=P(i)+S(m);
        end
    end
    E(k) = r(k)-yr(k);
    for i=1:R
        El(i)= (1-alpha^i)*E(k)-P(i);
    end
    deltau(k) = KT*El';
    if k==1 
        u(k) = deltau(k);
    else 
        u(k) = u(k-1) +deltau(k);
    end
    %state space system 
    x(:,k+1)=A1*x(:,k)+B1*u(k);
    yr(k+1)=C1*x(:,k)+D1*u(k);   
end
figure; subplot(211); plot(yr, 'linewidth', 2); hold on; plot(r,'k--'); hold on;
plot(x(2,:),'r--','linewidth',2); legend('output', 'set point' ,'Capacitor voltage');
title('output'); ylabel('y'); subplot(212); plot(u,'linewidth',2); title('control signal')

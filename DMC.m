clear all
close all
clc

% Definição do sistema

m=2;                        % Massa (kg)   
l=1;                        % Comprimento do pêndulo (m)    
grav=9.8;                   % Gravidade (m/s^2)
b=5;                        % Coeficiente de atrito viscoso Ns/m^2  
u0=m*grav*l*sin(50*pi/180)  % Valor do degrau aplicado
umax=1.2*u0                 % Saturação
Ts=0.01;                    % Tempo de amostragem    
N=300;                      % Número de pontos da resposta ao degrau 
Nu = 100;                     % Horizonte de Controle

% Iniciaização das Variáveis

x1(1)=0;
x2(1)=0;
u(1:N)=u0;

% Cálculo de resposta ao degrau
for t=1:N
    x1(t+1)=Ts*(-b*x1(t)-(grav/l)*sin(x2(t))+u(t)/(m*(l^2)))+x1(t);
    x2(t+1)=Ts*x1(t)+x2(t);
    h(t+1)=x2(t+1);
end

h=h*180/pi

lamb = 100;
tempo = 0:Ts:(N)*Ts;
plot(tempo,h)
ylabel('Posição (graus)')
xlabel('Tempo')
title('Resposta ao Degrau')
grid

g = h/u0;
Du0 = [u0;zeros(N-1,1)];
G = g(2:end);
x(1) = 0;
for i =1:N
    x(i+1) =  G*Du0;
    Du0 = [0;Du0(1:end-1)];
end
hold on
plot(tempo,x,'r--','linewidth',4)

legend('Resposta ao degrau', 'Resposta ao degrau modelado')

H(1) = G(1);
for i = 2:N
    H(i)= G(i)-G(i-1);
end

% matriz que multiplica o controle passado

Gf = zeros(N-1);
for linha = 1:N-1
    for coluna = 1:N-1
        if coluna ==1 
            Gf(linha,coluna) = G(linha+1);
        else
            if coluna+linha >N
                Gf(linha,coluna);
            else
                Gf(linha,coluna) = H(coluna+linha);
            end
        end
    end
end
Gf = [Gf;[G(N),zeros(1,N-2)]];

% Matriz dinâmica Gm

Gm = zeros(N);
for i = 1:N
    for j = 1:N
        if j>i
        Gm(i,j) = 0;
        else
            Gm(i,j) = G(i-j+1);
        end
    end
end

Gm = Gm(:,1:Nu);
M = inv(Gm'*Gm+lamb*eye(Nu))*Gm';
Km = M(1,:);

% Simula a resposta em malha fechada
tsim = 198;

clear y
clear x

y(1) = 0;              % saida inicial processo
x(1) = 0;              % saida inicial processo

up = [zeros(N-1,1)];  % controle passado

for t = 1:tsim+1;
    r(t) = 50;         % Referência Futura
    if t == tsim-10
        xx = 1;
    end
    w = ones(N,1)*r(t);  % trajetoria de referência futura = referência atual
    n(t) = y(t)-x(t);
    Mn = n(t)*ones(N,1);
    livre = Gf*up + Mn;  % Resposta livre   
    Du = Km*(w-livre);
    if t==1
        u(t) = Du;
    else
        u(t) = u(t-1) + Du;
    end
    
    % saturação
    
    if u(t)>umax
        u(t)=umax;
    end
    
    %------------
    x11(1)=0;
    x22(1)=0;

    % Simulação da planta
    
    x11(t+1)=Ts*(-b*x11(t)-(grav/l)*sin(x22(t))+u(t)/(m*(l^2)))+x11(t);
    x22(t+1)=Ts*x11(t)+x22(t);
    y(t+1)=x22(t+1)*180/pi;

    aux = 0;
   for i =1:N
       if t+1-i <=0
       else
       aux = aux+H(i)*u(t+1-i);
       end
   end
    x(t+1) = aux;
    
    up = [u(t);up(1:end-1)];
    
end
figure
subplot(2,1,1)
plot(y)
ylabel('Posição (graus)')
hold on

plot(r,'k--')
plot(x,'r--','linewidth',2)
title('Saídas')
legend('saída processo', 'saída modelo')
xlabel('amostras')
grid                                                                                                                                                                                                                                                                                                                                               
subplot(2,1,2)
plot(u,'linewidth',2)
title('Sinal de Controle')
ylabel('Torque (Nm)')
xlabel('amostras')
grid



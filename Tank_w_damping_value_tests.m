%% Simulation of tank
% Simulation of tank using euler forward, no delay or damping. Just a first
% order ODE.
%%
close all
clear all

% Importing last 3 tank flow data from dropbox folder, created and exported
% from makeCellMeassurements.m
%Inflowdata = importdata('Measurement data\SumInflowData_VEAS_EV_and_EV_VA.mat');
Inflowdata=load(fullfile('Measurement data','SumInflowData_VEAS_EV_and_EV_VA.mat'));
PASL_VA=load(fullfile('Measurement data','PASL25_VA.mat'));
FT_VA=PASL_VA.PASL25_VA/1000;

SUM_PASL_EV_VA = Inflowdata.SUM_PASL_EV_VA/1000;
d_v=SUM_PASL_EV_VA/1000;
SUM_PASL_VEAS_EV = Inflowdata.SUM_PASL_VEAS_EV/1000;
SUM_DIRECT_IPU = Inflowdata.SUM_DIRECT_IPU/1000;

IPU_FT=(load(fullfile('Measurement data','IPU_FT.mat')));
IPU_FT=IPU_FT.IPU_FT/1000;

IPU_FTsmooth = smooth(IPU_FT,60,'moving' );
IPU_FTsmooth = IPU_FTsmooth-0.9;        %scaling for plotting purposes
IPU_FT = IPU_FT-0.9;                    %scaling for plotting purposes
%%
% instantiating relevant varables
h=60;                %Step length
%N=length(d_m);       %Number of steps
N=60*6;
N=5500;

A=4000;             %Area of tank [m^2]
U=zeros(1,N);       %Control input u, outflow [m^3/s]
D=zeros(1,N);       %distubance
tVec=zeros(1,N);    %Time vector
R=zeros(1,N);       %Reference

%Tank dynamics, pure integrator
xDot=@(u,w) (w-u)/A;    % ODE
xSol=zeros(1,N);        % Solution of ODE

%Damping equation
%damping = 4000;
T=linspace(100,1100,11);
wDot=@(w,d,damping)(-w+d)/damping; 
wSol=zeros(1,N);

tests = 3;
shift = 210;
% Simulation of tank w. euler forward
for j=1:tests
    %Controller
    Kp=35;
    Ti=3000;
    %initial values
    startk = 12000;
    x=4;
    w=IPU_FTsmooth(startk+shift);
    t=0;
    u=3;
    r=4;
    z=0;
    u=0;
    if j==1
        damping=3780; % Winner of test 1 (k=1500:(1500+5500)
    end
    if j==2
        damping = 5920; % Winner of test 2 (k=5500:(5500+5500)
    end
    if j==3
        damping = 5240; % Winner of test 3 (k=1500:(1500+2*5500)
    end 
    %damping=5000+j*10;
    for i=1:N
        %store variables
        if x<0
            x=0;
        end
        tVec(i)=t;
        xSol(i)=x;
        wSol(j,i)=w;
        R(i)=r;
        D(i)=FT_VA(i+startk);
        U(i)=u;
        d=D(i);

        %solve ODEs
        w=w+wDot(w,d,damping)*h;
        x=x+xDot(u,w)*h;
        
        e=r-x;
        u=-Kp*(e+z);
        z=z+(e/Ti)*h;
        if u<=0
            u=0;
        end
        t=t+h;
    end
    
    figure(1)
    plot(tVec,IPU_FT(startk+shift:startk+N-1+shift,1),tVec,IPU_FTsmooth(startk+shift:startk+N-1+shift,1),tVec,wSol); hold on;
    legend('Scaled IPU flow','Smoothed IPU flow','Damper value 3780', 'Damper value 5920', 'Damper value 5240')
    grid minor, grid on
end

sipu = (IPU_FTsmooth(startk+shift:startk+N-1+shift,1))';

hold off
%% Check if the
error = length(sipu);
err = zeros(tests, N);
mse = zeros(1, tests)
for j=1:tests
    error = (abs([sipu-wSol(j,:)])).^2;
    err(j,:) = error;
    mse(:,j) = (sum(err(j,:)))/length(sipu);
end

[m, i] = min(mse);
fprintf('Best match is: ');
m
fprintf('Index number : ');
i
fprintf('Damping value is :');
%5000+i*10
    if i==1
        damping=3780 % Winner of test 1 (k=1500:(1500+5500)
    end
    if i==2
        damping = 5920 % Winner of test 2 (k=5500:(5500+5500)
    end
    if i==3
        damping = 5240 % Winner of test 3 (k=1500:(1500+2*5500)
    end 
    
%%
figure(2)
subplot(2,1,1)
plot(tVec,xSol)
legend('Tunnel Height [m]')
%set(gca,'ylim',[3.9 4.1]);
grid minor, grid on
subplot(2,1,2)
plot(tVec,U)
legend('Control input [m^3]')
grid minor, grid on
% subplot(2,1,2)
% plot(tVec,D,tVec,wSol,tVec,IPU_FT(startk+220:startk+N-1+220,1),tVec,IPU_FTsmooth(startk+220:startk+N-1+220,1))
% legend('Disturbance [m^3]','Damped disturbance [m^3]','Measured IPU flow', 'Smoothed measured IPU flow')
% grid minor, grid on








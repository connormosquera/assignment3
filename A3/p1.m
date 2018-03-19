C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.q_0 = 1.60217653e-19;             % electron charge

nElectrons = 10000;
nPlot=20; % number of electrons to actually plot
T = 300;
L = 200e-9;
W = 100e-9;
dt = 1e-15; % since 1/100 of 200nm is 2nm, smallest step allowed is 2nm/vth ~= 1e-14s
TStop = 1e-12; % 1000 timesteps
Vth = sqrt(2*C.kb*T/(C.m_0*0.26)); % using 2 degrees of freedom
time = 0;
Temp = T; % temperature variable that updates in TempCalc
taumn = 0.2e-12; % average time between collisions
sigmaMB = sqrt(C.kb*T/(C.m_0*0.26)); % standard deviation on vth
cc = jet(nPlot); % colorscale used to plot different electron colors

delVx = 0.8; % voltage difference in x direction
delVy = 0;

Ex = delVx/L;
Ey = delVy/W;

collisionT = zeros(200,nElectrons); % matrices for tracking collision time and velocities
collisionV = zeros(200,nElectrons);
collisionIndex = ones(1,nElectrons);
collisions = 0;

x = rand(1, nElectrons)*L; % assigning random initial particle positions
y = rand(1, nElectrons)*W;

Theta = rand(1, nElectrons)*2*pi; % selecting Vx and Vy from Gaussian centered at vth
Vx = cos(Theta).*(Vth + sigmaMB*randn(1, nElectrons));
Vy = sin(Theta).*(Vth + sigmaMB*randn(1, nElectrons));

avgV = sum(sqrt(Vx.^2+Vy.^2))/nElectrons % calculation of initial average velocity

currentD = zeros(1,1000);
currentDTime = linspace(1,1000,1000);

figure(2)
hFig2 = figure(2);
set(hFig2, 'Position', [0 0 1200 1000])

q=0;

for i=0:dt:TStop
    time = i;
    q=q+1;

%     subplot(2,1,1); % plotting electron positions
    hold on
%     for j=1:nPlot
%         plot(x(j), y(j), 'o','markers', 1, 'Color', cc(j,:));
%     end
    plot(x(1:30), y(1:30), 'bo','markers', 1);
    axis([0 L 0 W]);
    
    V2tot=Vx.*Vx+Vy.*Vy; % calculated temp based on total velocities
    KE = mean(V2tot)*0.5*(C.m_0*0.26);
    Temp = KE/C.kb;
    
%     subplot(2,1,2); % plotting temp vs. time
%     plot(time,Temp, 'ro', 'markers',1,'MarkerFaceColor', 'b');
%     hold on
%     axis([0 TStop 350 550]);
%     xlabel('Time(s)');
%     ylabel('Temp (K)');
%     grid on
%     title('Temperature vs. Time');
    
    x = x - dt * Vx; % moving the particles in one time step
    y = y - dt * Vy;
    
    for j=1:nElectrons % specular and periodic boundaries
        if x(j) > L
            x(j) = x(j) - L;
        elseif x(j) < 0
            x(j) = x(j) + L;
        end
        
         if y(j) > W
             Vy(j) = -Vy(j);
         elseif y(j) < 0
             Vy(j) = -Vy(j);
         end
    end
    
    Vx = Vx + Ex*C.q_0/C.m_0*dt; % second term is acceleration x time
    Vy = Vy + Ey*C.q_0/C.m_0*dt;
    
    for j=1:nElectrons % collision, mfp, and mean time between collisions tracking
        if (1-exp(-dt/taumn)) > rand()
            collisions = collisions+1;
            collisionT(collisionIndex(j)+1,j) = time;
            collisionV(collisionIndex(j)+1,j) = sqrt(Vx(j)^2+Vy(j)^2);
            collisionIndex(j)=collisionIndex(j)+1;
            
            Theta = rand(1, 1)*2*pi; % rethermalizing after collision
            Vx(j) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
            Vy(j) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
        end
    end
    
    currentD(q) = 10^15*C.q_0*mean(Vx);
    
    pause(0.000001)
    
end

figure(1) % plotting current density vs time
hold on
plot(currentDTime(1:1000), currentD(1:1000))
title('Current Density vs. Time');
hold off

MFP=0;
TBC=0;

for i=1:nElectrons % calculation of time between collisions (TBC) and mean free path (MFP)
    for j=1:collisionIndex(i)
        if j ~= 1
            TBC = TBC + collisionT(j,i)-collisionT(j-1,i);
            MFP = MFP + (collisionT(j,i)-collisionT(j-1,i))*(collisionV(j,i)-collisionV(j-1,i));
        end
    end
end

TBC = TBC/collisions
MFP = MFP/collisions

figure(3) % plotting electron density map in a 50x50 grid
hold on
n=hist3([x',y'],[50 50]);
pcolor(n');
colorbar;
title('Electron Density Map');
hold off

V5050 = zeros(50);

for h=1:nElectrons % calculating velocities for temperature calculation
    for i=1:50
        for j=1:50
            if x(h)>((i-1)/50*L) && x(h)<(i/50*L) && y(h)>((j-1)/50*W) && y(h)<(j/50*W)
                V5050(i,j)=Vx(h)^2+Vy(h)^2;
            end
        end
    end
end

for i=1:50 % taking average velocity per cell
    for j=1:50
       if n(i,j)~=0
          V5050(i,j) = V5050(i,j)/n(i,j);  
       else
          V5050(i,j) = 0; 
       end
    end
end

figure(4) % plotting temperature density
hold on
m=V5050.*0.5*0.26*C.m_0/C.kb;
pcolor(m');
colorbar;
title('Temperature Map');
hold off
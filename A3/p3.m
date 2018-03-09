clearvars
clearvars -GLOBAL
close all

L = 200;
W = 100;
Lb = 40;
Wb = 40;
Vo = 10;
maxI = 200;
delta = 1;

Sigma = ones(W,L);

for i=1:Wb
    for j=round(L/2-Lb/2):round(L/2+Lb/2)
        Sigma(i,j) = 0.01;
    end
end
for i=round(W-Wb):W
    for j=round(L/2-Lb/2):round(L/2+Lb/2)
        Sigma(i,j) = 0.01;
    end
end

G = sparse(L*W,L*W);
B = zeros(1,L*W);

for i=1:W
    for j=1:L
        n = j + (i-1)*L;
        
         if j==1
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = Vo;
        elseif j==L
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==1 && i==1
            G(n,:) = 0;
            G(n,n) = -(Sigma(i+1,j) + Sigma(i,j+1));
            G(n,n+1) = Sigma(i,j+1);
            G(n,n+L) = Sigma(i+1,j); 
        elseif j==L && i==1
            G(n,:) = 0;
            G(n,n) = -(Sigma(i+1,j) + Sigma(i,j-1));
            G(n,n-1) = Sigma(i,j-1);
            G(n,n+L) = Sigma(i+1,j);
        elseif j==1 && i==W
            G(n,:) = 0;
            G(n,n) = -(Sigma(i-1,j) + Sigma(i,j+1));
            G(n,n+1) = Sigma(i,j+1);
            G(n,n-L) = Sigma(i-1,j); 
        elseif j==L && i==W
            G(n,:) = 0;
            G(n,n) = -(Sigma(i-1,j) + Sigma(i,j-1));
            G(n,n-1) = Sigma(i,j-1);
            G(n,n-L) = Sigma(i-1,j);    
        elseif i==1
            G(n,:) = 0;
            G(n,n) = -(Sigma(i+1,j) + Sigma(i,j+1) + ...
                Sigma(i,j-1));
            G(n,n+1) = Sigma(i,j+1);
            G(n,n+L) = Sigma(i+1,j); 
            G(n,n-1) = Sigma(i,j-1);
        elseif i==W
            G(n,:) = 0;
            G(n,n) = -(Sigma(i,j+1) + Sigma(i,j-1) + ...
                Sigma(i-1,j));
            G(n,n-1) = Sigma(i,j-1);
            G(n,n+1) = Sigma(i,j+1); 
            G(n,n-L) = Sigma(i-1,j); 
        else
            G(n,:) = 0;
            G(n,n) = -(Sigma(i,j+1) + Sigma(i,j-1) + ...
                Sigma(i-1,j) + Sigma(i+1,j));
            G(n,n-1) = Sigma(i,j-1);
            G(n,n+1) = Sigma(i,j+1); 
            G(n,n-L) = Sigma(i-1,j); 
            G(n,n+L) = Sigma(i+1,j); 
        end
        
    end
end

F=G\B';

Vmap = zeros(W,L);

for i=1:W
    for j=1:L
        n = j + (i-1)*L;
        Vmap(i,j) = F(n);
    end
end

% figure(1)
% hold on
% surf(Vmap);
% title('Voltage Distribution')
% hold off
% pause(0.001);

[Ex,Ey]=gradient(-Vmap);

% figure(2)
% hold on
% quiver(Ex,Ey);
% axis([0 L 0 W])
% title('Electric Field (V/m)')
% hold off
% pause(0.001);

Ex=Ex'/(1e-9);
Ey=Ey'/(1e-9);

%%%%%%%%%%%%%%%%%%%%

C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.q_0 = 1.60217653e-19;             % electron charge

nElectrons = 10000;
nPlot=20; % number of electrons to actually plot
T = 300;
Lp = L*1e-9;
Wp = W*1e-9;
dt = 1e-15; % since 1/100 of 200nm is 2nm, smallest step allowed is 2nm/vth ~= 1e-14s
TStop = 1e-12; % 1000 timesteps
Vth = sqrt(2*C.kb*T/(C.m_0*0.26)); % using 2 degrees of freedom
time = 0;
Temp = T; % temperature variable that updates in TempCalc
taumn = 0.2e-12; % average time between collisions
sigmaMB = sqrt(C.kb*T/(C.m_0*0.26)); % standard deviation on vth
cc = jet(nPlot); % colorscale used to plot different electron colors

collisionT = zeros(200,nElectrons); % matrices for tracking collision time and velocities
collisionV = zeros(200,nElectrons);
collisionIndex = ones(1,nElectrons);
collisions = 0;

x = rand(1, nElectrons)*Lp; % assigning random initial particle positions
y = rand(1, nElectrons)*Wp;

for i=1:nElectrons % ensuring particles do not start in boxed boundaries
   while(1) 
      if ( x(i)<1.2e-7 && x(i)>0.8e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)) 
          x(i) = rand*Lp;
          y(i) = rand*Wp;
      else
        break
      end
   end
end

Theta = rand(1, nElectrons)*2*pi; % selecting Vx and Vy from Gaussian centered at vth
Vx = cos(Theta).*(Vth + sigmaMB*randn(1, nElectrons));
Vy = sin(Theta).*(Vth + sigmaMB*randn(1, nElectrons));

avgV = sum(sqrt(Vx.^2+Vy.^2))/nElectrons

figure(4)
hFig4 = figure(4);
set(hFig4, 'Position', [200 0 900 1000])

hold on
plot([0.8,0.8]*1e-7,[0,0.4]*1e-7, 'r-')
plot([0.8,0.8]*1e-7,[0.6,1]*1e-7, 'r-')
plot([1.2,1.2]*1e-7,[0,0.4]*1e-7, 'r-')
plot([1.2,1.2]*1e-7,[0.6,1]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[0,0]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[0.4,0.4]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[0.6,0.6]*1e-7, 'r-')
plot([0.8,1.2]*1e-7,[1,1]*1e-7, 'r-')

axis([0 Lp 0 Wp]);

for i=0:dt:TStop
    time = i;
    
    for j=1:nPlot
        plot(x(j), y(j), 'o','markers', 1, 'Color', cc(j,:));
    end
    
    
    V2tot=Vx.*Vx+Vy.*Vy; % calculated temp based on total velocities
    KE = mean(V2tot)*0.5*(C.m_0*0.26);
    Temp = KE/C.kb;
    
    x = x - dt * Vx; % moving the particles in one time step
    y = y - dt * Vy;
    
    for i=1:length(Vx)
        xE = round((x(i)*1e9+1)*200/201);
        yE = round((y(i)*1e9+1)*100/101);
        if xE>L
            xE=Lp;
        end
        if yE>W
            yE=Wp;
        end
        if xE<1
            xE=1;
        end
        if yE<1
            yE=1;
        end
        Vx(i) = Vx(i) + Ex(xE,yE)*C.q_0/C.m_0*dt;
        Vy(i) = Vy(i) + Ey(xE,yE)*C.q_0/C.m_0*dt;     
    end
   
    for j=1:nElectrons % specular and periodic boundaries
        if x(j) > Lp
            x(j) = x(j) - Lp;
        elseif x(j) < 0
            x(j) = x(j) + Lp;
        end
        
         if y(j) > Wp
             Vy(j) = -Vy(j);
         elseif y(j) < 0
             Vy(j) = -Vy(j);
         end
    end
    

    
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
    
    % the following code is uncommented if diffusive boundaries are desired
    %%%%%%%%% BlockBordersDiffusive Begin %%%%%%%%%%%
    
%     for i=1:nElectrons % rethermalized when hit boundary, theta defines scattering angle so it reflects away from boundary
%         if Vy(i)<0 && y(i)>0.6e-7 && y(i)<0.61e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
%             Theta = rand*pi;
%             Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
%             Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
%         elseif Vy(i)>0 && y(i)<0.4e-7 && y(i)>0.39e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
%             Theta = rand*pi+pi;
%             Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
%             Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
%         elseif Vx(i)<0 && x(i)>0.8e-7 && x(i)<0.81e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
%             Theta = rand*pi-pi/2;
%             Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
%             Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
%         elseif Vx(i)>0 && x(i)<1.2e-7 && x(i)>1.19e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
%             Theta = rand*pi+pi/2;
%             Vx(i) = cos(Theta)*(Vth + sigmaMB*randn(1, 1));
%             Vy(i) = sin(Theta)*(Vth + sigmaMB*randn(1, 1));
%         end
%     end
    
    %%%%%%%%% BlockBordersDiffusive End %%%%%%%%%%%
    

    % the following code is uncommented if specular boundaries are desired
    %%%%%%%%% BlockBorders Begin %%%%%%%%%%%
    
    for i=1:nElectrons % conditions for meeting a boundary, specular reflection by inverting x or y velocity
        if Vy(i)<0 && y(i)>0.6e-7 && y(i)<0.61e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
            Vy(i)=-Vy(i);
        elseif Vy(i)>0 && y(i)<0.4e-7 && y(i)>0.39e-7 && x(i)<1.2e-7 && x(i)>0.8e-7
            Vy(i)=-Vy(i);
        elseif Vx(i)<0 && x(i)>0.8e-7 && x(i)<0.81e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
            Vx(i)=-Vx(i);
        elseif Vx(i)>0 && x(i)<1.2e-7 && x(i)>1.19e-7 && (y(i)<0.4e-7 || y(i)>0.6e-7)
            Vx(i)=-Vx(i);
        end
    end
    
    %%%%%%%%% BlockBorders End %%%%%%%%%%%
    
    pause(0.001)
    
end

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

figure(5) % plotting electron density map in a 50x50 grid
hold on
n=hist3([x',y'],[50 50]);
pcolor(n');
colorbar;
title('Electron Density Map');
hold off

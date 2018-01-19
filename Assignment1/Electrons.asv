clearvars
clearvars -GLOBAL
close all

global nElectrons T L W MarkerSize
global x y Vx Vy C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²
C.am = 1.66053892e-27;

nElectrons = 10;
T = 300;
L = 100000;
W = 100000;
MarkerSize = 1;
dt = 0.001;
TStop = 10;
Vth = sqrt(C.kb*T/C.m_0);

x = (rand(1, nElectrons)-0.5)*L; % assigning random initial particle positions
y = (rand(1, nElectrons)-0.5)*W;
Vx = (rand(1, nElectrons)-0.5)*Vth;
Vy = (rand(1, nElectrons)-0.5)*Vth;

figure(1)

for i=1:dt:TStop
    % PlotElectrons([-L/2 L/2 -W/2 W/2]);
    plot(x, y, 'bo', 'markers',...
    MarkerSize,'MarkerFaceColor', 'b');
    hold on
    axis([-L/2 L/2 -W/2 W/2]);
    title('Electrons')
    xlabel('X')
    ylabel('Y')
    
    x = x - dt * Vx; % moving the particles in one time step
    y = y - dt * Vy;
    
    for j=1:nElectrons
        if x(j) > L/2
            x(j) = x(j) - L;
        elseif x(j) < L/2
            x(j) = x(j) + L;
        end
        
         if y(j) > W/2
             Vy(j) = -Vy(j);
         elseif y(j) < -W/2
             Vy(j) = -Vy(j);
         end
    end
    pause(0.01)
    
end
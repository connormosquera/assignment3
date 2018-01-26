clearvars
clearvars -GLOBAL
close all

global nElectrons T L W MarkerSize
global x y Vx Vy C time Temp dt

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
TStop = 1;
Vth = sqrt(C.kb*T/(C.m_0*0.26));
time = 0;
Temp = T;

InitElectrons;

figure(1)

for i=0:dt:TStop
    time = i;
    %PlotAll;
    PlotElectrons;
    
    TempCalc();
    
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
    pause(0.0001)
    
end
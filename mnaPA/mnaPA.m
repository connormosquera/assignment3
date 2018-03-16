G1=1/1;
C=0.25;
G2=1/2;
L=0.2;
G3=1/10;
alpha=100;
G4=1/0.1;
G0=1/1000;

G = [1 0 0 0 0 0 0;
    -G2 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -alpha 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+G0];

C = [0 0 0 0 0 0 0;
    -C C 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];

% V = [V1
%       V2
%       IL
%       V3
%       I3
%       V4
%       V0]
Vin = 0;

F = [Vin; 0; 0; 0; 0; 0; 0];

%%% DC SWEEP %%%

vinvec = zeros(1,7);
v0vec = zeros(1,7);
v3vec = zeros(1,7);

for i=1:21
    F = [i-11; 0; 0; 0; 0; 0; 0];
    V = G\F;
    
    vinvec(i) = i-11;
    v0vec(i) =  V(7);
    v3vec(i) = V(4);
end

figure (1)
plot(vinvec,v0vec)
hold on
grid on
plot(vinvec,v3vec)
title('DC Vin Sweep')
xlabel('Vin (V)')
ylabel('(V)')
legend('v0','v3')

%%% AC SWEEP %%%

omegavec = logspace(0,10,30);
v0vec = zeros(1,30);

F = [1; 0; 0; 0; 0; 0; 0];

for i=1:30
    
    V = (G+i*omegavec(i)*C)\F;
    
    v0vec(i) =  V(7);
end

figure (2)
loglog(omegavec,v0vec)
hold on
grid on
title('Frequency Sweep, showing V1')
xlabel('Omega')
ylabel('V1 (V)')

figure (3)
semilogx(omegavec,20*log10(v0vec))
hold on
grid on
title('Frequency Sweep showing Gain in dB')
xlabel('Omega')
ylabel('Gain (dB)')

%%% PERTURBATIONS %%%

dist = normrnd(.25,0.05,80,1)
v0vec = zeros(1,80);

for i=1:30
    C(2,1) = -dist(i);
    C(2,2) = dist(i);
    
    V = (G+i*pi*C)\F;
    
    v0vec(i) =  V(7);
end

figure (4)
histogram(real(20*log10(v0vec)),15)
hold on
grid on
title('Capacitor Variations')


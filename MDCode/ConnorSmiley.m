doPlot = 1;
dt = 5e-15;
TStop = 3000 * dt;
InitDist = 0.0;
Method = 'VE'; % VE -- verlot; FD -- Forward Difference

Mass0 = 14 * C.am; % Silicon
Mass1 = 100 * C.am; % Argon

AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing / 2^(1 / 6);
LJEpsilon = 1e-21;

PhiCutoff = 3 * AtomSpacing * 1.1;

T = 300;

AddHalfCircAtomicArray(10, 0, 0, 0, 0, 0, T, 1);
AddCircAtomicArray(3, -3*10^-9, 4*10^-9, 0, 0, 0, T, 1);
AddCircAtomicArray(3, 3*10^-9, 4*10^-9, 0, 0, 0, T, 1);
% vy0 = -sqrt(0.02*Ep/Mass1);
% AddRectAtomicArray(4,4,0,12*AtomSpacing,0,vy0,0,T,1);
Ep = 1;
AddParticleStream(5, 2, 10, -2 * pi / 3, 0, Ep * C.q_0, 7);


Size = 13*AtomSpacing;
Limits = [-Size +Size -Size +Size]; % square is good
PlDelt = 5 * dt;

PlotFile = 'BlockSt.gif';
PlotPosOnly = 1;
doPlotImage = 0;
PlotSize = [100, 100, 1049, 1049];

ScaleV = .02e-11;
ScaleF = 10;
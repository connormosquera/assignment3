function [x,y,Vx,Vy] = InitElectrons()
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)

%------------- BEGIN CODE --------------
global nElectrons T L W
global x y Vx Vy C

x = (rand(1, nElectrons)-0.5)*L; % assigning random initial particle positions
y = (rand(1, nElectrons)-0.5)*W;

Vth=67431;

Theta = rand(1, nElectrons)*2*pi;
Vx = cos(Theta)*Vth;
Vy = sin(Theta)*Vth;



%------------- END OF CODE --------------
%Please send suggestions for improvement of the above template header 
%to Denis Gilbert at this email address: gilbertd@dfo-mpo.gc.ca.
%Your contribution towards improving this template will be acknowledged in
%the "Changes" section of the TEMPLATE_HEADER web page on the Matlab
%Central File Exchange
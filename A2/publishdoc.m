%% ELEC 4700 Assignment 2 - Finite Difference Method
%% 1a - 1D Electrostatic Potential
% Below is a solution to the 2D potential, but represented as a 3D
% surface with a constant potential in the y direction. One side of the
% rectangle is held at a potential of 1, the other side is held at 0, and
% dV/dy is set to 0 for the other two sides to allow for the flat 1D-like
% potential.
% 

clearvars
clearvars -GLOBAL
close all

L = 100;
W = 3/2*L;
Vo = 1;

G = sparse(L*W,L*W);
B = zeros(1,L*W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        
        if i==L 
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = Vo;
        elseif i==1 
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==1
            G(n,:) = 0;
            G(n,n) = -1;
            G(n,n+1) = 1;     % dV/dy = 0      
        elseif j==W
            G(n,:) = 0;
            G(n,n) = -1; 
            G(n,n-1) = 1;     % dV/dy = 0    
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,n-1) = 1;
            G(n,n+1) = 1;
            G(n,n-W) = 1;
            G(n,n+W) = 1;
        end
        
    end
end

spy(G);

F=G\B';

Vmap = zeros(L,W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        Vmap(i,j) = F(n);
    end
end


surf(Vmap);
pause(0.001);

%% 1b - 2D Electrostatic Potential (Finite Difference Method)
% Below is the solution to the potential system where instead of 
% having dV/dy = 0 for the last two sides, they are instead fixed at 0.

clearvars
clearvars -GLOBAL
close all

W = 100;
L = 3/2*W;
Vo = 1;

G = sparse(L*W,L*W);
B = zeros(1,L*W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        
        if i==1 || i==L
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = Vo;
        elseif j==1
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==W
            G(n,:) = 0;
            G(n,n) = 1;
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,n-1) = 1;
            G(n,n+1) = 1;
            G(n,n-W) = 1;
            G(n,n+W) = 1;
        end
        
    end
end

F=G\B';

Vmap = zeros(L,W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        Vmap(i,j) = F(n);
    end
end

surf(Vmap);
pause(0.001);

%% 1b - 2D Electrostatic Potential (Analytic Method)
% Below is the same potential system as before, but solved using the 
% analytic solution provided by the lab document. 

clearvars
clearvars -GLOBAL
close all

W = 100;
L = 3/2*W;
Vo = 1;
Terms = 60;

V = zeros(L+1,W);

for k=1:Terms
    for j=1:W
        for i=1:L+1
            V(i,j) = V(i,j) + 4*Vo/pi*1/(2*k-1)*cosh((2*k-1)*pi*(i-L/2-1)/W)/cosh((2*k-1)*pi*(L/2)/W)*sin((2*k-1)*pi*j/W);
        end
    end
    surf(V);
    pause(0.05);
end

%% 1b - Discussion
% 60 terms were used for the analytic series solution method, as it was
% determined that using larger values of the index n resulted in an
% inf/inf calculation, resulting in a NaN output from Matlab. Both the
% analytic and finite difference methods produced approximately the same
% potential distributions. One of the downsides of the finite difference
% method is the rough edges on the corners of the map, which the analytic
% method seemed not to have. However, the analytic method experiences the
% Gibbs phenomenon due to the termination of the series solution, resulting
% in overshoots of the solution at the map boundaries.
%
%



%% 2 - Non-Uniform Conductivity Map
%
%
%

clearvars
clearvars -GLOBAL
close all

L = 100;
W = 3/2*L;
Lb = round(L/3);
Wb = round(W/3);
Vo = 1;
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

figure(1)
hold on
surf(Sigma);
axis([0 L 0 W])
hold off

G = sparse(L*W,L*W);
B = zeros(1,L*W);

for i=1:W
    for j=1:L
        n = j + (i-1)*L;
        
        if i==1 && j==L/2
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = Vo;
        elseif i==W && j==L/2
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==1 && i==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n+L) = 1/Sigma(i+1,j); 
        elseif j==L && i==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j-1));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+L) = 1/Sigma(i+1,j);
        elseif j==1 && i==W
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i-1,j) + 1/Sigma(i,j+1));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n-L) = 1/Sigma(i-1,j); 
        elseif j==L && i==W
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i-1,j) + 1/Sigma(i,j-1));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n-L) = 1/Sigma(i-1,j);    
        elseif j==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1) + ...
                1/Sigma(i-1,j));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n+L) = 1/Sigma(i+1,j); 
            G(n,n-L) = 1/Sigma(i-1,j);
        elseif j==L
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j-1) + ...
                1/Sigma(i-1,j));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+L) = 1/Sigma(i+1,j); 
            G(n,n-L) = 1/Sigma(i-1,j);
        elseif i==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1) + ...
                1/Sigma(i,j-1));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n+L) = 1/Sigma(i+1,j); 
            G(n,n-1) = 1/Sigma(i,j-1);
        elseif i==W
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i,j+1) + 1/Sigma(i,j-1) + ...
                1/Sigma(i-1,j));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+1) = 1/Sigma(i,j+1); 
            G(n,n-L) = 1/Sigma(i-1,j);  
        else
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i,j+1) + 1/Sigma(i,j-1) + ...
                1/Sigma(i-1,j) + 1/Sigma(i+1,j));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+1) = 1/Sigma(i,j+1); 
            G(n,n-L) = 1/Sigma(i-1,j); 
            G(n,n+L) = 1/Sigma(i+1,j); 
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

figure(2)
hold on
surf(Vmap);
hold off
pause(0.001);

[Ex,Ey]=gradient(-Vmap);

figure(3)
hold on
quiver(Ex,Ey);
axis([0 L 0 W])
hold off
pause(0.001);

Jx=Sigma.*Ex;
Jy=Sigma.*Ey;

Cin = sqrt(Jx(1,L/2)^2+Jy(1,L/2)^2);
Cout = sqrt(Jx(W,L/2)^2+Jy(W,L/2)^2);
Curr = (Cin+Cout)*0.5;

figure(4)
hold on
quiver(Jx,Jy);
axis([0 L 0 W])
% axis([0 W 0 L])
hold off
pause(0.001);





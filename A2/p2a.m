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
title('Conductivity Mapgit')
colorbar
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
title('Voltage Distribution')
colorbar
hold off
pause(0.001);

[Ex,Ey]=gradient(-Vmap);

figure(3)
hold on
quiver(Ex,Ey);
axis([0 L 0 W])
title('Electric Field (V/m)')
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
title('Current Density (A/m^2)')
hold off
pause(0.001);

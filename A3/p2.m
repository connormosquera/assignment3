clearvars
clearvars -GLOBAL
close all

L = 200;
W = 100;
Lb = 40;
Wb = 40;
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

figure(1)
hold on
surf(Vmap);
title('Voltage Distribution')
hold off
pause(0.001);

[Ex,Ey]=gradient(-Vmap);

figure(2)
hold on
quiver(Ex,Ey);
axis([0 L 0 W])
title('Electric Field (V/m)')
hold off
pause(0.001);
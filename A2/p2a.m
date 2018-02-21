clearvars
clearvars -GLOBAL
close all

L = 100;
W = 100;
Lb = round(L/3);
Wb = round(W/3);
Vo = 1;
maxI = 200;
delta = 1;

V = zeros(L,W);
V(L/2,W)=1;
Sigma = ones(L,W);

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

surf(Sigma);

G = sparse(L*W,L*W);
B = zeros(1,L*W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        
        if i==1 && j==W/2
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = Vo;
        elseif i==L && j==W/2
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==1 && i==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n+W) = 1/Sigma(i+1,j); 
        elseif j==W && i==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j-1));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+W) = 1/Sigma(i+1,j);
        elseif j==1 && i==L
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i-1,j) + 1/Sigma(i,j+1));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n-W) = 1/Sigma(i-1,j); 
        elseif j==W && i==L
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i-1,j) + 1/Sigma(i,j-1));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n-W) = 1/Sigma(i-1,j);    
        elseif j==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1) + ...
                1/Sigma(i-1,j));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n+W) = 1/Sigma(i+1,j); 
            G(n,n-W) = 1/Sigma(i-1,j);
        elseif j==W
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j-1) + ...
                1/Sigma(i-1,j));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+W) = 1/Sigma(i+1,j); 
            G(n,n-W) = 1/Sigma(i-1,j);
        elseif i==1
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i+1,j) + 1/Sigma(i,j+1) + ...
                1/Sigma(i,j-1));
            G(n,n+1) = 1/Sigma(i,j+1);
            G(n,n+W) = 1/Sigma(i+1,j); 
            G(n,n-1) = 1/Sigma(i,j-1);
        elseif i==L
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i,j+1) + 1/Sigma(i,j-1) + ...
                1/Sigma(i-1,j));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+1) = 1/Sigma(i,j+1); 
            G(n,n-W) = 1/Sigma(i-1,j);  
        else
            G(n,:) = 0;
            G(n,n) = -(1/Sigma(i,j+1) + 1/Sigma(i,j-1) + ...
                1/Sigma(i-1,j) + 1/Sigma(i+1,j));
            G(n,n-1) = 1/Sigma(i,j-1);
            G(n,n+1) = 1/Sigma(i,j+1); 
            G(n,n-W) = 1/Sigma(i-1,j); 
            G(n,n+W) = 1/Sigma(i+1,j); 
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

figure(2)
hold on
surf(Vmap);
hold off
pause(0.001);

[Ex,Ey]=gradient(Vmap);

% Ex=zeros(L,W);
% Ey=zeros(L,W);
% 
% for i=1:L
%     for j=1:W
%         if i==1
%             Ex(i,j) = Vmap(i+1,j) - Vmap(i,j);
%         elseif i==L
%             Ex(i,j) = Vmap(i,j) - Vmap(i-1,j);
%         else
%             Ex(i,j) = Vmap(i+1,j) - Vmap(i-1,j);
%         end
%         if j==1
%             Ey(i,j) = Vmap(i,j+1) - Vmap(i,j);
%         elseif j==W
%             Ey(i,j) = Vmap(i,j) - Vmap(i,j-1);
%         else
%             Ey(i,j) = Vmap(i,j+1) - Vmap(i,j-1);
%         end
%     end
% end

figure(3)
hold on
quiver(Ex,Ey);
hold off
pause(0.001);

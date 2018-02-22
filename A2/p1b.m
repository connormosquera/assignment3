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

%figure(1)
%hold on
surf(Vmap);
%axis([0 W 0 L]);
pause(0.001);
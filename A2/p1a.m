clearvars
clearvars -GLOBAL
close all

L = 100;
W = 100;
Vo = 1;
maxI = 20000;

V = zeros(L,W);
V(:,1)=ones(L,1);

Vtemp=V;

% for h=1:maxI
%     for i=2:L-1
%         for j=2:W-1
%             Vtemp(i,j) = 1/4*(V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1));
%         end
%         Vtemp(1,i) = Vtemp(2,i);
%         Vtemp(W,i) = Vtemp(W-1,i);
%     end
%     V=Vtemp;
% end

G=sparse(L*W,L*W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        
        if i==1 || i==L || j==1 || j==W
            G(n,:) = 0;
            G(n,n) = 1;
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,n-1) = 1;
            G(n,n+1) = 1;
            G(n,n-ny) = 1;
            G(n,n+ny) = 1;
        end
        
    end
end

spy(G);

[E,D] = eigs(G,9,'SM');

surf(V);
pause(0.001);

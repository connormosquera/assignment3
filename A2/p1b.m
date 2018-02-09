clearvars
clearvars -GLOBAL
close all

L = 100;
W = 100;
Vo = 1;
maxI = 20000;

V = zeros(L,W);
V(:,1)=ones(L,1);
V(:,W)=ones(L,1);

Vtemp=V;

for h=1:maxI
    for i=2:L-1
        for j=2:W-1
            Vtemp(i,j) = 1/4*(V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1));
        end

    end
    V=Vtemp;
end

surf(V);
pause(0.001);
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
    for j=round(L/2-Lb):round(L/2+Lb)
        Sigma(i,j) = 0.01;
    end
end
for i=round(W-Wb):W
    for j=round(L/2-Lb):round(L/2+Lb)
        Sigma(i,j) = 0.01;
    end
end

Vtemp=V;

for h=1:maxI
    for i=2:L-1
        for j=2:W-1
            Vtemp(i,j) = 1/(4*(Sigma(i,j)*delta^2))*(Sigma(i,j)*(V(i-1,j) ...
            + V(i+1,j) + V(i,j-1) + V(i,j+1))+ ...
            (Sigma(i+1,j)-Sigma(i-1,j))*(V(i+1,j)-V(i-1,j))+ ...
            (Sigma(i,j+1)-Sigma(i,j-1))*(V(i,j+1)-V(i,j-1)));
        end
    end
    V=Vtemp;
    surf(V);
    pause(0.001);
end



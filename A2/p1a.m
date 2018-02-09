clearvars
clearvars -GLOBAL
close all

% global 

C.q_0 = 1.60217653e-19;             % electron charge

L = 100;
W = 100;
Vo = 1;
maxI = 20000;

V = zeros(L,W);
V(:,1)=ones(L,1);

Vtemp=V;

for h=1:maxI
    for i=2:L-1
        for j=2:W-1
            Vtemp(i,j) = 1/4*(V(i-1,j) + V(i+1,j) + V(i,j-1) + V(i,j+1));
        end
        Vtemp(1,i) = Vtemp(2,i); 
        Vtemp(W,i) = Vtemp(W-1,i); 
    end
    V=Vtemp;
end

surf(V);
pause(0.001);

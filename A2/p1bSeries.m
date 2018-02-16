clearvars
clearvars -GLOBAL
close all

W = 100;
L = 3/2*W;
Vo = 1;
Terms = 30;

V = zeros(L,W);

for i=1:L
    for j=1:W
        for k=1:Terms
            V(i,j) = 4*V0/pi*1/k...;
        end
    end
end

surf(Vmap);
pause(0.001);
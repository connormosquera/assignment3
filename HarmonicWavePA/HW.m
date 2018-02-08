% (E(m+1,n) - 2*E(m,n) + E(m-1,n))/dx^2 + (E(m,n+1) - 2*E(m,n) + E(m,n-1))/dy^2 = a*E

clear all

size=50;
nx=size;
ny=size/2;
G=sparse(nx*ny,nx*ny);
B=zeros(1,nx*ny);

maxI = 500;
d=1;

for i=1:nx
    for j=1:ny
        n = j + (i-1)*ny;
        
        if i==1 || i==nx || j==1 || j==ny
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


Emap = zeros(nx,ny)
for i=1:nx
    for j=1:ny
        n = j + (i-1)*ny;
        Emap(i,j) = E(n,8);
    end
end

figure(1)
surf(Emap);




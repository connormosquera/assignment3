function [] = ConnorTraj2(nTraj)

    dt=1;

    x = zeros(1, nTraj);
    x1 = 0;
    t = 0;
    t1 = 0;
    Vx = zeros(1, nTraj);
    Vx1 = 0;
    Drv = 0;
    Drv1 = 0;
    
    dx = zeros(1, nTraj);
    
    n=200;

    for c=1:n

        dvx = dt;
        Vx(1, :) = Vx(1, :) + dvx;
        for i=1:nTraj 
            if rand(1,1)>0.95
                Vx(1,i)=0;
            end
        end
        dx(1, :) = Vx(1, :)* dt + dt^2 / 2;
            
        x(1, :) = x(1, :) + dx(1, :);
        t = t + dt;
            
        x1=[x1 sum(x)];
        t1=[t1 t];
        Vx1=[Vx1 sum(Vx)];
        
        Drv=sum(Vx1)/numel(Vx1);
        Drv1=[Drv1 Drv];
        
        subplot(2,1,1)
        plot(t1,x1,'r');
        hold on
        xlabel('t');
        ylabel('x');
        grid on
        title('Position vs. Time');
        
        subplot(2,1,2)
        plot(t1,Vx1,'r');
        hold on
        plot(t1,Drv1,'co:');
        xlabel('t');
        ylabel('Vx');
        grid on
        title('Velocity and Drift Velocity vs. Time');
        
        pause(0.01)
        
    end

end

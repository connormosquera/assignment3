function [] = ConnorTraj(nTraj)

    dt=1;

    x=0;
    x1=0;
    t=0;
    t1=0;
    Vx=0;
    Vx1=0;
    Drv=0;
    Drv1=0;
    
    n=200;

    for c=1:n

        dvx = dt;
        Vx = Vx + dvx;
        if rand(1,1)>0.95
            Vx=0;
        end
        dx = Vx * dt + dt^2 / 2;
            
        x = x + dx;
        t = t + dt;
            
        x1=[x1 x];
        t1=[t1 t];
        Vx1=[Vx1 Vx];
        
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

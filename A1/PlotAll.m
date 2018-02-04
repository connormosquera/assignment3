function [ output_args ] = PlotAll(Limits)
global nElectrons T L W MarkerSize
global x y Vx Vy time Temp dt

subplot(2,1,1)
hold on
PlotElectrons;
        
subplot(2,1,2)
plot(time,Temp, 'ro', 'markers',MarkerSize,'MarkerFaceColor', 'b');
hold on
axis([0 time+dt 250 1050]);
xlabel('Time(s)');
ylabel('Temp (K)');
grid on
title('Temperature vs. Time');


end
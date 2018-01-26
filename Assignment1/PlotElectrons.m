function [ output_args ] = PlotElectrons(Limits)
global nElectrons T L W MarkerSize
global x y Vx Vy

plot(x, y, 'bo', 'markers',MarkerSize,'MarkerFaceColor', 'b');
hold on
axis([-L/2 L/2 -W/2 W/2]);
xlabel('X');
ylabel('Y');
title('Electron Position');


end

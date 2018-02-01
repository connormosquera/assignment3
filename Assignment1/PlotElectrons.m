function [ output_args ] = PlotElectrons(cc)
global nElectrons T L W MarkerSize
global x y Vx Vy

plot(x, y, 'bo', 'markers',MarkerSize,'MarkerFaceColor', 'b');
%hold on
axis([-L/2 L/2 -W/2 W/2]);
%xlabel('X');
%ylabel('Y');
%title('Electron Position');


% color code

% for i=1:nElectrons
%     plot(x(i), y(i), 'o','markers', 1, 'Color', cc(i,:));
% end
% 
% axis([-L/2 L/2 -W/2 W/2]);


end

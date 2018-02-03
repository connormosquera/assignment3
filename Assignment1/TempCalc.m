function [] = TempCalc()
global nElectrons T L W MarkerSize
global x y Vx Vy C Temp

V2tot=Vx.*Vx+Vy.*Vy;
KE = mean(V2tot)*0.5*(C.m_0*0.26);
Temp = KE/C.kb;

end
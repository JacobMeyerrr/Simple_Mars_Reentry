function qAvg = HeatingRate(States)

R = vecnorm(States(:,1:3),2,2);

V = vecnorm(States(:,4:6),2,2);

rho = AtmDensityMars(R-3390);

Temp = Martian_Temp(R);

nu = KinVisc(Temp,rho);

S = Mars_SpeedofSound(Temp);

M = Mach_Number(V,S);

Re = Reynolds_Number(V,5,nu);

Cf = BA_SkinFric(Re,M);

qAvg = BA_HeatRate(rho,V,Cf);
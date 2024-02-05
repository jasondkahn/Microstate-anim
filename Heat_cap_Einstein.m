% Heat capacity plot
particles = 20000000;
numpts = 1000;
frospace = 1;
lowspace = 10;
midspace = 1000;
hispace = 35000;
hcvfrozen = micro_mod_collect(numpts,1,frospace,particles,'fast');
deltaTvf = hcvfrozen(2:numpts,1,10)-hcvfrozen(1:numpts-1,1,10);
% Makes a difference to use the midpoint temp as the x-axis!
Tavgf = hcvfrozen(1:numpts-1,1,10) + 0.5*deltaTvf;
hcvcold = micro_mod_collect(numpts,1,lowspace,particles,'fast');
deltaTvc = hcvcold(2:numpts,1,10)-hcvcold(1:numpts-1,1,10);
Tavgc = hcvcold(1:numpts-1,1,10) + 0.5*deltaTvc;
hcvmid = micro_mod_collect(numpts,1,midspace,particles,'fast');
deltaTvm = hcvmid(2:numpts,1,10)-hcvmid(1:numpts-1,1,10);
Tavgm = hcvmid(1:numpts-1,1,10) + 0.5*deltaTvm;
hcvhot = micro_mod_collect(numpts,1,hispace,particles,'fast');
deltaTvh = hcvhot(2:numpts,1,10)-hcvhot(1:numpts-1,1,10);
Tavgh = hcvhot(1:numpts-1,1,10) + 0.5*deltaTvh;
figure
hold on
plot(hcvhot(1:numpts,1,1),hcvhot(1:numpts,1,10),'mo')
plot(hcvmid(1:numpts,1,1),hcvmid(1:numpts,1,10),'rsquare')
plot(hcvcold(1:numpts,1,1),hcvcold(1:numpts,1,10),'gx')
plot(hcvfrozen(1:numpts,1,1),hcvfrozen(1:numpts,1,10),'b+')
set(gca,'FontSize',16)
xlabel('Energy(units of q)','FontSize',18)
ylabel('Temperature (units of thetaE)','FontSize',18)
figure
hold on
% ...uncomment to see the error from using the endpoint T
%plot(hchot(2:numpts,1,10),hispace./deltaTvh,'m+')
%plot(hcvmid(2:numpts,1,10),midspace./deltaTvm,'r+')
%plot(hcvcold(2:numpts,1,10),lowspace./deltaTvc,'g+')
%plot(hcvfrozen(2:numpts,1,10),frospace./deltaTvf,'b+')
plot(Tavgh,hispace./deltaTvh,'mo')
plot(Tavgm,midspace./deltaTvm,'rsquare')
plot(Tavgc,lowspace./deltaTvc,'gx')
plot(Tavgf,frospace./deltaTvf,'bo')
set(gca,'FontSize',16)
title("Model generates Einstein's Heat Capacity vs. temperature curve")
xlabel('Temperature T','FontSize',18)
ylabel('Heat Capacity dq/dT','FontSize',18)

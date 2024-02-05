function [entropyvecL,entropyvecR,tempervecL,tempervecR] = micro_mod_2Dcollect(QN,PN,pline)
% This function runs micro_movr3 and collects entropy and temperature
% values for 1:QN quanta being placed in the system on the left, with pline
% particles, with the remainder of the energy going to the right side, with
% PN-pline particles. Running micro_movr instead of just eval'ing the
% factorial is wasteful and slow, but one won't run this function many
% times.
% The list of values in the 3rd dimension is below. #3,5,7,and 9 will be
% blank or NaN, since the 'fast' option is used.
% 1  outstruct.q    quanta
% 2  outstruct.p    particles
% 3  outstruct.Wd   # of microstates - direct calc
% 4  outstruct.WE   # of microstates - from Einstein solid (p+q-1)!/[(p-1)!q!]
% 5  outstruct.Sd   entropy - direct calc
% 6  outstruct.SE   entropy - from Einstein solid model
% 7  outstruct.NCE  # of configurations - direct calc or from integer
%                   partition function
% 8  outstruct.N0E  # of particles with zero energy in the ensemble - from
%                 the area of an exponential given known TE from the thermo dE/dS
% 9  outstruct.Td   T from the fit to an exponential to direct calc occupancies
% 10 outstruct.TE   T from from the thermodynamic (1/T) = dE/dS
%
% The qlist is in the outv structure but is never collected
% outstruct.ql   the qlist, the list of configurations, is provided if q < p 
%                for possible use in a subsequent calculation, otherwise it is 0
clear S_vs_qL
clear S_vs_qR
S_vs_qL = zeros(QN-1,1,10);
S_vs_qR = zeros(QN-1,1,10);
%speed = 'fast';
%if strcmp(speed,'fast')
    ndisp = 1;
%elseif strcmp(speed,'slow')
%    ndisp = 0;
%elseif strcmp(speed,'quiet') % Not implemented yet, may not make any difference
%    ndisp = 3;
%else
%    ndisp = 1;
%end
    
%outv = 0;
for qnum = 1:QN-1
    clear outv
%    disp(['Working on q = ',num2str(qnum*qspacer),' p = ',num2str(pnum)])
% The evalc commands pipe all the disp output from micro_movr3 into the
% dummy res variable so it doesn't clutter up the screen. If you want to
% see the output, use the line below instead
    res = evalc('outv = micro_movr3(qnum,pline,''pc'','''',ndisp);');
%    outv = micro_movr3(qnum,pline,'pc','',ndisp);
%    disp(['Temp = ',num2str(outv.TE)])
    S_vs_qL(qnum,pline,:) = [outv.q  outv.p  outv.Wd  outv.WE ...
            outv.Sd outv.SE outv.NCE  outv.N0E  outv.Td ...
            outv.TE];
    res = evalc('outv = micro_movr3(QN-qnum,PN-pline,''pc'','''',ndisp);');
%    outv = micro_movr3(QN-qnum,PN-pline,'pc','',ndisp);
    S_vs_qR(qnum,pline,:) = [outv.q  outv.p  outv.Wd  outv.WE ...
            outv.Sd outv.SE outv.NCE  outv.N0E  outv.Td ...
            outv.TE];
end
% This evals the combined system
disp('Complete system results:');
comboout = micro_movr3(QN,PN,'pc','',ndisp);
entropyvecL = S_vs_qL(:,pline,6);
tempervecL = S_vs_qL(:,pline,10);
entropyvecR = S_vs_qR(:,pline,6);
tempervecR = S_vs_qR(:,pline,10);
[peakent,qmax] = max(entropyvecL+entropyvecR);
f = figure;
f.Position(4) = 2*f.Position(4) ;
subplot(2,1,1)
hold on
plot(1:QN-1,entropyvecL,1:QN-1,entropyvecR,1:QN-1,entropyvecL+entropyvecR)
plot(qmax,peakent,'o',qmax,comboout.SE,'x','MarkerSize',14)
box on
title('Entropy vs. energy distribution','FontSize',14)
xlabel('# of quanta on the left side')
ylabel('Eyrtpon')
legend({'Left side','Right side','Total','','Combined systems'},'Location','east')
t1 = ['Partition is at position ',num2str(pline,'%d'),' of ',num2str(PN,'%d')];
t2 = ['Total quanta = ',num2str(QN,'%d')];
%sprintf('\n');
%disp(newline);
t3 = ['Peak total S of ',num2str(peakent),' at q = ',num2str(qmax),' on the left'];
disp([newline t3])
disp(['Component entropies are ',num2str(entropyvecL(qmax)),' and ',num2str(entropyvecR(qmax))])
text(0.05,0.65,t1,'Units','normalized')
text(0.05,0.59,t2,'Units','normalized')
text(0.05,0.53,t3,'Units','normalized')
[~,qmax] = min(abs(tempervecL-tempervecR));
disp(['Left temp = ',num2str(tempervecL(qmax)),', right temp = ',num2str(tempervecR(qmax))]);
subplot(2,1,2)
hold on
plot(1:QN-1,tempervecL,'o',1:QN-1,tempervecR,'o')
plot(qmax,comboout.TE,'ks','MarkerSize',14,'LineWidth',2)
legend('Left side','Right side','Combined systems')
box on
title('Temperature vs. energy distribution','FontSize',14)
xlabel('# of quanta on the left side')
ylabel('Temperature')
%legend
%plot(1:QN-1,([1:QN-1]/pline)-((QN-[1:QN-1])/(PN-pline)))

    
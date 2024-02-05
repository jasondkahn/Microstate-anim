function S_vs_q2 = micro_mod_collect(QN,PN,qspacer,pspacer,speed)
% This function runs micro_movr3 and collects the results in a 3-D array called
% S_vs_q2, where S_vs_q2(qnum,pnum,n) contains the nth value from the list below,
% from the output structure from micro_movr3
% Typical values for input parameters:
% QN or PN = 50; %Careful here or it will take forever!
% qspacer = 5
% pspacer = 10
% speed = 'fast'  % or set to 'slow', or anything else defaults to nodispflag = 2
% 'fast' does not do any direct enumeration. It does slow down if it tries
% to do the partition function calculation to estimate the number of
% configurations.
% The first two indices of S_vs_q2 correspond to qnum and pnum
% So, for example, to plot Entropy vs. q at the 7th value of p, one could do
% plot(S_vs_q2(:,7,1),S_vs_q2(:,7,6))
% The list of values in the 3rd dimension is below. #3,5,7,and 9 will be
% blank or NaN if the 'fast' option is used.
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
%S_vs_q = zeros(N+1,3);
%deltaS = zeros(N+1,1);
clear S_vs_q2
S_vs_q2 = zeros(QN,PN,10);
if strcmp(speed,'fast')
    ndisp = 1;
elseif strcmp(speed,'slow')
    ndisp = 0;
elseif strcmp(speed,'quiet') % Not implemented yet, may not make any difference
    ndisp = 3;
else
    ndisp = 1;
end
    
%outv = 0;
for qnum = 1:QN
    clear outv
    pnum = 1;
    disp(['Working on q = ',num2str(qnum*qspacer),' p = ',num2str(pnum*pspacer)])
    outv = micro_movr3(qnum*qspacer,pnum*pspacer,'pc','',ndisp);
    disp(['Temp = ',num2str(outv.TE)])
    S_vs_q2(qnum,pnum,:) = [outv.q  outv.p  outv.Wd  outv.WE ...
            outv.Sd outv.SE outv.NCE  outv.N0E  outv.Td ...
            outv.TE];
    %qnum*spacer
    outv_old = outv;
    close all
    for pnum = 2:PN
    disp(['Working on q = ',num2str(qnum*qspacer),' p = ',num2str(pnum*pspacer)])
    outv = micro_movr3(qnum*qspacer,pnum*pspacer,'pc','',ndisp,outv_old);
    disp(['Temp = ',num2str(outv.TE)])
    S_vs_q2(qnum,pnum,:) = [outv.q  outv.p  outv.Wd  outv.WE ...
            outv.Sd outv.SE outv.NCE  outv.N0E  outv.Td ...
            outv.TE];
    outv_old = outv;
%    deltaS(qnum+1) = S_vs_q(qnum+1,2) - S_vs_q(qnum,2);
    close all
    end
end

    
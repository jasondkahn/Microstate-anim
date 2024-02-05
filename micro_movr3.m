function outstruct = micro_movr3(varargin)
% Usage: outputstructure = micro_movr3(quanta,particles,pcflag,toppc,nodispflag,inputstructure);
% Illustrates and animates configurations and microstates for a toy model
%   of distributing energy among particles, inspired by Nash's ChemThermo,
%   essentially a 1-D simplified version of the Einstein model for a solid.
%   Calculates W, T, S, and other results.
% Direct enumeration is SLOW for big q and p numbers! Set the nodispflag
%   flag to 1 to forgo looking at individual configurations and simply
%   output values calculated from the Einstein solid model formulae. The
%   program will automatically and progressively scale back its ambitions
%   for large p and q to avoid calculating numbers that exceed realmax.
% Arguments: micro_movr3(quanta,particles,pcflag,toppc,nodispflag,inputstructure);
%   All arguments are optional but must be given in order if used, some
%     conflict with each other, some can be omitted. Ace programmer, here!
%   Default >> micro_movr3; is equivalent to >> micro_movr3(3,5,'',8,0,''); 
%   Set quanta and particles and allow animations as in >> micro_movr3(6,8);
%   set pcflag to 'pc' to show only predominant configurations, as in
%     >> micro_movr3(13,20,'pc');
%   set toppc to change the number of configurations to display and average
%     over in 'pc' view (default toppc = 8). Takes a numerical argument
%     (best args are 8, 15, 31, 6n-1 for n >= 6) as in 
%     >> res3 = micro_movr3(17,24,'pc',15,'',''); 
%     In animation mode, toppc controls how many configurations to iterate
%     through, as in >> micro_movr3(12,8,[],31);
%     In pc mode, or if p and q are too large and the program goes
%     into pc mode automatically, it will use the value of toppc. So for example 
%     >> micro_movr3(17,24,'',15,'',''); and micro_movr3(17,24,'pc',15,'','');
%     will give the same result.
%   set nodispflag = 0 or omit it or enter [] to allow animations if possible
%   set nodispflag = 1 to just get a calculation of W, S, N0, and T as in
%     >> micro_movr3(47,20,'','',1); . If nodispflag is 1, then we still try to count the
%     total number of configurations in lieu of enumerating each of them,
%     if there are not too many configurations. However, we do not generate graphs of 
%     configurations, so nodispflag = 1 overrides 'pc' or '' in the previous argument. 
%   Set nodispflag = 2 to get a plot of one representative predominant
%     configuration. This routine assumes the Boltzmann distribution is
%     true and calculates a configuration that is usually very close to the
%     most populated. Works for very large numbers of quanta and particles.
%     For example, >> micro_movr3(57,30,'','',2)
%   [Should implement nodispflag = 3 for shutting off all "disp" statements.
%     If desired, could add more logic that uses the actual value of
%     nodispflag.]
%   set inputstructure to be the output structure (outstruct) from a previous
%     run, to save time by reusing a partition list (qlist) under some
%     conditions. For example try the lines below. 
%     res1 = micro_movr3(37,42,'','',0); % But note if nodispflag is set to 1 here,
%        there will be no qlist generated
%     res2 = micro_movr3(37,45,'','',0,res1);  % This will find and re-use the
%        qlist from the res1 structure, if it still exists
% (RFE : the partition list generation is the slowest step for big problems. Should 
%   filter it and write out qlists for all smaller problems. All of the configs for
%   q = 50 and p = 30 are contained in the q = 50 and p = 40 list. So for
%   every q,p that the user enters could just generate ql's for all smaller
%   values of p. Somehow they would have to get stored for re-use, which would
%   require a lot of storage and bookkeeping.)
% RFE: implement an absolutely-no-plots output for use within scripts. For
%   now, be careful to close plots after making too many of them or the system 
%   will bog down, using >> close all .
% Jason Kahn University of Maryland College Park
%   4/11/2021  version 1.6  fixed help text above, minor changes and formatting,
%     improved messages to user, improved plot labeling, error checking for
%     compatibility with pre-2020 versions.
%
% The program outputs a structure with the following fields:
% outstruct.q    quanta
% outstruct.p    particles
% outstruct.Wd   # of microstates - direct calc
% outstruct.Sd   entropy - direct calc
% outstruct.NCd  # of configurations - direct calc
% outstruct.N0d  direct count of particles with zero energy
% outstruct.Td   T from the fit to an exponential to directly calculated occupancies
% outstruct.WE   # of microstates - from Einstein solid formula (p+q-1)!/[(p-1)!q!]
% outstruct.SE   entropy - from Einstein solid model formula
% outstruct.NCE  # of configurations from integer partition function
% outstruct.N0E  # of particles with zero energy in the ensemble - from
%                 the area of an exponential given known TE from the thermo dE/dS
% outstruct.TE   T from from the thermodynamic T = dE/dS
% outstruct.ql   the qlist, the list of configurations, is provided if q < p 
%                for possible use in a subsequent calculation, otherwise it is 0
% 
%   micro_movr3 calls JDK's functions placequanta, calcdist, partitionnum, partitionnumrec, 
%     partitionfct, pentagonalnums, and RoundSum -- these must all be in the Matlab Path as well
%   Also requires "partitions" function from the Matlab File Exchange:
%   This one: http://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer
%   Try/catches are there to use cprintf if it is available, to print color to the command window.
%
%   Shows all of the configs and microstates for "quanta" units of energy
%       randomly allocated into "particles" boxes.
%   Animations pause often -- hit spacebar in the plot window to continue, if you hit
%       the spacebar too many times will bring the front window back and you 
%       may need to bring plot window to the front again to continue. You
%       can always control-C out from the command window.
%   Can generate a lot of plots -- get rid of them with "close all"
%   Caution: We are counting microstates! Will take forever for large
%       numbers of particles! You're asking for trouble if you enumerate
%       when either quanta or particles is  > 50. In this case the program
%       may decide to take shortcuts on its own.
%   For just calculating W, uses expression for the Einsten model of a solid:
%   Weinsolid = factorial(particles+quanta-1)/(factorial(quanta)*factorial(particles-1));
%       which will be non-infinite for particles+quanta-1 < 171, otherwise
%       uses some tricks to avoid evaluating the factorial and still gives W, or just 
%       uses the gammaln function gammaln(n) = log(factorial(n-1)) to report
%       just the entropy, proportional to ln(W).
%   T is in units of the Einstein Temperature (ThetaE). T is thus 1/(the exp
%       decay constant) from the Boltzmann distribution.
%
% Defaults just so something will happen if you simply type "micro_movr3":
quanta = 3;
particles = 5;
% turnoffzoom turns off the auto-zooming that I find so annoying in plots 
turnoffzoom = 1;
% Set up numbers of plots and configurations to work up.
% The numbers below can be changed -- crank them up with care!
toppc = 8;  % toppc is the max number of configs to work up and display in pc mode 
% If toppc > 8, may get tiny subplots -- see subplot arrangements below.
% The best values for toppc are 8, 15, 31, or 6n-1, n >= 6. It may be
% really slow and crowded for toppc > 65. Should put logic in to question
% toppc > 100.
% Might need to mess with the size of the config graph below. This will
% depend on details of monitor size.
config_graph_position = [2 2 18 13];
% # of microstates to show before going into fast animation
num_ustates_stepthrough = 5;
% Should not need to mess with stuff below unless you are trying to push
% the limits of characterizing large models.
% topcon is actually set at around line 155 below once we have the final
% value for toppc
%topcon = max(17,toppc); % Topcon is the max number of rows of qlist to
% display in any display mode
% The default "topcon" # is easily adjusted with this one line
animconfiglim = 30;
% Program will query before showing >= animconfiglim animated configs
% The program will refuse to list or animate more than pconlyconfiglim
% configs, where the default limit is 30 or toppc, whichever is larger.
% If you feel limited by this, change the code to output the qlist and work with it.
% Similarly pconlyconfiglim is set below at around line 155 with final
% value of toppc
%pconlyconfiglim = max(30,toppc);
% Must have toppc <= topcon <= pconlyconfiglim and animconfiglim < pconlyconfiglim
% Sumlim is set to the max argument for the factorial function
% If quanta+particles > sumlim, just return calculated values without going
%through the time-consuming task of partitioning to generate configurations
% Sumlim is the largest number whose factorial is non-infinite + 1
sumlim = 171;
% Configlim is the limit for direct enumeration of configs. If evals take
% forever under direct enumeration, decrease this number.
% So far the only way to display the predominant configurations is to
% enumerate them and pick the winners. I should implement a way to find
% the toppc best without enumeration, but they will all just look like the
% Boltzmann distribution anyway.
configlim = 1e7;
%
%-----------There should not be any need to change anything below----------
%
% The program defaults to showing microstate animations as long as there
% are not too many configurations
pcflag = '';  
% The inputql may be replaced by an output structure from a previous run
inputql = 0;
% When nosdispflag is nonzero it suppresses all graphs e.g. for collecting W vs q stats
% (Controlled from the command line)
nodispflag = 0;
% For grabbing command line inputs:
% Put defaults into the optargs cell array
optargs = {quanta particles pcflag toppc nodispflag inputql};
maxnumopts = 6;
% Boltzmann's constant (also def in WS_Einstein_solid function but because
% it is not a nested function the value is not available to the main
% program.
% kb is now exact as of 2018 recommendations. K(elvin) is derived
kb = 1.380649e-23 ; % Joule / K
%
% Read in varargin to change default values if user provides new values
%
numvarargs = length(varargin);
if numvarargs > maxnumopts
    error('micro_mov:TooManyInputs', ...
        'micro_mov allows at most %d optional inputs',maxnumopts);
end
%
% Now put the new input values from varargin into the optargs cell array, 
% overwriting the default values.
%
% But if the user has entered '' then we do not want to replace the default
% There is probably a more elegant way to do this.
aretheythere = cellfun(@isempty,varargin);
for ivarg = 1:length(varargin)
  if ~aretheythere(ivarg)
    optargs(ivarg) = varargin(ivarg);
  end
end
%optargs(1:numvarargs) = varargin;
[quanta, particles, pcflag, toppc, nodispflag, inputql] = optargs{:};
% Reset config lims here based on final value of toppc:
if ~isnumeric(toppc) || (length(toppc) > 1) || (rem(toppc,1) ~= 0)
 disp('The value for toppc must be a single integer, no quotes')
 return
end
topcon = max(17,toppc);
pconlyconfiglim = max(30,toppc);
% Logic to provide variables used below to decide which plots to create
% and to check the value of pcflag before we dive into calculations
if ~isempty(pcflag)   
    if strcmp(pcflag,'pc')
% The pf flag means the pc flag has been set, so no animations just configs
        pf = 1;
% [This logic has evolved over time, should re-simplify. Can always
% control-C to bail out of animations]
    else
        disp('The only permissible arguments for pcflag, the third entry, are ''pc'' or '' '' or []')
        return
    end
else
    pf = 0;
end
% Could add more interpretation of strings like the following:
if strcmp(nodispflag,'nd')
    nodispflag = 1;
end
% Baseline checking to bail gracefully and informatively
if (particles == 0)
    disp('If particles = 0, there is no there, there.');
    return
end
if (quanta == 0)
    disp('If quanta = 0, then T = 0 K, S = 0, and there isn''t much more to say.');
    return
end
if (particles < 3)
    disp('Warning: If particles = 1 or 2, results are non-physical.');
end
%
% Now start calculations
%
[W_einsolid, S_einsolid] = WS_Einstein_solid(particles,quanta);
% Above 1000 the partitionfct calculation is unreliable.
% If nodispflag is 1 we do not bother to calculate the number of partitions
% This logic has not been fully checked out.
% partitionfct is recursive and can be very slow for big numbers...although
% not nearly as slow as actually enumerating the partitions!
if (nodispflag == 1)
    numconfigscalc = NaN;
elseif (quanta+particles) < 1000
    numconfigscalc = partitionfct(quanta,particles);
else
    numconfigscalc = NaN;
    if nodispflag ~= 3
        disp('For "big" systems, # of configs is unreliable, not reported.')
    end
end
toppc = min(numconfigscalc,toppc);
% Adjust topcon downward if toppc is more than the number of configs
if toppc <= 8
    spx = 3;
    spy = 3;
    spavg = 9;
elseif toppc <= 15
    spx = 4;
    spy = 4;
    spavg = 16;
elseif toppc <= 31
    spx = 4;
    spy = floor(toppc/4)+1;
    spavg = 4*spy;
else
    spx = 6;
    spy = floor(toppc/6)+1;
    spavg = 6*spy;
end
% T_calc from the analytical expression for 1/T = (dS/dE) eval at const N
%   (particles). Using Stirling's approximation:
% T_stirling = 1./log((quanta+particles-1)/quanta)
% Using Psi function = d/dx(gammaln(x)) where gamma(x) = factorial(x-1)
T_calc = 1./(psi(quanta+particles)-psi(quanta+1));
% N0E just comes from the integrated area of an exponential assuming T_calc
% The total number of particles must be W*particles
% This overestimates NOE because the quanitization of energy
% throws off the integration of the area.
N0E = W_einsolid*particles/T_calc;
%N0E_stirling = W_einsolid*particles/T_stirling
if nodispflag ~= 3
disp(' ') % skip a line
disp(['Distributing ',num2str(quanta),' quanta among ',num2str(particles),' particles'])
disp(['T calculated from T = (dE/dS)|const N is ',num2str(T_calc),...
    ' (units of Einstein Temp ThetaE)'])
disp(['N0 (particles at zero energy) from normalization to area of exponential is ',...
    num2str(N0E)])
disp(['Number of configs from partition fxn calc is '...
    , num2str(numconfigscalc),' configurations'])
disp(['Number of microstates from Einstein solid calc is ', num2str(W_einsolid)])
disp(['Entropy from Einstein solid calc is ', num2str(S_einsolid),' J/K'])
end
outstruct = struct('q',quanta,'p',particles,'WE',W_einsolid,...
        'SE',S_einsolid, 'NCE',numconfigscalc, 'N0E',N0E,'TE', T_calc,'ql',0);
bail = 0;
% This would be a good place to put a graph of the calculated Boltzmann
% distribution PC.
% Burying the code for nodispflag == 2 here is terrible programming.
if ((particles + quanta) > sumlim || numconfigscalc > configlim || nodispflag == 2)
    if ((particles + quanta) > sumlim || numconfigscalc > configlim)
        if nodispflag ~= 3
        disp('Problem is beyond direct enumeration.')
        end
    end
% If N0E is infinite then we don't plot the Boltzmann distributions
    if (nodispflag == 2)
        theox1 = 0:1:quanta+1; % remove semi for debug
        bfit_pq = [num2str(quanta),' quanta distributed among ',...
                num2str(particles),' particles'];
        if (N0E ~= Inf)
            cfig3 = figure;
            if (~verLessThan('matlab','9.5') && turnoffzoom)
                disableDefaultInteractivity(gca);
            end
            cfig3.Units = 'inches';
            calcboltzpos = cfig3.Position;
            hold on
            box on
% Axis goes to quanta+1 to make sure actual # of particles with
% that energy goes to zero
            plot(theox1,N0E*exp(-theox1/T_calc),'r','LineWidth',2)
%   Line below is for Bose-Einstein statistics and probably wrong  
%    plot(theox,expA_BE.*(1./(exp(theox./T_BE) - 1 )),'m','LineWidth',1)
            title('Calculated Boltzmann Distribution','FontSize',14)
            xlabel('Number of quanta (E)','FontSize',12)
% Plot twenty reasonable xticks
            if quanta <= 25
                xtickinc = 1;
            elseif quanta <= 50
                xtickinc = 2;
            elseif quanta <= 125
                 xtickinc = 5;
            else
                 xtickinc = 10*ceil(ceil(quanta/20)/10);
            end     
            xticks(0:xtickinc:quanta+1)
            ylabel('Number N of particles with energy E','FontSize',12)
            lgdB = legend('Continuous model best at high T');
            lgdB.FontSize = 14;
            V = axis;
            bfit_11 = [num2str(W_einsolid,6),' total microstates'];
            bfit_11a = [num2str(numconfigscalc),' configurations'];
            bfit_12 = ['N = ' num2str(N0E,6) ' * exp(-E/',num2str(T_calc),')'];
            bfit_13 = ['Entropy is ',num2str(S_einsolid),' J/K'];
            text(.4*quanta,.75*V(4),bfit_pq,'FontSize',14);
            text(.4*quanta,.67*V(4),bfit_11,'FontSize',14);
            text(.4*quanta,.59*V(4),bfit_11a,'FontSize',14);
            text(.4*quanta,.51*V(4),bfit_12,'FontSize',14);
            text(.4*quanta,.43*V(4),bfit_13,'FontSize',14);
        end
% Repeats code developed below with modifications. Should create a separate
% function?
% plot config differently
% Can't let axes get out of hand but these limits on q and p are arbitrary
% This draws a representative PC. Doesn't work very well for q > (p^2/2)
% (roughly) because then every particle can have a unique energy and there
% may be a very large number of such configurations, which means that
% averages are not very good. Should program a way to average over all of
% the configs with exactly the same W.
        if (quanta < 200000 && particles < 2000)
            idealqvsp = (particles/T_calc)*exp(-theox1/T_calc); % remove semi for debug
% idealconfignums is the number of particles with each value of q from
% quanta+1 down to zero (for consistency with distribm)
            tempconfig = RoundSum(idealqvsp);
            % renormalize. This can be off because calculated temperatures
            % from enumeration and the einstein model do not agree
            % perfectly
            % remove semi for debug
            idealconfignums = fliplr(RoundSum(idealqvsp*particles/sum(tempconfig))); 
%            initidealcfnums = idealconfignums; %for debugging of how much
%            it changes in renormalization of the energy below
%           I have not observed the renormalization of particle #s to fail
            if (sum(idealconfignums) ~= particles)
                disp('Renormalization of ideal config particle # is off.')
            end
% Total energy may be a bit off also and the code below is a very ad hoc way
% of fixing it, so some populations may be slightly non-optimal, but where it
% can be compared to an enumerated set of configs it does very well.
% It takes one of the zero-energy particles and pops it up as high as
% possible without introducing gaps, then takes the highest-energy energy
% state in the next cohort of abundance and pops it up as much as possible,
% etc. At the end it adds a quantum to the highest levels so the only gaps
% are at the highest energy and they are rare.
            total_energy = dot(idealconfignums,quanta+1:-1:0);
            err_energy = quanta - total_energy;
            if (err_energy ~= 0)
                if (err_energy > 0)
                    % maxqni = quanta + 2 - find(idealconfignums,1);
                    % remove semi to debug
                    promonum = find(idealconfignums == max(idealconfignums),1,'first');
%                    promonum = quanta + 2;
                    rem_eng = err_energy; %remove semi to debug
                    while (rem_eng > 0)
                        %promonum;
                        %remove semi to debug
                        idealconfignums(promonum) = idealconfignums(promonum) - 1; 
%                        if maxqni >= rem_eng
                            %ediff = promonum - rem_eng; %uncomm to debug
                            %remove semi to debug
                            pop_of_energy_to_go_to = idealconfignums(promonum-rem_eng); 
                            %remove semicolon from next line to debug
                            level_to_go_to = ...
                                find(idealconfignums==pop_of_energy_to_go_to,1,'last'); 
                            if level_to_go_to == promonum
                                idealconfignums(promonum) = idealconfignums(promonum) + 1;
                                % We want the last of the contiguous zeros?
                                %remove semi to debug
                                highest_level = find(idealconfignums==0,1,'last'); 
                                idealconfignums(highest_level) = 1; 
                                idealconfignums(highest_level+1) = ...
                                    idealconfignums(highest_level+1) - 1;
                                rem_eng = rem_eng - 1; %remove semi to debug
                                promonum = find(idealconfignums == ...
                                    max(idealconfignums),1,'first');
%                                promonum = quanta + 2;
                            elseif level_to_go_to <= promonum
                                idealconfignums(level_to_go_to) = ...
                                    idealconfignums(level_to_go_to) + 1;
                                % Should never do anything to increase the
                                % rem_eng
                               
                                %rem_eng = rem_eng - (promonum - level_to_go_to); 
                                %remove semis from next three lines to debug
                                rem_eng = rem_eng - (promonum - level_to_go_to); 
                                pop_of_new_energy_to_come_from = ...
                                    max(1,idealconfignums(promonum - 1)); 
                                promonum = find(idealconfignums == ...
                                    pop_of_new_energy_to_come_from,1,'first'); 
                            else
                                idealconfignums(promonum) = idealconfignums(promonum) + 1;
                                last_nonzero = find(idealconfignums,1,'last');
                                last_zero_to_use = ...
                                    find(idealconfignums(1:last_nonzero)==0,1,'last');
                                idealconfignums(last_zero_to_use) = ...
                                    idealconfignums(last_zero_to_use) + 1;
                                idealconfignums(last_zero_to_use+1) = ...
                                    idealconfignums(last_zero_to_use+1) - 1;
                                rem_eng = rem_eng - 1;
                            end
                            % if block below is for when we are near the high-energy end
                            if (idealconfignums(promonum) == 1 || rem_eng > 0) 
                                idealconfignums(promonum-1) = ...
                                    idealconfignums(promonum-1) + 1;
                                idealconfignums(promonum) = idealconfignums(promonum) - 1;
                                rem_eng = rem_eng - 1;
                                %remove semi to debug
                                promonum = ...
                                    find(idealconfignums == max(idealconfignums),1,'first'); 
%                                promonum = quanta + 2; % reset
                            end
                      % Promote the highest level that is after the first zero-pop level
                            if rem_eng == 1 
                                
                                last_nonzero = find(idealconfignums,1,'last');
                                last_zero_to_use = ...
                                    find(idealconfignums(1:last_nonzero)==0,1,'last');
                                idealconfignums(last_zero_to_use) = ...
                                    idealconfignums(last_zero_to_use) + 1;
                                idealconfignums(last_zero_to_use+1) = ...
                                    idealconfignums(last_zero_to_use+1) - 1;
                                rem_eng = 0;
                            end
                    end
                end
                if (err_energy < 0) % This should never happen
                    idealconfignums(quanta+2) = idealconfignums(quanta+2) + 1;
                    idealconfignums(quanta+2+err_energy) = ...
                        idealconfignums(quanta+2+err_energy) - 1;
                end
                new_total_energy = dot(idealconfignums,quanta+1:-1:0);
%                initidealcfnums - idealconfignums
            if (new_total_energy ~= quanta)
                disp(['Renormalization of ideal config energy is off. Energy is ',...
                    num2str(new_total_energy)])
            end
                
            end
            maxqni = quanta + 2 - find(idealconfignums,1);
            cfig4 = figure;
            if (~verLessThan('matlab','9.5') && turnoffzoom)
                disableDefaultInteractivity(gca);
            end
            cfig4.Units = 'inches';
            if (N0E ~= Inf)
                cfig4.Position = [calcboltzpos(1)+1,calcboltzpos(2)-1,...
                    calcboltzpos(3),calcboltzpos(4)];
            end
            hold on
            box on
% each iteration of the for loop plots a row of green dots
%            weightavg = weightavg + distribm(ifig,:)*nummicros(currentconfig);
%            runsummicro = runsummicro + nummicros(currentconfig);
            for qni = find(idealconfignums,1):quanta+2
                % qz is the number of particles with zero quanta, counting
                % down by subtraction
%                qz = qz - distribm(ifig,qn);
                if (idealconfignums(qni) > 0)
                    dots = (quanta+2-qni+0.25)*ones(idealconfignums(qni),1);
                    plot(1:idealconfignums(qni),dots,'o',...
                        'MarkerSize',12,...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor',[0,0.75,0]);       
                end
            end
% This plots the green dots at zero energy
%            dots = (0.25)*ones(qz,1);
%            plot(1:qz,dots,'o',...
%                'MarkerSize',12,...
%            	'MarkerEdgeColor','k',...
%            	'MarkerFaceColor',[0,0.75,0]);
            %axis([0.5 particles+0.5 -.25 quanta+0.5])
            axis([0.5 particles+0.5 -.25 min(quanta+0.5,maxqni+2.5)])
% Now plot a random microstate from the config as an example.
% Swap in the line below to show the first microstate
%        plot(1:particles,partsvsnumq,'ro',...
%
        mqni = 1;
        idealpartsvsnumq = zeros(particles,1);
        for kqni = 1:quanta+2   % kqni is the number of quanta in a particle,
            % counting down from quanta+1
            for l = 1:idealconfignums(kqni)  % l is the number of times it shows up
% Looks like we are calculating this for all configs and we
% don't need to, since we are only printing the top ones, and
% that only if we are displaying?
% But the profiler suggests that the time wasted is immaterial
                idealpartsvsnumq(mqni) = quanta-kqni+2; 
                mqni=mqni+1;
            end
        end

            plot(randperm(length(1:particles)),idealpartsvsnumq,'ro',...
                'MarkerSize',12,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0.75,0,0]);
%        axis([0.5 particles+0.5 -.25 quanta+0.5])
%            set(idcfg,'YTick',0:quanta) 
%            set(idcfg,'Xtick',0:particles)
            xlabel('Particle Number','FontSize',12)
            ylabel('Quanta per Particle','FontSize',12)
%            if (mod(ifig-1,spy) == 0 || ~pf)
       tipc = title('Predominant Configuration (green) with random microstate (red)');
%                if pf
            tipc.FontSize = 14;
%                else
%                    tipc.FontSize = 10;
%                end
            text(particles*0.15+1,min(quanta+0.5,maxqni+2.5)*0.85,...
                bfit_pq,'FontSize',14)
        	t3a = text(particles*0.15+1,min(quanta+0.5,maxqni+2.5)*0.95,...
                'Representation of a PC or neighboring config');
            t3a.FontSize = 14;
            %idealconfignums
            %pzero = idealconfignums(quanta+2);
            % Calculate W for the PC. Use "prod" so it will calculate the
            % N!/(N0!) without calculating either one individually
            %numdivfactors = min(particles-pzero,quanta+2);
%numer_numideal = (prod((particles:-1:(particles-numdivfactors))./...
%   factorial(idealconfignums(quanta+1:-1:);
            %denom_numideal = prod(factorial(idealconfignums(2:quanta+1)))
            %if numer_numideal = Inf
            %    t4a = text(particles*0.15+1,min(quanta+0.5,maxqni+2.5)*0.75,...
            %    ['W for this PC is probably ',num2str(numer_numideal)]);
            %numidealPCmicros = prod(particles:-1:pzero+1)/...
            %     prod(factorial(idealconfignums(2:quanta+1)));
           %if ~isnan(numidealPCmicros)
        	%t4a = text(particles*0.15+1,min(quanta+0.5,maxqni+2.5)*0.75,...
            %    ['W for this PC is ',num2str(numidealPCmicros)]);
           % t4a.FontSize = 14;
%disp(['Displaying sample PC comprising ', num2str(numidealPCmicros),' microstates'])
           % else
                %try harder to calculate the true number, for now just say
                %it is infinite. Problem is that Inf/Inf comes out NaN
                % With this logic I am reasonably sure that Inf is really >
                % 10^306
                [~,~,nzi] = find(idealconfignums(1:quanta+1));
                snzi = sort(nzi,'descend');
                %length(snzi)
                wfac = particles;
                numidealPCmicros = 1;
                for jsnzi = 1:length(snzi)
                        lnum = snzi(jsnzi);
                        factor = prod((wfac:-1:wfac-lnum+1)./(lnum:-1:1));
                        numidealPCmicros = numidealPCmicros*factor;
                        wfac = wfac-lnum;
                end
                %Catch possible NaN if something has escaped me
                if isnan(numidealPCmicros)
                    numidealPCmicros = Inf;
                end
               t4a = text(particles*0.15+1,min(quanta+0.5,maxqni+2.5)*0.75,...
                ['W for this PC is ',num2str(numidealPCmicros)]);
                t4a.FontSize = 14; 
                disp(['Displaying sample PC comprising ', ...
                    num2str(numidealPCmicros),' microstates'])
           % end
        end
%         t3 = text(particles*0.15,quanta*0.85,['Config ', num2str(i), ' of ',...
%             num2str(numconfigsdirect)]);
%         t2 = text(particles*0.35,quanta*0.65,[num2str(nummicros(i),'%4.3e'),...
%             ' microstates']);
    end
    bail = 1;
end
% Could change this logic to enumerate here, but it is likely to be very slow and
% it should not give different answers than the Einstein calculation
if (nodispflag)
    if nodispflag ~= 3
disp('Problem is supra-enumerable or user chose no display, so do not enumerate NC or W.')
disp('  The output structure has NC and W values calculated from Einstein solid model.')
    end
    bail = 1;
end
if bail == 1
    outstruct.Wd = NaN;
    outstruct.Sd = NaN;
    outstruct.NCd = NaN;
	outstruct.N0d = NaN;
	outstruct.Td = NaN;
    % not reporting a qlist
    return
end
% Calculate the number of configurations
% If the mumber is too large, evaluate whether to bail on the animated graphs.
% The partitions function usage is
% plist = partitions(total_sum,candidate_set,max_count,fixed_count)
% We need all partitions with a fixed_count that ranges from 1 up
%  to the number of quanta or particles, whichever is smaller.
% For large quanta this saves a lot of time not generating useless
% partitions, vs. generating all partitions and then filtering.
% If q < p all partitions will have at least one particle with no quanta,
% and if this is the case any larger value of p will generate the same
% set of partitions, just with more zero-energy particles, so provide an option to
% save the qlist for re-use.
% Flip the order returned by partitions so that the first column has the largest
% number of quanta in one particle ...the rot90 with the 2 flips by 180
if (isa(inputql,'struct') && isfield(inputql,'ql'))
    disp('Attempting to re-use qlist')
    if (inputql.q == quanta && inputql.p >= quanta && particles > quanta ...
            && (size(inputql.ql,2) == quanta))
        qlist = inputql.ql;
        try
            cprintf('blue', 'Re-using supplied qlist, hooray.\n')
        catch
            disp('Re-using supplied qlist, hooray.')
        end
    else
       disp('Supplied qlist is unsuitable, sorry.')
    end
else
    %preallocate qlist. If the prealloc is wrong the program still works
	qlist = zeros(numconfigscalc,quanta);
    % The first config is all the energy in one particle. Do this line
    % first for debugging and setting dimensions of qlist
    qlist(1,:) = rot90(partitions(quanta,1:quanta,[],1),2); 
% This next loop is the exponentially increasing part of the problem
% because the partitions function is recursive. [Maybe there is some smarter
% way to save and re-use partitions?] Because of the way this loop is set
% up, the predominant configurations tend to emerge near the end but not at
% the very end of the qlist -- in these configs, energy is distributed over a
% large number of particles without piling every particle into the bottom
% levels.
%
% j ranges over the mumber of particles that receive energy.
% Therefore partitions returns an empty matrix if j > quanta. So just avoid
% calling it and save a few microseconds.
    qk = 2;
    for j = 2:min(quanta,particles)
        qlfixedj = rot90(partitions(quanta,1:quanta,[],j),2);
        qlist(qk:qk+size(qlfixedj,1)-1,:) = qlfixedj;
        qk = qk+size(qlfixedj,1);
    end
end
%
% We already know the number of configs but we gather it explicitly here to
% make sure all methods agree.
%
numconfigsdirect = size(qlist,1); % the number of columns in qlist is always = quanta
disp(['There are ',num2str(numconfigsdirect),' directly enumerated configurations'])
%
% Logic to provide variables used below to decide which plots to create
%
if (numconfigsdirect > pconlyconfiglim)
   % nodispflag
    if (~nodispflag && pf ~= 1)
        try
            cprintf('blue','Executive override: too many graphs, showing PCs\n')
        catch
            disp('Executive override: too many graphs, showing PCs')
        end
    end
    pf = 1;
end
if ~nodispflag
    %disp(['There are ',num2str(numconfigsdirect),' configurations'])
    % We are displaying something, so let's find out about the default
    % figure position
    defpos = num2cell(get(0,'defaultfigureposition')); 
    [defx,defy,defsizex,defsizey] = defpos{:};
    if (numconfigsdirect > animconfiglim && pf == 0)
prompt = 'Do you really want to display all of the microstate animations, Y/N [N]? ';
        cresp = input(prompt,'s');
        if isempty(cresp)
            cresp = 'N';
        end
        if (cresp == 'Y' || cresp == 'y')
        	disp('Okay, you asked for it.')  % go ahead
        end
        if (cresp ~= 'Y' && cresp ~= 'y')
        	pf = 1;
        	disp('Probably a wise choice.')
        end
    end
end
%
% Now directly enumerating the number of microstates
%
nummicros = zeros(size(numconfigsdirect));
boltz = zeros(size([qlist(1,:),0]));
% Work through configurations
for iconf=1:numconfigsdirect
% Count how many microstates there are in the configuration
% And collect data on how many particles have each energy
% Distrib is a row vector with k entries like for example
%    [ 0 1 1 0 0 ] for 1x4quanta and 1x3quanta and no other excitations
    distrib = qlist(iconf,:);
%     try
    nummicros(iconf) = factorial(particles)/...
        (prod(factorial(distrib))*factorial(particles-sum(distrib)));
%     catch
%         distrib
%         particles
%         sum(distrib)
%         factorial(distrib)
%     end
    distwzero = [distrib , particles-sum(distrib)];
    boltz = boltz + nummicros(iconf)*distwzero;
end
%sum(boltz)
Wdirect = sum(nummicros);
disp(['There are ',num2str(Wdirect),' directly enumerated microstates'])

%[W_einsolid, S_einsolid] = WS_Einstein_solid(particles,quanta);
%Weinsolid = factorial(particles+quanta-1)/(factorial(quanta)*factorial(particles-1));
% At the end we will plot out the development of the Boltzmann
% distribution. Do fitting here to extract "T"
% The built-in exponential fit works a lot better than fminsearch
fitexp=fit((quanta:-1:0)',boltz','exp1');
expA = fitexp.a;
T_direct = -1/fitexp.b;
disp(['T from direct exponential fit is ',num2str(T_direct), ' ThetaE'])
% Bose-Einstein fit. This is here for no good reason, since the particles
% in the Einstein solid are distinguishable. 
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'StartPoint',[expA -fitexp.b]);
% %               'Lower',[0,-Inf],...
% %               'Upper',[Inf,0],...
% bose_einstein = fittype('a*(1/(exp(b*x)-.999))','options',fo);
% [fitexp2,~] = fit((quanta:-1:1)',boltz(1:end-1)',bose_einstein)
% expA_BE = fitexp2.a
% T_BE = 1/fitexp2.b
% outstruct = struct('q',quanta,'p',particles,'W',W_einsolid,...
%         'S',S_einsolid, 'NC',numconfigscalc, 'N0',expA,'T', T_calc,'ql',0);
%
% If q < p then adding more particles does not change the available
%   configurations, so output the qlist for subsequent use. 
% Could just do this automatically but it might use a lot of memory.
%
outstruct.Wd = Wdirect;
outstruct.Sd = kb*log(Wdirect);
outstruct.NCd = numconfigsdirect;
outstruct.N0d = expA;
outstruct.Td = T_direct;
if quanta < particles
    outstruct.ql = qlist;
    if nodispflag ~= 3
disp('Can re-use the qlist in the future, with the same q, via micro_movr3(q,p,'''','''',0,ans_of_this_run)')
    end
end
%
% Now plot configurations
%
% select the top "topcon" configurations to limit display length in pc mode
% pcis is a list of the configs with the most microstates
% We don't know which are the top ones until we sort (could we do this above while
% collecting the nummicros list?)
%
[topclist,pcis] = sort(nummicros,'descend');
% cut down the list
topconfigs = pcis(1:min(topcon,length(pcis)));
displayconfigs = pcis(1:min(pconlyconfiglim,length(pcis)));
maxqn = 1;
if ~nodispflag
disp('Configuration List: Each line gives the number of particles that have each number of quanta in');
    disp(['  descending in q, so 0 3 2 4 means no particles with 3 quanta, 3 p with 2 q,'...
        ' 2 with 1, 4 with none']);
    disp(['Top-ranked configurations are listed in order of which configs ',...
        'have the most microstates (ustates)']);
    disp(['A configuration''s ID number, e.g. # ',num2str(topconfigs(1)),...
        ', is not meaningful except that because of the way the'])
    disp('  partitions routine works, later configs in the overall set tend to be more populated');
    % 
    if verLessThan('matlab','9.1')
        %if 1  % To test the avoidance of the strings command
        disp('Display of configs requires "strings" command, disabled for pre-R2016b')
        disp('I have not tested the consequences thoroughly.')
    else
    tcfarray = strings(min(topcon,length(pcis)),1);
    Warray = strings(min(topcon,length(pcis)),1);
    qlarray = zeros(min(topcon,length(pcis)),quanta+1);
    for dj = 1:min(topcon,length(pcis))
        tcfarray(dj,1) = string(num2str(topconfigs(dj)));
        Warray(dj,1) = string(num2str(topclist(dj)));
        qlarray(dj,1:quanta) = qlist(topconfigs(dj),:);
        pz = particles - sum(qlist(topconfigs(dj),:));
        qlarray(dj,quanta+1) = pz;
    end
    tcfpadded = pad(tcfarray,'left');
    Warraypadded = pad(Warray,'left');
    qlarraypadded = strings(min(topcon,length(pcis)),quanta+1);
    for qj = 1:quanta+1
        qlarraypadded(:,qj) = pad(string(qlarray(:,qj)),'left');
    end
    for dj = 1:min(topcon,length(pcis))
% Need to pad the array of config number strings for visual alignment. The %6u flag does
% not seem to work
% pz is particles with zero energy
%        pz = particles - sum(qlist(topconfigs(dj),:));
%         disp(['Config # ', char(tcfpadded(dj)),', with ',char(Warraypadded(dj)),...
%             ' ustates, has  P vs Q: ', num2str(qlist(topconfigs(dj),:)),...
%             '  ',num2str(pz) ]);
        disp(['Config # ', char(tcfpadded(dj)),', with ',char(Warraypadded(dj)),...
            ' ustates, has  P vs Q: ', char(strjoin(qlarraypadded(dj,:)))]);
        % Calculating maxqn to do y-axis scaling on plots
        
    end
    end
    %move the line below outside the version control block.
    for dj = 1:min(topcon,length(pcis))
        maxqn = max(maxqn,quanta+1-find(qlist(topconfigs(dj),:),1));
    end
    partW = sum(topclist(1:min(topcon,length(topclist))));
    disp(['These top ', num2str(min(topcon,length(pcis))),' configurations out of ',...
        num2str(numconfigsdirect),', or ',...
        num2str(100*min(topcon,length(pcis))/numconfigsdirect),' %, comprise '])
    disp(['  ',num2str(partW),' microstates out of a total W of ',num2str(Wdirect),...
        ' enumerated microstates, or ',...
        num2str(100*partW/Wdirect),' %.'])
    disp(['Total W calculated by Einstein solid method is ',num2str(W_einsolid)])
    disp(['The entropy S = kb*ln(W) is ',num2str(S_einsolid),' J/K'])
% Shorten the pcis list of configs to be no more than the first toppc elements
    toppclist = pcis(1:min(toppc,length(pcis)));
% Distribm is a matrix for the distributions among top configs -- needed to
% set consistent axes for config plots and for other plots.
% Pull out the best pconlyconfiglim entries in the qlist matrix
    distribm = zeros(min(numconfigsdirect,pconlyconfiglim),size(qlist,2));
            runsummicro = 0;
        weightavg = 0;
    if (~nodispflag && pf)
        % cfig is the handle for the multiplot PC figure
        cfig = figure;
        cfig.Units = 'inches';
        cfig.Position = config_graph_position;
    end
    %[numconfigsdirect,pconlyconfiglim,length(toppclist)]
    numfigs = min([numconfigsdirect,pconlyconfiglim,length(toppclist)]);
    % Sigh -- previously we were iterating through pconlyconfiglim whether
    % we needed to or not -- no reason to go beyond the toppclist
    for ifig=1:numfigs
% Distrib is a row vector with k entries like [ 0 1 1 0 0 ] for 1x4 and 1x3 quanta
% Preserved here to make microstate plots work rather than rewriting to just use
% distribm and indexing in as needed...laziness
        distribm(ifig,:) = qlist(displayconfigs(ifig),:);
        currentconfig = displayconfigs(ifig);
        distrib = distribm(ifig,:);  
        m = 1;
% partvsnumq is an array of how many quanta for particle 1, particle 2 etc.
% sorted in descending order of energy e.g. it would be [4 3 0 0 0] for
% above distrib
        partsvsnumq = zeros(particles,1);
        for k = 1:quanta   % k is the number of quanta in a particle, counting down
            for l = 1:distribm(ifig,k)  % l is the number of times it shows up
% Looks like we are calculating this for all configs and we
% don't need to, since we are only printing the top ones, and
% that only if we are displaying?
% But the profiler suggests that the time wasted is immaterial
                partsvsnumq(m) = quanta-k+1; 
                m=m+1;
            end
        end
% Uncomment below to print out partsvsnumq for debugging
%partsvsnumq
%
% ~isempty(toppclist(toppclist==currentconfig)) checks whether the config
% we are working on is on the list of top pc's
% toppclist % 
        if (~nodispflag && (~pf || ~isempty(toppclist(toppclist==currentconfig))))
            if (~pf)
                % This will be the left side of a two-subplot figure
                microfigpos = [50+ifig*50 defy-(ifig-1)*50 defsizex defsizey];
                figure('Position',microfigpos);
                gcb = subplot(1,2,1);
                if (~verLessThan('matlab','9.5') && turnoffzoom)
                    disableDefaultInteractivity(gca);
                end
            else
                figure(cfig)
%                spx
%                spy
                gcb = subplot(spx,spy,ifig);
                if (~verLessThan('matlab','9.5') && turnoffzoom)
                    disableDefaultInteractivity(gca);
                end
            end
% plot config differently
            hold on
            box on
            qz = particles;
% each iteration of the for loop plots a row of green dots
            weightavg = weightavg + distribm(ifig,:)*nummicros(currentconfig);
            %ifig
            runsummicro = runsummicro + nummicros(currentconfig);
            for qn = 1:quanta
                % qz is the number of particles with zero quanta, counting
                % down by subtraction
                qz = qz - distribm(ifig,qn);
                if (distribm(ifig,qn) > 0)
                    dots = (quanta+1-qn+0.25)*ones(distribm(ifig,qn),1);
                    plot(1:distribm(ifig,qn),dots,'o',...
                        'MarkerSize',12,...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor',[0,0.75,0]);       
                end
            end
% This plots the green dots at zero energy
            dots = (0.25)*ones(qz,1);
            plot(1:qz,dots,'o',...
                'MarkerSize',12,...
            	'MarkerEdgeColor','k',...
            	'MarkerFaceColor',[0,0.75,0]);
            %axis([0.5 particles+0.5 -.25 quanta+0.5])
            axis([0.5 particles+0.5 -.25 min(quanta+0.5,maxqn+2.5)])
% Now plot a random microstate from the config as an example.
% Swap in the line below to show the firat microstate
%        plot(1:particles,partsvsnumq,'ro',...
            plot(randperm(length(1:particles)),partsvsnumq,'ro',...
                'MarkerSize',12,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0.75,0,0]);
%        axis([0.5 particles+0.5 -.25 quanta+0.5])
            set(gcb,'YTick',0:quanta) 
            set(gcb,'Xtick',0:particles)
            xlabel('Particle Number','FontSize',12)
            ylabel('Quanta per particle','FontSize',12)
            if (ifig == spy || ~pf || (ifig == length(toppclist) && length(toppclist) < spy))
                    tipc = title([num2str(quanta),' quanta distributed among ',...
                            num2str(particles),' particles']);
                if pf
                    tipc.FontSize = 14;
                    % The axes-level alignment commands are new in R2020b
                    % Newer versions of Matlab want you to use
                    % "isMatlabReleaseOlderThan"...but that command has only
                    % worked sicne 2020b
                    if verLessThan('matlab','9.9') %#ok<*VERLESSMATLAB>
                        % Do nothing here -- right justification too much
                        % of a pain
                        %tipc.HorizontalAlignment = 'right';
                    else
                        ax3 = gca;
                        ax3.TitleHorizontalAlignment = 'right';
                    end
                else
                    %don't think there is a path to this statement.
                    tipc.FontSize = 11;
                end
            end
            %if (mod(ifig-1,spy) == 0 || ~pf)
            if (ifig == 1 || ~pf)
                if pf
                    tipc = title('Configurations in green with a random microstate for each in red');
            % The axes-level alignment commands are for R2020b (9.9) and later
                    if verLessThan('matlab','9.9')
                        tipc.HorizontalAlignment = 'left';
                    else
                        axrand = gca;
                        axrand.TitleHorizontalAlignment = 'left';
                    end
                    tipc.FontSize = 14;
                else
                    tipc = title('Configuration in green with random microstate in red');                  
                    tipc.FontSize = 11;
                end
            end
%         t3 = text(particles*0.15,quanta*0.85,['Config ', num2str(i), ' of ',...
%             num2str(numconfigsdirect)]);
%         t2 = text(particles*0.35,quanta*0.65,[num2str(nummicros(i),'%4.3e'),...
%             ' microstates']);
            t3 = text(particles*0.15+1,min(quanta+0.5,maxqn+2.5)*0.85,...
                ['Config ', num2str(currentconfig), ' of ', num2str(numconfigsdirect)]);
            if toppc <=8
                t2 = text(particles*0.35+1,min(quanta+0.5,maxqn+2.5)*0.65,...
                    [num2str(nummicros(currentconfig),'%4.3e'), ' microstates']);
                t2.FontSize = 16;
                t3.FontSize = 16;
            else
                t2 = text(particles*0.30+1,min(quanta+0.5,maxqn+2.5)*0.65,...
                    [num2str(nummicros(currentconfig),'%4.3e'), ' ustates']); 
                t2.FontSize = 12;
                t3.FontSize = 12;
            end
        end
        if (~pf && ~nodispflag)
            
% Plot microstates here
% precalculate matrix of values for all the microstates
%
% If nummicros(i) exceeds 10^4, need to specifically declare it as an
% integer. Seems odd.
            qineachpart = zeros(int32(nummicros(currentconfig)),particles);
%placemat = [];
%cdnum = 1 ;
            microindex = 1;
%    gca = subplot(1,2,2);
%for z=1:nummicros(i) % z enumerates microstates
%global avail_places
% avail_places = ones(nummicros(i),particles);
        % work down from q to q-1 to q-2 to zero quanta
        % Placequanta should create a matrix of the positions for the k'th
        % quanta number and then iteratively drill down for each of those
        % positions. So it needs to be passed a vector of available
        % positions or an index into a list of available positions
%        distrib
%        1:partic
%        zeros(1,length(1:partic))
            [qineachpart, ~, ~, ~] = ...
                placequanta(1,distrib,particles,1:particles,...
                zeros(1,length(1:particles)),qineachpart,microindex);
%qineachpart
% recursive function to place the quanta position k of the distrib vector 
            gcd = subplot(1,2,2);
            particlerow = 1:particles;
            rowq = qineachpart(1,:);
            h = plot(particlerow,rowq,'ro',...
                'MarkerSize',12,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0.75,0,0]);
            axis([0.5 particles+0.5 -.25 quanta+0.5])
            if (~verLessThan('matlab','9.5') && turnoffzoom)
                disableDefaultInteractivity(gca);
            end
            set(gcd,'YTick',0:quanta) 
            set(gcd,'Xtick',0:particles)
            title('All Microstates','FontSize',14)
            xlabel('Particle Number','FontSize',12)
            ylabel('Quanta per particle','FontSize',12)
            h.XDataSource = 'particlerow';
            h.YDataSource = 'rowq';
            t1 = text(particles*0.50,quanta*0.65,['ustate # ',num2str(1)]);
            t1.FontSize = 16;
            t2 = text(2,quanta-0.5,'Spacebar to animate');
            t2.FontSize = 14;
            if (ifig ==1)
disp('If you are staring at this line in the command window and there is no >> below it,');
disp('  find the active plot and hit the spacebar to animate.');
disp('To terminate the animations, Control-C in this command window.');
            end
            pause
            for z = 2:min(num_ustates_stepthrough,nummicros(currentconfig))
% Matlab compiler throws an incorrect error on the next line 
                rowq = qineachpart(z,:); %#ok<NASGU>
                t1.String = '';
                t1.String = ['ustate # ',num2str(z)];
                refreshdata(h,'caller') ;
                drawnow
                pause
                lastz = z;
            end   
            for z = lastz:nummicros(currentconfig)
% Compiler throws incorrect error here too
                rowq = qineachpart(z,:); %#ok<NASGU>
                t1.String = '';
                t1.String = ['ustate # ',num2str(z)];
                refreshdata(h,'caller') ;
                drawnow
            end
        end
    end
%
% This plots a weighted average configuration among the top PC's. Either
% include it in the PC multiplot or display separately
%
if(pf) % We are doing PC's only so make this a subplot
    figure(cfig)
    subf = 1;
    subplot(spx,spy,spavg)
else % We are doing animations so make this a separate figure
    figure
    subf = 0;
end
    if (~verLessThan('matlab','9.5') && turnoffzoom)
    	disableDefaultInteractivity(gca);
    end
    numzeroavg = particles - sum(weightavg)/runsummicro;
    weightavg = weightavg/runsummicro;
    % ifig % for debugging
    bh = barh([quanta+1-find(weightavg ~= 0) 0],[weightavg(weightavg ~= 0) numzeroavg]);
    bh.FaceColor = [0,0.75,0];
    bh.BarWidth = 0.4;
    ylim([-0.3 min(quanta+0.5,maxqn+2.5)]);
    set(gcb,'YTick',0:quanta)
    xl = xlim;
    xlim([xl(1) particles])
    xticks(0:particles)
    yticks(0:quanta)
    % Use the next line to fiddle with right justification for older versions as below
    %tipc2 = title(['Wghtd Avg Pops Among ',num2str(ifig),' PCs'],'FontSize',14);
    title(['Wghtd Avg Pops Among ',num2str(ifig),' PCs'],'FontSize',14);
    % The axes-level alignment commands are for R2020b and later
    if (subf)
        if (verLessThan('matlab','9.9'))
            % Do nothing here -- right justification looks like a lot of
            % fiddling with setting the position, calculating string
            % length, etc.
            %tipc2.HorizontalAlignment = 'right';
        else
            ax1 = gca;
            ax1.TitleHorizontalAlignment = 'right';
        end
    end
    xlabel('Particle Number','FontSize',14)
    ylabel('Quanta per particle','FontSize',14)
    bfit_avg_1 = [num2str(quanta),' quanta dist. among '];
    bfit_avg_2 = [num2str(particles),' particles'];
    %maxqn
    text(.3*particles,.60*min(quanta+0.5,maxqn+1.5),bfit_avg_1,'FontSize',14);
    text(.33*particles,.45*min(quanta+0.5,maxqn+1.5),bfit_avg_2,'FontSize',14);
    if (~verLessThan('matlab','9.5') && turnoffzoom)
    	disableDefaultInteractivity(gca);
    end
    % This allows popping out the PC figure with a button click
    % Do it automatically if toppc is large as the figure is likely to be too small to see
    % (Could do this higher up, when we decide on figure vs subfigure based
    % on pf.)
    if (toppc >= 59)
        hh = copyobj(gca,figure);
%resize the axis to fill the figure
        set(hh, 'Position', get(0, 'DefaultAxesPosition'));
    else
        set(gca,'ButtonDownFcn',@createnew_fig);
    end
%end
%
% This plots the Boltzmann distribution fit
%   Move it over and down a bit to emphasize and show that there might be
%   one behind it
%
    figure('Position',[defx+150 defy-150 defsizex defsizey])
    if (~verLessThan('matlab','9.5') && turnoffzoom)
    	disableDefaultInteractivity(gca);
    end
    hold on
    box on
    theox = 0:.1:quanta+1;
    xdots = 0:quanta;
    bar(xdots+.5,rot90(boltz,2),1,'EdgeColor','c','FaceColor','w','Linewidth',1);
    lb2 = plot(0:quanta,rot90(boltz,2),'bo','MarkerSize',12);
    lb3 = plot(theox,expA.*exp(-theox/T_direct),'b','LineWidth',2);
    lb4 = plot(theox,N0E*exp(-theox/T_calc),'r','LineWidth',2);
%   Line below is for Bose-Einstein statistics and probably wrong  
%    plot(theox,expA_BE.*(1./(exp(theox./T_BE) - 1 )),'m','LineWidth',1)
    title('Development of Boltzmann Distribution Using All Microstates','FontSize',14)
    xlabel('Number of quanta (E)','FontSize',12)
    xlim([0,Inf])
    xticks(0:1:quanta+1)
    ylabel('Number N of particles with energy E','FontSize',12)
    lgdB = legend([lb2,lb3,lb4],'Number N of particles',...
        'Direct Exponential Fit','Continuous model best at high T');
    lgdB.FontSize = 14;
    V = axis;
    bfit_1 = [num2str(quanta),' quanta distributed among ',...
        num2str(particles),' particles']; 
    bfit_2 = [num2str(Wdirect),' total microstates'];
    bfit_2a = [num2str(numconfigscalc),' configurations'];
    bfit_3 = 'Direct enumeration exp fit: ';
    bfit_4 = ['  N = ' num2str(expA,6) ' * exp(-E/',num2str(T_direct),')'];
    bfit_5 = 'Continuous fit: ';
    bfit_6 = ['  N = ' num2str(N0E,6) ' * exp(-E/',num2str(T_calc),')'];

    text(.4*quanta,.76*V(4),bfit_1,'FontSize',14);
    text(.4*quanta,.68*V(4),bfit_2,'FontSize',14);
    text(.4*quanta,.60*V(4),bfit_2a,'FontSize',14);
    text(.4*quanta,.52*V(4),bfit_3,'FontSize',14);
    text(.4*quanta,.44*V(4),bfit_4,'FontSize',14);
    text(.4*quanta,.36*V(4),bfit_5,'FontSize',14);
    text(.4*quanta,.28*V(4),bfit_6,'FontSize',14);
    hold off
end
end

function [W, S] = WS_Einstein_solid(particles,quanta)
% kb is exact as of 2018. K(elvin) is derived
kb = 1.380649e-23 ;
% In case the limits of the factorial function (now 170) change someday or different
% data type is used.
sumlim = 171;
if (particles+quanta) <= sumlim
    W = factorial(particles+quanta-1)/(factorial(quanta)*factorial(particles-1));
    S = kb*log(W);
% The calculation below is can be much slower, but at worst it should return
% Inf and not crash. If an integer becomes much larger than 1/eps then the array
% in the numerator gives elements that are not guaranteed to differ in units of 1,
% so the number of elements in the numerator may not equal quanta, so it crashes.
% [Can't find a reference documenting this specifically but it is easy to demo with
% >> size(2/eps - 1 + 100 - 1 : -1 : 2/eps -1 ) % ]
% If quanta is too large the continued product prod calculation below needs a memory-busting
% array. In that case give up on calculating W. For very large q and p < 50 there
% may be times when this is wrong, when W could be calculated directlt and is actually < realmax
% but the error in back-calculating from S should be small.
% Can push q to 1e8 if you want to for some reason.
elseif ((particles+quanta) < 1/eps && quanta <= 1e7)
    W = prod( ((particles+quanta-1):-1:particles)./(quanta:-1:1));
    if W ~= Inf
        S = kb*log(W);
    else
        S = kb*(gammaln(particles+quanta) - gammaln(quanta+1) - gammaln(particles));
    end
else
%    W = Inf;
    S = kb*(gammaln(particles+quanta) - gammaln(quanta+1) - gammaln(particles));
    W = exp(S/kb);
% else
%     W = Inf;
%     S = kb*(gammaln(particles+quanta) - gammaln(quanta+1) - gammaln(particles));
end
end

function createnew_fig(cb,~)
%cb is the handle of the axes that was clicked
%click on the whitespace within and axes and not on the line object
%copy the axes object to the new figure
hh = copyobj(cb,figure);
%for the new figure assign the ButtonDownFcn to empty
set(hh,'ButtonDownFcn',[]);
%resize the axis to fill the figure
set(hh, 'Position', get(0, 'DefaultAxesPosition'));
end

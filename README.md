# Microstate-anim
Matlab code to simulate a simplified version of the Einstein solid, to illustrate entropy

The guts of the framework is the program micro_movr3, which animates microstates and configurations.
External functions export_fig and partitions from the Matlab File Exchange of useful user-contributed programs are needed.
You may need to create a Matlab account to get these files. I am unsure of any protocol for providing them here,
so I am not doing it.

Then there are several routines that run micro_movr3 and collect its output, for example to build up curves of
entropy vs temperature or to calculate and graph heat capacity.

I also provide a very simple and specialized Matlab tutorial meant to enable students to run and use canned programs,
as opposed to the usual tutorials that lead students into programming Matlab.

Students will need to install the Matlab curvefit toolbox, and they need to know how to get all the programs and
scripts into their matlab Path.

Making this material widely available was supported and inspired by the Elevate Fellows program
run through the Teaching and Learning Transformation Center here at the University of Maryland College Park.

export_fig: http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

partitions: http://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer

From the comments at the top of the micro_movr3 file:

%  micro_movr3 calls JDK's functions placequanta, calcdist, partitionnum, partitionnumrec, 
%     partitionfct, pentagonalnums, and RoundSum -- these must all be in the Matlab Path as well
%   Also requires "partitions" function from the Matlab File Exchange:
%   This one: http://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer
%   Try/catches are there to use cprintf if it is available, to print color to the command window.
%
%   Shows all of the configs and microstates for "quanta" units of energy
%       randomly allocated into "particles" boxes.
%   Animations pause often -- hit spacebar in the plot window to continue, if you hit
%       the spacebar too many times it will bring the front window back and you 
%       may need to bring the plot window to the front again to continue. You
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

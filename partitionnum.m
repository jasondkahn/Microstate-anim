function PNK = partitionnum(z)
% Calculates the partition number of integer z, the number of ways to 
% partition z into a sum of smaller integers. (Can be confused with the
% paryition function, which is the number of ways to partition z into a set
% humber k of terms. [Not the stat mech partition function...]
% This routine is called by partitionfct.m
% This routine calls pentagonalnums.m
% It calls a recursive version of itself, partitionnumrec.m, provided as a
% separate file.
%
% The algorithm uses a recurrence relation
% from https://en.wikipedia.org/wiki/Partition_(number_theory)
% Due to double precision arithmetic in Matlab this does not work forever...
% It works up to z = 296 and then starts to err, is very close up to z = 1000
% By 1200 it is off by 5%, after that it is entirely unreliable, based on
% comparison with http://www.numericana.com/data/partition.htm
% If you really need more precision than 5%, download their list up to
% 4096 and look it up. Beyond that, program smarter.
%
% Needs to calculate pentagonal numbers with values up to z, so the first
% k = sqrt(2/3 of n) of the generalized pentagonal numbers
% Calls pentagonalnums.m to get the list.
% Get a few extra pn's because of the way the routine uses recursion, just
% to be sure we can't run out of them.
    maxk = 2*ceil(sqrt(2*z/3)) + 4;
    pentnums = pentagonalnums(maxk);
if (z < 0)
	PNK = 0;
elseif (z == 0)
    PNK = 1;
elseif (z == 1)
    PNK = 1;
else
    i = 1;
    tz = 0;
    term2 = 0;
    term3 = 0;
    term4 = 0;
    PNK = 0;
    pvec = []; % recursove version builds a list stored in pvec
    while z >= pentnums(i)
        [term1, pnumvec1] = partitionnumrec(z - pentnums(i),pentnums,pvec);
        if ((z - pentnums(i+1)) == 0 )
            tz = 1;
        end
        if ((z - pentnums(i+3)) == 0 || (z - pentnums(i+2)) == 0)
            tz = -1;
        end
        if (z - pentnums(i+1)) > 0
            term2 = pnumvec1(z - pentnums(i+1));
        end
        if (z - pentnums(i+2)) > 0
            term3 = pnumvec1(z - pentnums(i+2));
        end
        if (z - pentnums(i+3)) > 0
            term4 = pnumvec1(z - pentnums(i+3));
        end
        PNK = PNK + tz+ term1 +  term2 - term3 - term4;
        i = i + 4;
        %term1 = 0;
        tz = 0;
        term2 = 0;
        term3 = 0;
        term4 = 0;
    end
end

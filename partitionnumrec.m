function [pnumz, pnumvecnew] = partitionnumrec(z,pentnums,pnumvec)
% Called by partitionnum.m
% This recursive calculation of the partition number of the integer z returns the vector it
% has built along the way, for improved speed.
% Typically this should only be called from itself or partitionnum, but if you
% provide enough pentnums = general pentagonal numbers, this can be run independently
%
if (z < 0)
	pnumz = 0;
    pnumvecnew = pnumvec;
    return
elseif (z == 0)
    pnumz = 1;
    pnumvecnew = pnumvec;
    return
elseif (z == 1)
    pnumz = 1;
    pnumvecnew(1) = 1;
    return
elseif z < length(pnumvec)
    pnumz = pnumvec(z);
    pnumvecnew = pnumvec;
    return
else
	pnumvecnew = pnumvec;
    pnumvecnew(1) = 1;
    i = 1;
    tz = 0;
    %term1 = 0;
    term2 = 0;
    term3 = 0;
    term4 = 0;
    pnumz = 0;
    while z >= pentnums(i)
        %z - pentnums(i)
        %pentnums
        %pnumvec
        %i
        %z - pentnums(i)
        [term1, pnumvec1] = partitionnumrec(z - pentnums(i),pentnums,pnumvecnew);
        pnumvecnew(1:length(pnumvec1)) = pnumvec1;
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
        % vpa works better but takes forever, and is still off by 500
        %pnumz = vpa(pnumz + tz + term1 +  term2 - term3 - term4);
        pnumz = pnumz + tz + term1 +  term2 - term3 - term4;
        tz = 0;
        %term1 = 0;
        term2 = 0;
        term3 = 0;
        term4 = 0;
        i = i + 4;
    end
    %pnumz;
    pnumvecnew(z) = pnumz;
end
end


function [qineachpart, placemat, cdnum, microindex] = placequanta(k,distrib,partic,posvec,qineachvec,qineachpart,microindex)
% Called by micro_mov
% Recursive routine to allocate quanta among particles, with the final
% result returned in qineachpart matrix
% Calls calcdist for making lists of places for identical particles
%global placemat
%global qineachpart
%global microindex
%global cdnum
%global partic
if sum(distrib) > partic  % Should never happen
    disp('placequanta: No room at the inn')
    return
end
if (k > length(distrib) || sum(distrib(k:end)) == 0)
    microindex = microindex+1;
    qineachpart(microindex,:) = qineachvec;
    placemat = [];
    cdnum = 1;
    return
else
    l = distrib(k);  % l is the number of times a q-level shows up in configuration k
% get a matrix of the allowed positions for the l identical quanta #'s among the available slots     
    if l == 0  % No particles with that value of q
%         if k == length(distrib)
%             k = k+1
%             return
%         else

            [qineachpart, placemat, cdnum, microindex] = placequanta(k+1,distrib,partic,posvec,qineachvec,qineachpart,microindex);
%         end
	else
        placemat = [];
%         distrib
%         k
%         l
%         posvec
        cdnum = 1;
        [placemat, cdnum] = calcdist(l,posvec,[],placemat,cdnum);  % returns a matrix of rows with l columns
% in the global placemat variable
% Then need to fill each row with the next quantum down
% Go through each set of positions recursively?
        placetemp = placemat;  % copy it because the next call will rewrite placemat
        qval = length(distrib)-k+1;
        for m = 1:size(placetemp,1)   % row by row
            posvec_new = posvec;
            qineachvec(qineachvec == qval) = 0;
            for n = 1:size(placetemp,2)
                posn = placetemp(m,n);
                qineachvec(posn) = qval;
                % take out the occupied positions
                posvec_new = posvec_new(posvec_new ~= posn);
            end
            if (k == length(distrib) || sum(distrib(k+1:end)) == 0)
%                 microindex
%                 qineachvec
                qineachpart(microindex,:) = qineachvec;
               % qineachvec = zeros(1,length(1:partic));
                microindex = microindex + 1;
            else
%                 k
%                 posvec_new
%                 qineachpart(microindex,:) = qineachvec                
                [qineachpart, placemat, cdnum, microindex] = placequanta(k+1,distrib,partic,posvec_new,qineachvec,qineachpart,microindex);
            end
        end
        % If something is broken consider uncommenting line below. I think
        % it was left behind from the age of global variables but maybe
        % not??
        %qineachvec(qineachvec == qval) = 0;
     %   qineachvec = zeros(1,length(1:partic));
        
    end
end
end


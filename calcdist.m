function [placemat, cdnum] = calcdist(number,position_vec,prev_vec,placemat,cdnum)   % for equivalent q numbers
% Called by placequanta to provide a list (in placemat matrix) of places
% for "number" of identical levels distributed among "position_vec"
% positions. Recursive function.
%global placemat
%global cdnum
%placemat = zeros(factorial(length(position_vec))/factorial(number),number)
if number == 0
    return
end
if number == 1
    for k=1:length(position_vec)
 %       prev_vec
 %       position_vec(k)
        placemat(cdnum,:) = cat(2,prev_vec,position_vec(k));
        cdnum = cdnum+1;
    end
    return
else
    for left = 1:length(position_vec)-number+1
        prev_vec_ext = cat(2,prev_vec,position_vec(left));
%        number - 1
%        position_vec(left+1:end)
        [placemat, cdnum] = calcdist(number-1,position_vec(left+1:end),prev_vec_ext,placemat,cdnum);
    end
end


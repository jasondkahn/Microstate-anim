function [integer_series] = RoundSum(Input_series)
%RoundSum Gives a series of integers that are as close as possible to the
% input series and whose sum is as close as possible to the sum of the input series
% Don't know if this is a unique or even best solution, though
init_sum = sum(Input_series);
len = length(Input_series);
integer_series = floor(Input_series);
tempsum = sum(integer_series);
diff = Input_series - integer_series;
% Sort the differences and return the indexes as well
[dsort,dind] = sort(diff,'descend');
idiff = 1;
% Round up starting with the numbers that are as close as possible to the
% next integer up
if (dsort(end) < init_sum-tempsum) % bail out in case it is fed integers ?
    while ((idiff <= len) && ((init_sum - tempsum) > 0.5))
%        dind(idiff)
        integer_series(dind(idiff)) = integer_series(dind(idiff)) + 1;
        tempsum = tempsum + 1;
        idiff = idiff+1;
    end
    return
else
    return
end


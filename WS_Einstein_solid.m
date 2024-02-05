function [W, S] = WS_Einstein_solid(particles,quanta)
k = 1.380649e-23 ; % apparently exact as of CODATA 2018 recommendations. K is derived
if (particles + quanta) < 170
    W = factorial(particles+quanta-1)/(factorial(quanta)*factorial(particles-1));
    S = k*log(W);
else
    W = Inf;
    S = k*(gammaln(particles+quanta) - gammaln(quanta+1) - gammaln(particles));
end
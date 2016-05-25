%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes the log-PDF of Poisson distribution
%
% Inputs:
% x:   Obserations
% par: Parameters of the PDF (location, scale and shape)
%
% Outputs:
% out: Logarithm of the PDF evaluated at each x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[out] = poisson(x, par)
    part1 = -log( factorial( x ) );
    part2 = x * log( par );
    part3 = -par;
    out = part1 + part2 + part3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
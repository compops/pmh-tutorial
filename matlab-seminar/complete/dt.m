%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes the log-PDF of a non-central Student's t distribution
%
% Inputs:
% x:   Obserations
% par: Parameters of the PDF (degrees of freedome, location, spread)
%
% Outputs:
% out: Logarithm of the PDF evaluated at each x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[out] = dt(x, par)
    nu    = par(1);
    mu    = par(2);
    sigma = par(3);
    
    part1 = gammaln( 0.5 * ( nu + 1 ) );
    part2 = - log(sigma) - 0.5 * log( nu * pi ) - gammaln( 0.5 * nu );
    part3 = -0.5 * ( nu + 1 ) .* log( nu + ( (x - mu) ./ sigma ).^2 );
    part4 =  0.5 * ( nu + 1 ) .* log( nu );

    out = part1 + part2 + part3 + part4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
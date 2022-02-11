function C = c_m(T,par)
% chemical kinetic equation for isotope m (major) subject to multiple reactions

% Syntax
% T: fluid age (or residence time)
% par: Matlab structure with all the parameters

% compute useful terms
Lambda_m=sum(par.lambda); %overall kinetic constant for isotope m
CLIM_m=(par.lambda*par.Clim')/Lambda_m; %overall limiting concentration for isotope m

%overall reaction for isotope m
C = CLIM_m+(par.C0-CLIM_m)*exp(-Lambda_m*T);

end
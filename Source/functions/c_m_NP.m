function C = c_m_NP(T,par)
% chemical kinetic equation for isotope m (major) subject to reaction 1
% only (basically, no mineral precipitation)

% Syntax
% T: fluid age (or residence time)
% par: Matlab structure with all the parameters

% compute useful terms
Lambda_m=par.lambda(1); %overall kinetic constant for isotope m
CLIM_m=par.Clim(1); %overall limiting concentration for isotope m

%overall reaction for isotope m
C = CLIM_m+(par.C0-CLIM_m)*exp(-Lambda_m*T);

end
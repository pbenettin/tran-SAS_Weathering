function C = c_n_FF(T,par)
% chemical kinetic equation for isotope n (rare) subject to multiple
% reactions. Implements a Fixed Fractionation using Rayleigh Distillation

% Syntax
% T: fluid age (or residence time)
% par: Matlab structure with all the parameters

% compute useful terms
Lambda_nf=sum(par.lambda); %overall kinetic constant for isotope n without fractionation
CLIM_nf=(par.lambda*(par.rlim.*par.Clim)')/Lambda_nf; %overall limiting concentration for isotope n without fractionation
Lambda_af=sum(par.alph_f.*par.lambda); %overall kinetic constant for isotope n with fractionation

% little adjustment if size(T) = 1 (need at least 2 values to compute C_n
if length(T) == 1
    T = [T-1,T];
end

% reaction for the major isotope
C_m = c_m(T,par);

% overall reaction for the rare isotope (non fractionating)
C_n_nf=CLIM_nf+(par.C0*par.r0-CLIM_nf)*exp(-Lambda_nf*T); %overall reaction for species n wthout fractionation
C = C_n_nf; %preallocate C equal to non-fractionating reaction

% apply fractionation as soon as precipitation starts
p_start = find(C_n_nf > par.rlim(2)*par.Clim(2),1,'first');
p_start = max(p_start,2); %to avoid errors when p_start = 1
for i = p_start:length(C_m)
    CLIM_rf = (par.alph_f(1)*par.lambda(1)*par.rlim(1)*par.Clim(1)+par.alph_f(2)*par.lambda(2)*(C(i-1)/C_m(i-1))*par.Clim(2))/Lambda_af;
    C(i)=CLIM_rf+(C(i-1)-CLIM_rf).*exp(-Lambda_af*(T(i)-T(i-1)));
end

end
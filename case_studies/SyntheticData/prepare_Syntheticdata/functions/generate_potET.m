function [ET0] = generate_potET(data,par)

% syntax
% ET0=ETmodel(data,par) %choose the particular model type

% variables
% data: table with fields 'time', 'Rg', 'T', 'Hr', 'u'
% par: structure with the parameters of the model


% remember the units:
% 'Rg' = global radiation (W/m2)
% 'T'  = air temperature (degC)
% 'Hr' = relative humidity (%)
% 'u'  = wind speed (m/s)
    

% Define the constants: 
epsil=0.622; % ratio molecular weight of water vapour/dry air [-]
R=0.287; %specific gas constant [kJ/kg/degK]
rho = 1.246; % kg/m3
cp = 1013; % J/kg/degC
gamma = 0.0652; % kPa/degC 
sigma = 5.67e-08; %Stefan-Boltzmann constant [W/(m2 K^4)];
Gsc = 1360; %W/m2
%lambda = 2.45e6; % J/kg
%rs=100/(0.5*24*par.h); %surface resistance (rs=rl/LAIactive) %s/m

% give some defaults if altitude and latitude are not defined
if ~isfield(par,'elev')
    elev = 0;
end
if ~isfield(par,'lat')
    lat = 45;
    lat_rad = lat*pi/180;
end

% Compute other variables based on weather data
lambda = 2.501*10^6-2361*data.T; % J/kg
es = 0.611*exp(17.27*data.T ./ (237.3+data.T)); % vapor pressure at saturation kPa
ea = es.*data.Hr/100; % vapor pressure at saturation kPa
VPD = es-ea; % vapor pressure deficit kPa
Delta = 4098*es ./ (data.T+237.3).^2; %gradient of saturation vapor pressure kPa/degC

% compute net radiation (Rn)
% dayn = day(data.time,'dayofyear');
% P = 101.3 * ((293 - 0.0065 * elev)/293)^5.26;
% gamma = 0.00163 * P/lambda;
% d_r2 = 1 + 0.033 * cos(2 * pi/365 * dayn);
% delta2 = 0.409 * sin(2 * pi/365 * dayn - 1.39);
% w_s = acos(-tan(lat_rad) * tan(delta2));
% R_a = (1440/pi) .* d_r2 .* Gsc .* (w_s .* sin(lat_rad) .* ...
%         sin(delta2) + cos(lat_rad) .* cos(delta2) .* ...
%         sin(w_s)); %W/m2
% R_so = (0.75 + (2 * 10^-5) * elev) * R_a; %W/m2
% f_cloudiness = data.Rg./R_so;
% Rln = sigma * (0.34 - 0.14 * sqrt(ea)) .* (data.T+273.2).^4  .* (1.35 * f_cloudiness - 0.35); %net longwave radiation
% Rsn = data.Rg*(1-par.albedo); %net solar (shortwave) radiation
% Rn = Rsn - Rln; % net radiation W/m2 

% compute net radiation (Rn) VERY SIMPLIFIED
k = 0.8; % reduction factor to take into account the emitted radiation
Rn = k*data.Rg*(1-par.albedo); %simplified approach

% start with the implementations
switch par.ETmodel
        case 'Penman'
        fu=1.313+1.381*data.u; % wind function, Penman 1956 [units are days?]
        ET0 = (Rn.*Delta./lambda + gamma*fu.*VPD/86400) ./ (Delta+gamma)*3600; %kg/m2/h

    case 'PM'
        % Penman-Monteith according to Allen et al. (1998) [i.e. FAO Irrigation and Drainage Paper 56]
        ra=1./(data.u*par.kappa^2) * log((par.z - 0.67*par.h) / (0.123*par.h)) * log((par.z - 0.67*par.h) / (0.1*0.123*par.h)); % aerodynamics resistance s/m
        rs=par.rs;
        ET0 = (Rn.*Delta + rho*cp./ra.*VPD) ./ (lambda.*(Delta+gamma*(1+rs./ra)))*3600; %kg/m2/h

    case 'PMref'
        % Penman-Monteith for the reference crop [i.e. FAO Irrigation and Drainage Paper 56]
        %ra=208./data.u; % aerodynamics resistance s/m
        %rs=70; %surface resistance [s/m]
        %ET0 = (Rn.*Delta + gamma*epsil.*lambda./(1.01*(data.T+273).*R.*ra).*VPD) ./ (lambda.*(Delta+gamma*(1+rs./ra)))*3600; %kg/m2/h
        ET0 = (0.408.*86400/10^6*Rn.*Delta + gamma*900./(data.T+273).*data.u.*VPD) ./ ((Delta+gamma*(1+0.34.*data.u)))/24; %kg/m2/h
end




% include some check
if min(ET0)<-0.1 %this is a potentially serious problem
    warning('min ET0 was < -0.1')
elseif min(ET0)>=-0.1 && min(ET0)<-0.01   %this is considered very minor
    fprintf('\n min ET0 was lower than 0 (%.4f) and was rounded to 0\n',min(ET0))
    ET0(ET0<0)=0;
elseif  min(ET0)>=-0.01 && min(ET0)<0 %so small that we set to zero withouth even displaying a message
        ET0(ET0<0)=0;
end

end
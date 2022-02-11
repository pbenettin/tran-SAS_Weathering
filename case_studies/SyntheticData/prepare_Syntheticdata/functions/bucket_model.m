function data = bucket_model(data,param,dtgen)
% implement a bucket model to make a groundwater balance and compute total
% streamflow from soil and gw contributions

% preallocate
Nt = length(data.J);
Sgw=zeros(Nt,1);  %gw storage [mm]
Qgw=zeros(Nt,1);  %gw flow [mm/d]
Sgw(1)=0.5*param.Sgwmax;   %initial condition for gw storage [-]

% go for the for loop
for i=1:Nt
    
    % groundwater input and output
    Lin=param.frech*data.L;  %leakage input to gw [mm/d]
    Qgw(i) = param.Kgw*(Sgw(i)/param.Sgwmax)^param.cgw;
    
    % groundwater balance (bounded within [0,param.Sgwmax])
    Sgw(i+1) = max(0, Sgw(i) + (Lin(i)-Qgw(i))*dtgen);
    Qexcess = max(0,(Sgw(i+1)-param.Sgwmax)/dtgen);
    Sgw(i+1) = min(Sgw(i+1),param.Sgwmax);
    Qgw(i) = Qgw(i)+Qexcess;
        
end

% add results to the table
data.Q = data.R + (1-param.frech)*data.L + Qgw;
data.f = (data.R+(1-param.frech)*data.L)./data.Q;




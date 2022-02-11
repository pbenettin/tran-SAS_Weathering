function data = sbalance_EF(s0,dtgen,J,E0,T0,param)
% integrate the moisture balance equation using a Euler Forward scheme
% compute all outflows and export to a datatable

% preallocate
Nt = length(J);
R=zeros(Nt,1);  %overland flow [mm/d]
L=zeros(Nt,1);  %leakage [mm/d]
E=zeros(Nt,1);  %actual E [mm/d]
T=zeros(Nt,1);  %actual T [mm/d]
s=zeros(Nt,1);  %soil moisture content [-]
s(1)=s0;        %initial condition for soil moisture [-]

% go for the for loop
for i=1:Nt
    
    % compute fluxes based on timestep i
    E(i)=max(0,min(E0(i),E0(i)*((s(i)-param.sh)/(param.sw-param.sh))));
    T(i)=max(0,min(T0(i),T0(i)*((s(i)-param.sw)/(param.sstar-param.sw))));
    L(i)=max(0,param.Ksat*((s(i)-param.sfc)/(1-param.sfc))^param.c);
    
    % update the balance at time i+1
    s(i+1) = s(i) + (J(i)-E(i)-T(i)-L(i)-R(i))*dtgen/(param.n*param.Z);
    
    % check if there is any overflow
    if s(i+1)>1
        R(i) = (s(i+1)-1)/dtgen*(param.n*param.Z);
        s(i+1) = 1;
    end
end
s(end)=[];

% insert results into a table
data = array2table([s,J,E,T,L,R],'VariableNames',{'s','J','E','T','L','R'});



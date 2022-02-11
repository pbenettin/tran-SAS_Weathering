function [J,dt] = generate_rainfall(Nday,lambda_P,gamma_P,ndown)

% options
dtunits = 1; %if dtunits=1, then units are mm/d, if dtunits=24, then units are mm/h
max_Pdaily = 100; %discard values that are higher than max_Pdaily

% poisson process for the generation of daily rainfall
P_daily=zeros(Nday,1);
for day=1:Nday
    if rand<=lambda_P
        P_daily(day)=min(max_Pdaily,-log(rand)/gamma_P);
    end
end

% random cascade to downscale precipitation from daily to dt data
dt=dtunits/2^ndown;   %time step of the generated series [units depend on dtunits]
subd=2^ndown;
Xrc=zeros(Nday*subd,1); %rainfall depth every d/subd
for i=1:Nday
    if P_daily(i)>0
        app=1;
        for iter=1:ndown
            app2=zeros(2*length(app),1);
            for s=1:length(app)
                %r=rand;
                r = 0.2+rand*(0.8-0.2); %generate random number between 0.2 and 0.8;
                app2(2*(s-1)+1)=app(s)*r;
                app2(2*s)=(1-r)*app(s);
            end
            app=app2;
        end
        Xrc((i-1)*subd+1:i*subd)=app*P_daily(i);
    end
end
J=Xrc/dt; %rainfall intensity at dt time step [units depend on dtunits]

end
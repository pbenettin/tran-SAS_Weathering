% implementation of the age Master Equation using a modified Euler Forward 
% scheme (EF*) that accounts for the presence of 'event' water in streamflow
% and evapotranspiration. 
% Details on the model variables can be found at the bottom of the file

function [varargout] = SAS_EFs_isoweathering(Pars,data)


%--------------------------------------------------------------------------
% prepare the run
%--------------------------------------------------------------------------

% assign CALIBRATION parameters
parQ = Pars(1:data.SASQl); %SASQ parameters
parET = Pars(data.SASQl+1:data.SASQl+data.SASETl); %SASET parameters
S0 = Pars(data.SASQl+data.SASETl+1); %initial storage parameter

% parameters for the MAJOR isotope species (dissolution and precipitation)
par.Clim = Pars([data.SASQl+data.SASETl+2, data.SASQl+data.SASETl+3]);
par.lambda = Pars([data.SASQl+data.SASETl+4, data.SASQl+data.SASETl+5]); % (k/Clim in the paper)

% parameters for the RARE isotope (dissolution and precipitation with fractionation)
% alph = [1, 1]; %when reactions are undisturbed (no fractionation)
par.alph_f = [1, Pars(data.SASQl+data.SASETl+6)]; %when fractionation occurs for the second reaction
par.rlim = Pars([data.SASQl+data.SASETl+7, data.SASQl+data.SASETl+8]);

% set a few constants
NN=length(data.J); %length of the timeseries
ndistr=length(data.index_datesel); %number of age distributions that will be saved

% preallocate variables
S_T=zeros(NN,1); %rank storage (function of age T)
C_Qm=zeros(NN,1); %stream concentration for major isotope (function of time t)
C_Qm_NP=zeros(NN,1); %stream concentration for major isotope without precipitation (function of time t)
C_Qn=zeros(NN,1); %stream concentration for rare isotope with fractionation (function of time t)
r_f=zeros(NN,1); %isotope ratio (fractionating only upon precipitation)
d_f=zeros(NN,1); %isotope delta value (fractionating only upon precipitation)
age_matr=zeros(NN,ndistr); %matrix to store the desired age distributions
med=zeros(NN,1); %median discharge age
Fyw=zeros(NN,1); %young water fraction
Fbp=zeros(NN,1); %fraction of water that does not fractionate (because C<Clim2 for the rare isotope)

% initial conditions
par.C0=0; %rainfall concentration for the main isotope specie
par.r0=par.rlim(1); %rainfall isotope ratio
length_s=1;        %rank storage vector length 
S_T(length_s)=S0;  %initial rank storage [mm]
Omega_Q=feval(data.SASQName,S_T(1:length_s)/S_T(length_s),parQ,data.wi(1)); %[-]
Omega_ET=feval(data.SASETName,S_T(1:length_s)/S_T(length_s),parET,data.wi(1)); %[-]

% automatic: compute useful variables 
Lambda_m=sum(par.lambda); %overall kinetic constant for isotope m
CLIM_m=(par.lambda*par.Clim')/Lambda_m; %overall limiting concentration for isotope m
CLIM_m_NP=par.Clim(1); %overall limiting concentration for isotope m
CMAX_n = max(data.frac_reaction(linspace(0,10/par.lambda(2),1/par.lambda(2)),par)); %maximum concentration value for rare isotope
dMAX = (CMAX_n/CLIM_m-par.rlim(1))/par.rlim(1)*1000; %estimate of maximum delta value

% initial conditions for chemical equations (little importance if the system warms up quickly)
C_Qm(length_s) = CLIM_m;  %initial streamflow concentration major specie
C_Qm_NP(length_s) = CLIM_m_NP;  %initial streamflow concentration major specie
C_Qn(length_s) = CMAX_n;  %initial streamflow concentration rare specie with fractionation

% compute the solution for each species (do it here because they do not
% depend on time but only on age in this formulation)
c_s_m = c_m(data.dt*(0:NN-1)',par); %'batch' concentration (major isotope) of a water parcel
c_s_m_NP = c_m_NP(data.dt*(0:NN-1)',par); %'batch' concentration (major isotope) of a water parcel with no precipitation
c_s_n = data.frac_reaction(data.dt*(0:NN-1)',par); %'batch' concentration (rare isotope) of a water parcel

% find the age at which precipitation starts
index_prec=find(c_s_n>par.rlim(2)*par.Clim(2),1,'first');

%--------------------------------------------------------------------------
% MODEL LOOPS
%--------------------------------------------------------------------------

% let's go
for j=1:NN-1   

    %------------------------------------------------------------------
    % SOLVE THE AGE BALANCE and evaluate the rank storage concentration
    %------------------------------------------------------------------
    % 0) define the domain for the SAS function evaluation (basically a shifted S_T with new water addition)
    age1=max(0,data.dt*(data.J(j)-data.Q(j)*Omega_Q(1)-data.ET(j)*Omega_ET(1))); %estimate of resident water with age 1
    dom=([0;S_T(1:length_s)]+age1)/(S_T(length_s)+age1); %rescaled domain for SAS evaluation

    % 1) evaluate the SAS functions Omega over the domain 'dom'
    Omega_Q=feval(data.SASQName,dom,parQ,data.wi(j)); %[-]
    Omega_ET=feval(data.SASETName,dom,parET,data.wi(j)); %[-]
    
    % 2) solve the master equation balance
    S_T(1:length_s+1)=max(0,[0;S_T(1:length_s)]... 
        +data.dt*data.J(j)... 
        -data.dt*(data.Q(j)*Omega_Q+data.ET(j)*Omega_ET));
    for i=2:length_s+1
        S_T(i)=max(S_T(i),S_T(i-1)); %ensure that S_T is not decreasing
    end 

	
    % 3) compute solute concentration for each parcel in storage    
    % (this is done offline at lines ~60-64 because the reactions only
    % depend on age and not on time)
    
	% 4) make space for the new time step
    length_s=length_s+1;
    
    %----------------------------------------
    % COMPUTE output: stream concentration
    %----------------------------------------
    % compute discharge age distribution (pQ) and concentration (C_Q)
    pQ=diff([0;Omega_Q]);  %[-] this is pQ(T)*dT and it is equivalent to omegaQ(S_T)*dS_T
    C_Qm(j+1)=[c_s_m(1:length_s-1);CLIM_m]'*pQ;   %streamflow modeled concentration
    C_Qm_NP(j+1)=[c_s_m_NP(1:length_s-1);CLIM_m_NP]'*pQ;   %streamflow modeled concentration no precipitation
    C_Qn(j+1)=[c_s_n(1:length_s-1);CMAX_n]'*pQ;   %streamflow modeled concentration

    % compute the isotope ratio and delta values within the stream
    r_f(j+1)=C_Qn(j+1)./C_Qm(j+1); %fractionating only upon precipitation
    d_f(j+1) = (r_f(j+1)-par.rlim(1))/par.rlim(1)*1000; %fractionating only upon precipitation
        
    %-------------------------
    % COMPUTE output: other
    %-------------------------
    % for the selected dates, store discharge age distributions in a matrix
    if any(data.index_datesel-1==j)
        age_matr(1:length_s,data.index_datesel-1==j)=pQ;
    end
    
    % compute some percentile age
    pp=0.5; %percentile [-], (note pp=0.5 is the median)
    med(j+1)=find(Omega_Q>=pp,1,'first');
    
    % compute some young water fractions
    %ywt=30; %young water threshold [days]
    ywt = [...
        14,... %just a small threshold
        75,... %typical 2-3 months as in Kirchner 2016
        1./Lambda_m/24,... %inverse of the effective kinetic constant
        1./Lambda_m/24+index_prec/24*data.dt,... %inverse of the effective kinetic constant + time before starting precipitation
        3*(1./Lambda_m/24+index_prec/24*data.dt)... %just a larger threshold
        ]; %young water threshold(s) [days]
    for i=1:length(ywt)
        Fyw(j+1,i)=Omega_Q(min(length(pQ),round(ywt(i)*24/data.dt)));
    end
    
    % compute fraction of streamflow that carries no fractionation
    Fbp(j+1) = Omega_Q(min(length_s,index_prec));
    
end

%--------------------------------------------------------------------------
% define function output
%--------------------------------------------------------------------------

varargout{1} = []; %just a default

% first check if data.outputchoice exists or not
if isfield(data,'outputchoice') == 1

    % all model values
    if strcmp(data.outputchoice,'C_Qmodel')==1 
        %C_Q(data.Q==0)=NaN;
        varargout{1}=C_Q;
    end

end



%--------------------------------------------------------------------------
% save output
%--------------------------------------------------------------------------

% option to save some output
if isfield(data,'save_output')
    if data.save_output==1
        
        % define the output file name
        outfilename = fullfile(data.case_study_path,'results',...
            data.outfilename);  %set the output filename
        if strcmp(data.outfilename(1:7),'sprintf') %the name is created here
            outfilename = eval(data.outfilename);
        end
        
        % save the variables to output
        varlist={... %list here the variables that you want to save
            'data',...
            'Pars',...
            'Lambda_m','CLIM_m',...
            'CMAX_n','dMAX',...    
            'C_Qm','C_Qm_NP','C_Qn',...
            'r_f',...
            'd_f',...
            'age_matr',...
			'pp','med',...
            'ywt','Fyw',...
            };  
        save(outfilename,varlist{:}) %save selected variables as Matlab file
    end
end



end



%--------------------------------------------------------------------------
% notation details:
% T: age
% t: time
% S: total system storage
% pS: storage age distribution
% pQ: discharge age distribution
% Ps: cumulative age distribution of the system storage
% S_T=S*P_S: rank storage
% Omega: cumulative StorAge Selection function

% all the 'diff' functions return the derivative (df/dx) multiplied by some increment:
% diff(T) represents dT
% diff(Ps) represents (dPs/dT)*dT=pS*dT
% diff(S_T) represents (dS_T/dT)*dT=S*pS*dT
% diff(Omega) represents omega(S_T)*dS_T or equally omega(Ps)*dPs or pQ(T)*dT

% so if one wants the 'classic' variables, some conversion is needed:
% Ps = S_T/S
% omega(S_T) = diff(Omega(S_T))/diff(S_T)
% omega(Ps) = diff(Omega(Ps))/diff(Ps)
% pS(T) = diff(Ps)/diff(T)
% pQ(T) = diff(OmegaQ)/diff(T)
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function associated to a set of simulations

% choose the query among the various possibilities
q_type = 'transport1';

% write various possible queries here
if strcmp(q_type,'targeted_selection1')
    % choose 1-2 specific configurations
    % here a fast and dynamic responding VS a slow and more static responding system
    q = (T.z_Qmin==0.8 & T.z_Qmax == 1 & T.S0 == 300); %small but more static responding
    q = q | (T.z_Qmin==0.2 & T.z_Qmax == 1 & T.S0 == 3000); %large but dynamic
end
if strcmp(q_type,'transport1')
    % transport variability only due to SAS shape (not S0)
    q = true(N,1);
    q = q & T.S0==1000; %condition on transport
    q = q & T.z_Qmax == 0.8; %condition on transport
    q = q & T.chem3 == quantile(unique(T.chem3),.5); %condition on chemistry
end
if strcmp(q_type,'transport2')
    % transport variability only due to SAS shape (not S0)
    q = true(N,1);
    q = q & T.S0==1000; %condition on transport
    q = q & T.z_Qmax == 2; %condition on transport
    q = q & T.chem3 == quantile(unique(T.chem3),.5); %condition on chemistry
end
if strcmp(q_type,'chem1')
    % transport variability only due to SAS shape (not S0)
    q = true(N,1);
    q = q & T.S0 == 1000; %condition on transport
    q = q & T.z_Qmin == 0.5; %condition on transport
    q = q & T.z_Qmax == 2; %condition on transport
end

if sum(q)==0
    error('The query provided no output')
end

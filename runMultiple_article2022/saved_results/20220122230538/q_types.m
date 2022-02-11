% function associated to a set of simulations

% choose the query among the various possibilities
q_type = 'targeted_selection1';

if strcmp(q_type,'targeted_selection1')
    % choose 2-3 specific configurations
    q = (T.z_Qmin==0.3 & T.z_Qmax == 1 & T.S0 == 1000); %reference
    q = q | (T.z_Qmin==0.5 & T.z_Qmax == 0.8 & T.S0 == 300); %smaller, with similar young water contribution but less temporal variability
    q = q | (T.z_Qmin==0.2 & T.z_Qmax == 1.2 & T.S0 == 3000); %larger, with similar young water contribution but more temporal variability
end

if strcmp(q_type,'targeted_selection2')
    % choose 2-3 specific configurations
    q = (T.z_Qmin==0.3 & T.z_Qmax == 1 & T.S0 == 1000); %reference
    q = q | (T.z_Qmin==0.2 & T.z_Qmax == 1.2 & T.S0 == 300); %smaller, with similar young water contribution but more temporal variability
    q = q | (T.z_Qmin==0.5 & T.z_Qmax == 0.8 & T.S0 == 3000); %larger, with similar young water contribution but less temporal variability
end

if sum(q)==0
    error('The query provided no output')
end


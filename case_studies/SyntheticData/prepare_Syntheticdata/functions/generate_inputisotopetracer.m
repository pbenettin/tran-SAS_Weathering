function data = generate_inputisotopetracer(data,param)
% generation of virtual rainfall isotope composition: create a sinusoid and
% add noise to it

% generate the basic sinusoid for d18O
per = 1*365.25; %yearly period [d]
Cin0 = param.m18O+param.amp18O * sin(2*pi*(data.time-param.phase18O)/per);

% generate correlated normal errors
Nt = size(data,1);
err=zeros(Nt,1); %for d18O
err(1)=param.noiseAmplitude*randn;
for i=2:Nt
    err(i) = param.noisecorr*err(i-1)+...
        sqrt(1-param.noisecorr^2)*param.noiseAmplitude*randn; %for d18O
end

% combine the sinuoid with the error to get the final input timeseries
Cin_d18O = Cin0 + err;
Cin_d2H = Cin_d18O*param.lmwl1 + param.lmwl2 + 1*randn(Nt,1); %add some additional gaussian noise to d2H

% assign the output
data.d18O = Cin_d18O;
data.d2H = Cin_d2H;

end
function [actresp,extrapresp,kdiff]=extrapresponse(st,kvec,amat,algprms)
% Extrapolate the response of the community to change in conditions

spnum=length(st);      % number of species
rndwidth = algprms(5); % width of distribution for perturbation of k
exrptime = algprms(6); % length of time to use for extrapolation
exrpsnum = algprms(7); % number of frames to use for extrapolation

kdiff=(rand(size(kvec))-0.5)*rndwidth; % random change in carrying capacities enacted

tmpkvec=kvec+kdiff; % change carrying capacities

% run simulations, for short-term and long term response (same scenario, different time points taken)
shortst=get_ode_sol(tmpkvec,amat,[-1 exrptime exrptime/exrpsnum 0 st(:)']);
longst =get_ode_sol(tmpkvec,amat,[algprms(1:4) st(:)']);

% fit linear curve of short-term response for each species separately
xx=(0:exrpsnum)'*exrptime/exrpsnum;
for ii=1:spnum
    yy=shortst(:,ii);
    tmp=polyfit(xx,yy,1); % fitting a linear function
    outvec(ii)=tmp(1); % taking the slope of a linear fit
end;

actresp=longst-st;
extrapresp=outvec./(kvec.*st);

end

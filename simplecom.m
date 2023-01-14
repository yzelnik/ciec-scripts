function [kvec,amat,state]=simplecom(spnum,Ks,As,rndseed,algprms,settozero)
% Create a GLV community, following Barbier et al. 2018 PNAS 
% spnum is the number of species (not all of which are necessarily extant at equlibirum)
% Ks=[kavg kwdth] defines the carrying capacities distributions' average and standard-deviation
% As=[mu sigma gamma] defines the non-self interactions, as their average, standard deviation, and symmetry around diagonal.
% rndseed is used for randomization, and defines the specific instance of the community
% algprms is a vector of numerical parameters used for the ODE solver, see more info inside function get_ode_sol()
% settozero>0 sets a fraction of interaction matrix values to zero

% Allow for some parameters to not be defined
if(nargin<6) settozero=0; end;
if(length(As)<3) As = [As(:)' 0]; end;
if(length(Ks)<3) Ks = [Ks(:)' 0]; end;

% Randomize
rng(rndseed);

% Pick random values from a uniform distribution for carrying capacities
kvec = Ks(1)+(rand(1,spnum)-0.5)*Ks(2)*sqrt(12);

% Start with a gaussian matrix with zeros at the diagonal
tmpmat = randn(spnum);
tmpmat(logical(eye(spnum)))=0;

if(As(3)) % if there is symmetry/asymmetry (i.e., gamma is not 0)
    tmpmat = tmpmat+((1-sqrt(1-As(3)^2))/As(3))*tmpmat';
end;

% Normalize sigma to 1
tmpmat = tmpmat/std(tmpmat(:)); 
linknum= spnum*(1-settozero);
% Add in the effect of mu and sigma
amat = tmpmat*As(2)/sqrt(linknum)+As(1)/linknum; 
% If some interactions are to be set to zero
if(settozero) 
    % choose settozero fraction 
    tmpmat = triu(rand(spnum)<settozero); 
    % set both in and out going links to zero
    amat(logical(tmpmat+tmpmat')) = 0; 
end
% Set main diagonal (self-interaction) to -1
amat(logical(eye(spnum)))=-1;

% Find steady-state solution using the ode45 solver (for a GLV model)
state=get_ode_sol(kvec,amat,algprms);

end
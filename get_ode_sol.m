function state=get_ode_sol(kvec,amat,algprms)
% define some algorithmic parameters
thresh   = algprms(1);
maxtime  = algprms(2);
basetime = algprms(3);
warnflag = algprms(4);
% define start vector
if(length(algprms)<=5)
   start = repmat(algprms(5),size(kvec));
else
   start = reshape(algprms(5:end),size(kvec));
end;

if(thresh>0) % normal run, with a threshold
  maxchange = inf;
  % keep going while convergence is not good and we are not at max time
  while(maxchange/basetime>thresh && basetime*2<maxtime)
	% run ODE integration
	[~,Y] = ode45(@ode_glv,(0:4)*basetime/4,start,[],kvec,amat);
	% calculate change in last leg of simulation
    if(size(Y,1)<5)
        basetime=inf;
    else
        maxchange = max(abs(Y(4,:)-Y(5,:)));
        basetime=basetime*2;
        start = Y(5,:); % ready for next iteration
    end;
  end;

  % check if we converged
  if(maxchange/basetime>thresh) && warnflag
	warning('Solution does not appear to converge.')
  end;
  % Test if we got all the data we got (otherwise, this is a bad simulation)
  if(size(Y,1)<5)
	  state(1,:) = NaN;
  else
	  state = Y(5,:);
  end;
else % get frames within the predefined timeframe, not trying to reach equilibrium
    [~,state] = ode45(@ode_glv,0:basetime:maxtime,start,[],kvec,amat);
end;

end

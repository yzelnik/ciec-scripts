%% Assemble communities

randnum = 100; % number of randomized communities for a given value of sigma
spnum   = 50;  % number of (initial) species in the community
zerofrac= 0.8; % fraction of interactions set to zero
algprms = [1e-6 1e4 1e2 1 1]; % algorithmic parameters for the ODE solver
kwdth   = 0.2; % standard deviation of carrying capacities

sigvals=(1:50)/50; % define vector of different sigma values

alivethresh=1e-5;  % threshold for extant species
shanres = 0.001; % resolution for binning when calculating shannon diversity

for sigind=1:length(sigvals) % go over different values of sigma (y)
  for randind=1:randnum % go over different randomizations
    
    % save sigma value
    cursig=(sigvals(sigind)); 
    % assemble community (get back k-vector, interaction matrix, and abundnaces)
    [kvec,amat,st] = simplecom(spnum,[1 kwdth 0],[-1*cursig cursig -1],randind*100+sigind,algprms,zerofrac);
    
    % save all these for later
    allvec(:,sigind,randind)=kvec;
    allmat(:,:,sigind,randind)=amat;
    allsts(:,sigind,randind)=st;
    
    % also, which species are extant?
    alivenum(sigind,randind)=sum(st>alivethresh);
    
  end;
end;

%% calculate different measures on the set of communities

minspnum=20;

for sigind=1:length(sigvals) % go over different values of sigma (y)
   tic;
   for randind=1:randnum % go over different randomizations
      
      alive= allsts(:,sigind,randind)>alivethresh; % true/false vector for extant species
      amat = allmat(alive,alive,sigind,randind); % interaction matrix
      kvec = allvec(alive,sigind,randind)';      % k-vector
      st   = allsts(alive,sigind,randind)';      % abundnace vector
      actsp= length(st); % number of extant species in the community
  
        
      % Some general measures
      collec(sigind,randind) = max(abs(eig(amat+eye(actsp)))); % collectivity
      %actgam(sigind,randind) = calcgammatrix(amat);  % actual gamma value
      lnstab(sigind,randind) = max(real(eig(amat))); % linear stability

      % biotic contribution 
      relyld = st./kvec; % relative yield
      ryfo   = mean((relyld-1).^2); % mean distance of: relative yield from one
      ryfz   = mean(relyld.^2);     % mean distance of: relative yield from zero
      abdisc(sigind,randind) = ryfo./ryfz;
    
      % temporal unpredictability
      [bclong1,bcshort1]=extrapresponse(st,kvec,amat,[algprms(1:4) 0.5 0.5 10]); % extra prms: rnd-wdth of k-change, time-for extrap, point-num for extrap
      tmpcor   = corrcoef(bcshort1,bclong1);
      teminc(sigind,randind) = 1-tmpcor(1,2);
      
      % perturbation (removal) depth
      [res,extra]=removal_depth(st,kvec,amat,algprms);
      pdepth(sigind,randind)=res(1);
      p3d(sigind,randind,:)=res(3:5);
      p6d(sigind,randind,:)=res(6:8); % second normalization option, of sum over sum
      
      % shanon-diversity in direct and inverse matrix
      invmat = inv(amat);
      shan1(sigind,randind)=shanindex(amat(~eye(sum(alive))),shanres);
      shan2(sigind,randind)=shanindex(invmat(~eye(sum(alive))),shanres);
      shnrat(sigind,randind)=(shan2(sigind,randind)./shan1(sigind,randind));
      
    end;
    toc;
end;
save('community_data')

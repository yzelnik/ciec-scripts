%% assemble some communities and calculate some measures for them
algprms = [1e-6 1e4 1e2 1 1];  % algorithmic parameters for the ODE solver
commnum = 200; % number of communities
kwdth   = 0.2; % width of K distribution
alivethresh = 1e-5;  % threshold for extant species

comminitsz=randi(11,commnum,1)+4; % randomly choose community size between 5 and 15

% initilize
ind=0;
intmats={}; stvals={}; kkvals={}; alivenum=[];
collec=[]; complx=[]; lnstab=[];
% run through each community
for ii=1:commnum
    % three parameters for constructing the A matrix
    aprmvec=rand(1,3).*[-2 2 2]+[0 0.1 -1];
    % assemble a community
    [kvec,amat,st] = simplecom(comminitsz(ii),[1 kwdth 0],aprmvec,ii,algprms);
    % save info
    alive=st>alivethresh;
    intmats{ii}=amat(alive,alive);
    stvals{ii}=st(alive);
    kkvals{ii}=kvec(alive);
    alivenum(ii)=sum(alive);
    
    % if community is large enough, calculate collectivity, complexity, and stability
    if(alivenum(ii)>2)
      ind=ind+1;
      collec(ind)=max(abs(eig(intmats{ii}+eye(alivenum(ii))))); % collectivity
      complx(ind)=sqrt(sum(sum((intmats{ii}+eye(alivenum(ii))).^2)))./sqrt(alivenum(ii)); % complexity
      lnstab(ind)=max(real(eig(intmats{ii}))); % linear stability
    end;
    if(mod(ii,20)==0) disp(ii); end; % print out progress
end;

%% plot out results comparing collectivity with complexity and stability

% plot collectivity vs. complexity (in 2 versions of the latter)
subplot(1,2,1)
plot(log10(collec),log10(complx),'*',log10(collec),log10(complx.*sqrt(alivenum(alivenum>2)-1)),'o',[-1 1],[-1 1],'k','lineWidth',1.5)
xlabel('log10 collectivity \phi','fontSize',20);
ylabel('log10 complexity','fontSize',20);
legend('C','C*(S-1)^{0.5}','location','northwest')

% plot collectivity vs. stability
subplot(1,2,2)
scatter(collec,-lnstab,30,complx,'lineWidth',1.5)
hold on;
plot([0 1],[1 0],'k','lineWidth',2) % add line where collectivity cannot cross
hold off;
colorbar
xlabel('collectivity \phi','fontSize',20);
ylabel('stability','fontSize',20);


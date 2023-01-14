%% load sims
load community_data

% Below is the code to create the figures for the four main signatures (Figs. 2-5)
% Each one should be run separately

%% effective connectance (shannon-ratio)

% set up figure arrangment 
clf;
ha(1)=axes('Units','normalized', 'Position',[0.07 0.12 0.48 0.86]);
ha(2)=axes('Units','normalized', 'Position',[0.64 0.73 0.35 0.25]);
ha(3)=axes('Units','normalized', 'Position',[0.64 0.42 0.35 0.25]);
ha(4)=axes('Units','normalized', 'Position',[0.64 0.12 0.35 0.25]);

% define 3 specific points
sigchs  = [3 16 49]; % values along the y control parameter axis
randchs = [1 2 1];   % which randomization to choose

for jj=1:3 % calculate location of specific points in main panel
     tmpxys(jj,:)=[collec(sigchs(jj),randchs(jj)), shnrat(sigchs(jj),randchs(jj))];
end;

% plot effective connectance vs. collectivity
axes(ha(1))
set(gca, 'ColorOrder', jet(size(collec,1)), 'NextPlot', 'replacechildren');
plot(collec',shnrat','.','markerSize',8)
axis([0 2 1 4])
set(gca,'fontsize',14)
% plot specific points
hold on;
plot(tmpxys(:,1),tmpxys(:,2),'.k','markerSize',45);
hold off;
% write panel letters next to points
for jj=1:3
    text(tmpxys(jj,1)+0.07,tmpxys(jj,2)+0.05,char('a'+4-jj),'fontSize',24);
end;
text(0.05,3.9,'a','fontSize',36)

% labels and such
xlabel('collectivity','fontSize',24)
ylabel('effective connectance increase','fontSize',24)
set(gca,'xTick',0:1:3,'yTick',0:5)

% plot other panels
xlms=[ -0.1 0.1 ; -0.5 0.5; -1. 1.]*1.05;
for jj=1:3
    curjmp=diff(xlms(jj,:))/40;
    axes(ha(5-jj))
    % get the extant species, the interaction matrix, and its inverse
    alive= allsts(:,sigchs(jj),randchs(jj))>alivethresh;
    amat = allmat(alive,alive,sigchs(jj),randchs(jj)); 
    invmat =inv(amat);
    
    % get off-diagonal cell values of the direct and inverse matrix
    tmphist=[];
    tmphist(:,1)=amat(~eye(length(amat))); 
    tmphist(:,2)=invmat(~eye(length(amat))); 
    
    % plot historgram of direct and net interactions
    histogram((tmphist(:,1)),-2:curjmp:2,'FaceAlpha',0.5); 
    hold on; 
    histogram((tmphist(:,2)),-2:curjmp:2,'FaceAlpha',0.5,'FaceColor','r'); 
    hold off;
    set(gca,'fontsize',14)
    set(gca,'YScale','log')
    text(xlms(jj,1)*0.98+xlms(jj,2)*0.02,800,char('a'+(4-jj)),'fontSize',36)
    ylim([0 3000])
    xlim(xlms(jj,:))
end;

% labels and such
axes(ha(4))
xlabel('interaction strength','fontSize',24)
axes(ha(3))
ylabel('histogram','fontSize',24)
text(0.15,500,'direct int.','backgroundcolor',[0.35 0.6 1.0],'fontSize',18)
text(0.2,50,'net int.','backgroundcolor',[1.0 0.5 0.5],'fontSize',18)

%% figure for perturbation depth

% set up figure arrangment 
clf;
ha(1)=axes('Units','normalized', 'Position',[0.07 0.12 0.48 0.86]);
ha(2)=axes('Units','normalized', 'Position',[0.64 0.12 0.35 0.86]);

% plot perturbation depth vs. collectivity
axes(ha(1))
set(gca, 'ColorOrder', jet(size(collec,1)), 'NextPlot', 'replacechildren');
plot(collec',pdepth','.','markerSize',8)
set(gca,'fontsize',14)
xlim([0 2])

% labels and such
xlabel('collectivity','fontSize',24)
ylabel('perturbation depth','fontSize',24)
set(gca,'xTick',0:1:5,'yTick',0:0.2:2)
text(0.06,1.67,'a','fontSize',36)

% plot second panel
tmpcolors=[0. 0. 0.; 0.4 0.4 0.4; 0.8 0.8 0.8];%0.8 0.3 0.8; 0.4
axes(ha(2))
hold on;
for ii=1:3 % dummy plot (just to get legend right)
    plot(0,0,'.','color',tmpcolors(ii,:),'markerSize',20);
end;
for ii=1:3 % plot points of each grey shade (type/distance of link) 
    plot(collec,p6d(:,:,ii),'.','color',tmpcolors(ii,:))
end;
hold off;

% legends and such
set(gca,'fontsize',14)
xlim([0 2])
xlabel('collectivity','fontSize',24)
ylabel('relative impact','fontSize',24)
set(gca,'xTick',0:1:5,'yTick',0:0.2:10)
ll=legend('on direct links','on secondary links','on other links');
set(ll,'fontSize',20);
text(0.08,1,'b','fontSize',36)
ylim([0 1.05])



%% figure for temporal unpredictability

% set up figure arrangment 
clf;
ha(1)=axes('Units','normalized', 'Position',[0.07 0.12 0.48 0.86]);
ha(2)=axes('Units','normalized', 'Position',[0.64 0.57 0.35 0.41]);
ha(3)=axes('Units','normalized', 'Position',[0.64 0.12 0.35 0.41]);

% define 2 specific points
sigchs  = [6 50]; % values along the y control parameter axis
randchs = [7 15]; % which randomization to choose
largest = 6;
pertstr = 0.2;
frames  = 0:0.1:4;

% go over specific communities (in 2 specific points)
for jj=1:length(sigchs)
    % get interaction matrix, carrying capacity vector & abundance vector
    alive= allsts(:,sigchs(jj),randchs(jj))>alivethresh;
    amat = allmat(alive,alive,sigchs(jj),randchs(jj));
    kvec = allvec(alive,sigchs(jj),randchs(jj))';
    st   = allsts(alive,sigchs(jj),randchs(jj))';
   
    % save for later
    allk{jj}=kvec;
    alln{jj}=st;
    
     % define peturbation in k vector
    rng(jj);
    kshift{jj}=(rand(size(kvec))-0.5)*pertstr;
    % run simulations to see the dyamics unfold
    [~,out] = ode45(@ode_glv,frames,st,[],kvec+kshift{jj},amat);
    pertfrms{jj} = out;
    
    % calculate location of specific points in main panel
    tmpxys(jj,:)=[collec(sigchs(jj),randchs(jj)), teminc(sigchs(jj),randchs(jj))];
end;

% plot temporal unpredictability vs. collectivity
axes(ha(1))
set(gca, 'ColorOrder', jet(size(collec,1)), 'NextPlot', 'replacechildren');
plot(collec',teminc','.','markerSize',8)
axis([0 2 0 0.4])
set(gca,'fontsize',14)
% plot specific points
hold on;
plot(tmpxys(:,1),tmpxys(:,2),'.k','markerSize',45);
hold off;
% write panel letters next to points
for jj=1:length(sigchs)
    text(tmpxys(jj,1)+0.15,tmpxys(jj,2)-0.015,char('a'+1+length(sigchs)-jj),'fontSize',24);
end;
text(0.05,0.285,'a','fontSize',36)

% labels and such
ylim([0 0.3])
xlabel('collectivity','fontSize',24)
ylabel('temporal unpredictability','fontSize',24)
set(gca,'xTick',0:1:5,'yTick',0:0.2:2)

% setup for other panels
ylms=[-0.11 0.11; -0.11 0.11; -0.11 0.11];
yjmps=[0.02;0.1;0.1];
clrs = [0.9 0.2 0; 0.1 0.7 0.2; 0 0.2 0.9; 0.7 0.2 0.7; 0.2 0.7 1; 0.9 0.8 0.2; 0.2 0.2 0.2; 0.7 0.7 0.7; 0.2 0.5 0.9; 0.9 0.5 0.2];

% plot other panels
for jj=1:length(sigchs)
    % choose which species to show
    tmpdata = pertfrms{jj};
    [~,sortind]=sort(abs(diff(tmpdata([1 end],:))),'descend');
    spchs = sortind(1:largest);
    
    axes(ha(2+length(sigchs)-jj))
    if(jj>1) set(gca,'xTickLabel',[]); end;
    
    % get shift in k vector and abundances vector (saved earlier)
    finstate = kshift{jj};
    curri    = alln{jj};
    
    % plot curves of each small panel
    hold on;
    for kk=1:length(spchs) % for each species, plot the actual trajectory and the extrapolated one
      plot(frames,tmpdata(:,spchs(kk))-repmat(tmpdata(1,spchs(kk)),length(frames),1),'lineWidth',2,'color',clrs(kk,:));%,frames([1 end]),[1;1]*finstate(spchs),'k:','lineWidth',2);
      plot(frames,(1-exp(-frames*curri(spchs(kk))))*finstate(spchs(kk)),':','lineWidth',2,'color',clrs(kk,:));
    end;
    hold off;
    % finish up each panel
    set(gca,'xTick',0:2:10,'yTick',-1:yjmps(jj):2)
    set(gca,'fontsize',14)
    text(0.12,ylms(jj,1)*0.10+ylms(jj,2)*0.85,char('a'+(1+length(sigchs)-jj)),'fontSize',36)
    ylim(ylms(jj,:));
    box on;
    set(gca,'yTick',(-1:1)*0.05)
end;
axes(ha(3))
xlabel('time','fontSize',24)
axes(ha(3))
ylabel('                            abundance change','fontSize',24);


%% figure for biotic contribution to species realized niches

% set up figure arrangment 
clf;
ha(1)=axes('Units','normalized', 'Position',[0.07 0.12 0.48 0.86]);
ha(2)=axes('Units','normalized', 'Position',[0.64 0.71 0.17 0.27]);
ha(3)=axes('Units','normalized', 'Position',[0.64 0.415 0.17 0.27]);
ha(4)=axes('Units','normalized', 'Position',[0.64 0.12 0.17 0.27]);
ha(5)=axes('Units','normalized', 'Position',[0.82 0.71 0.17 0.27]);
ha(6)=axes('Units','normalized', 'Position',[0.82 0.415 0.17 0.27]);
ha(7)=axes('Units','normalized', 'Position',[0.82 0.12 0.17 0.27]);

% define 3 specific points
sigchs  = [2 22 48]; % values along the y control parameter axis
randchs = [1 5 2];   % which randomization to choose

for jj=1:3 % calculate location of specific points in main panel
     tmpxys(jj,:)=[collec(sigchs(jj),randchs(jj)), abdisc(sigchs(jj),randchs(jj))];
end;

% plot biotic contribution vs. collectivity
axes(ha(1))
set(gca, 'ColorOrder', jet(size(collec,1)), 'NextPlot', 'replacechildren');
plot(collec',abdisc','.','markerSize',8)
set(gca,'fontsize',14)
axis([0 2 0 1.3])
% plot specific points
hold on;
plot(tmpxys(:,1),tmpxys(:,2),'.k','markerSize',45);
hold off;
% write panel letters next to points
for jj=1:3
    text(tmpxys(jj,1)-0.05,tmpxys(jj,2)+0.08,char('a'+4-jj),'fontSize',24);
end;
text(0.05,1.24,'a','fontSize',36)

% labels and such
xlabel('collectivity','fontSize',24)
ylabel('biotic contribution','fontSize',24)
set(gca,'xTick',0:1:3,'yTick',0:0.5:2)

% plot other panels
clrs = [0.9 0.2 0; 0.1 0.7 0.2; 0 0.2 0.9; 0.7 0.2 0.7; 0.2 0.7 1; 0.9 0.8 0.2];
for jj=1:3
    axes(ha(5-jj))
    % get the interaction matrix, k vector, and abundance vector
    alive= allsts(:,sigchs(jj),randchs(jj))>alivethresh;
    amat = allmat(alive,alive,sigchs(jj),randchs(jj));
    kvec = allvec(alive,sigchs(jj),randchs(jj))';
    st   = allsts(alive,sigchs(jj),randchs(jj))';
    
    % plot kvector vs. abundance
    plot(kvec,st,'b*',[0 2],[0 2],'k--','lineWidth',1)
    set(gca,'xTick',0.:0.5:2,'yTick',0.:0.5:2)
    if(jj>1) set(gca,'xTickLabel',[]); end;
    text(0.05,1.15,char('a'+(4-jj)),'fontSize',36)
    set(gca,'fontsize',14)
    axis([0.0 1.4 0. 1.4])
    
    % plot the effect of interactions (sum of matrix: (1-A)^-1)
    axes(ha(8-jj))
    intereffect=transpose(sum(inv(amat),2));
    plot(-intereffect,st,'r*',[0 2],[0 2],'k--','lineWidth',1)
    set(gca,'xTick',0.:0.5:2,'yTick',0.:0.5:2,'yTickLabel',[])
    if(jj>1) set(gca,'xTickLabel',[]); end;
    set(gca,'fontsize',14)
    axis([0.0 1.4 0. 1.4])
end;
axes(ha(7))
xlabel('$\sum_j (1-A)^{-1}_{ij}$','fontSize',24,'Interpreter','latex')
axes(ha(4))
xlabel('$K_i$','fontSize',24,'Interpreter','latex')
axes(ha(3))
ylabel('equlibirum biomass','fontSize',24);


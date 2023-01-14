function [res,extra]=removal_depth(st,kvec,amat,algprms)
% Calculate the effect of removeing a species in the community

spnum=length(st);

% run over each species in the community
for ii=1:spnum
    start=st;
    start(ii)=0; % remove species ii
    % find the effect on other species
    tmpst=get_ode_sol(kvec,amat,[algprms(1:4) start(:)']);
    if(length(tmpst)<length(start))
        tmpst=nan(size(start));
    end;
    stin(:,ii) =start;
    stout(:,ii)=tmpst;
end;

% calculate the distances between species
for ii=1:spnum
    dists(:,ii)=marknetdist(amat,ii)'-1;
end;
% removal-depth, average over different removals: SUM(d*|in-out|)/SUM(|in-out|)
res = mean(sum(dists.*abs(stin-stout),2)./sum(abs(stin-stout),2));

% calculate some proxy (using inverse of interaction matrix)
bmat = abs(inv(amat));
bmat(logical(eye(spnum)))=0;
res(2) = mean(sum(dists.*bmat./repmat(sum(bmat,1),spnum,1),2));

% normalization -- maximal change
normval = max(abs(stin(:)-stout(:)));
% save 3 groups (direct links, secondary links, all other links)
ints{1}=(dists==1);
ints{2}=(dists==2);
ints{3}=(dists>2);

ints{4}=(dists>0); % calculate a different normalization
normval2 = sum(abs(stin(ints{4})-stout(ints{4})));

for ii=1:3 % for the 3 groups, save distance using two different normalizations
    res(2+ii) = mean(abs(stin(ints{ii})-stout(ints{ii})))/normval;
    res(5+ii) = sum(abs(stin(ints{ii})-stout(ints{ii})))/normval2; % sum divided by sum
end;

% return more info
extra = {stin,stout,dists};

end

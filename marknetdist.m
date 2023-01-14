function dists=marknetdist(mat,start)
% calculate the network distances from an initial node

% initilize
dists    = zeros(size(mat,1),1);
curnodes = false(size(mat,1),1);
curnodes(start)=1; % starting location within the network
binarymat= 0+logical(mat);

ind=1;
newnodes = curnodes; % first step
dists(newnodes)=ind;
while(sum(newnodes)) % run iteratively while we find new nodes
    ind=ind+1;
    actnodes = logical(binarymat*curnodes); % which nodes are reachable from our currently explored nodes?
    newnodes = (actnodes & (~curnodes));    % which are the newly found nodes?
    curnodes(newnodes)=1;                   % update which nodes we have gone through
    dists(newnodes)=ind;  % mark the distance (from the initial node) of these new nodes
end

end

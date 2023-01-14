function res=shanindex(data,datares)
% calculate shannon index

if(nargin<2) datares=0; end;
data=data(:);
if(datares) % keep value of cells only up to datares?
    data=round(data/datares)*datares;
end;
unq=unique(data);

% go over each unique value
for ii=1:length(unq)
    fracval=mean(data==unq(ii));
    parts(ii)=fracval.*log(fracval);
end;

% return the sum of all the differnet parts
res=-sum(parts);
end

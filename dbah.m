function [dbz] = dbah(mat)

% GIANMARCO PINTON
% CREATED: 2019-11-22
% LAST MODIFIED: 2019-11-22
% db scale with max at zero
% absolute value 
% zero padded hilbert

dbz=hilbert(mat,size(mat,1)*2);
dbz=abs(dbz);
dbz=db(dbz);
dbz=dbz-maxmax(dbz);


dim=ndims(mat);
if(dim==1)
    dbz=dbz(1:size(dbz,1)/2);
end
if(dim==2)
    dbz=dbz(1:size(dbz,1)/2,:);
end
if(dim==3)
    dbz=dbz(1:size(dbz,1)/2,:,:);
end
if(dim==4)
    dbz=dbz(1:size(dbz,1)/2,:,:,:);
end
if(dim>=5)
    disp('ERROR dimension too large')
end

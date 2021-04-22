function [idx]=ovalIdx(dims,cen,latrad,axialrad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: SEPT 29, 2020 (Oleksii Ostras)
% find indices of a 2D oval
% do not go out of bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx=zeros(1,ceil(latrad+1)^2); cc=1; %ceil rounds each element of X
%to the nearest integer greater than or equal to that element
if(length(dims)==2)
  for i=-ceil(latrad):ceil(latrad)
    for j=-ceil(axialrad):ceil(axialrad)
      if ((i^2)/(latrad^2))+((j^2)/(axialrad^2))<=1
        %idx(cc)=round(i+cen(1))*dims(1)+round(j+cen(2));
        idx(cc)=round(i+cen(1))+round(j+cen(2))*dims(1);
        cc=cc+1;
      end
    end
  end
else
  Disp('ERROR! length(dim)~=2!!')
end
idx=idx(1:cc-1);
idx=idx(find(idx>0));
idx=idx(find(idx<=dims(1)*dims(2)));
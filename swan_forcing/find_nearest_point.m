function [xip, etap]=find_nearest_point(varargin)
%
% [ix, iy]=find_nearest_point(arg1, arg2, arg3, ag4, arg5, arg6)
% where ix=xip and iy=etap
%
% list of arguments:
% arg1 and 2: x and y of the point (optional)
% arg3: name of X array
% arg4: name of Y array
% arg5: name of mask array
% arg6: option for automatic, non-interactively calls the script
% e.g.
% for known coordinates:
% [ix, iy]=find_nearest_point(-80.99, 32, lon_rho, lat_rho, mask_rho, 'auto')
% to pick a point on the figure:
% [ix, iy]=find_nearest_point(lon_rho, lat_rho, mask_rho)
%
% Zafer Defne 2015

if length(varargin{1})>1 %&& ~(nargin>5 &&  strcmpi(varargin(end),'auto')) % coordinates are picked interactively
    if size(varargin,2)>2
        if ~strcmpi(varargin{end}, 'auto')
            fprintf('Click on the image to select the point location\n');
            mask=varargin{end};
        else
        end 
        
    end
    aa=gca;
    X=varargin{1};
    Y=varargin{2};
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');
    x = point1(1,1);
    y = point1(1,2);
%     fprintf('(x,y)= (%5.4f,%5.4f)\n',x,y);
else % coordinates are given
    x=varargin{1};
    y=varargin{2};
        if ~strcmpi(varargin{end}, 'auto')
            fprintf('Click on the image to select the point location\n');
            mask=varargin{end};
        else
        end
        
    X=varargin{3};
    Y=varargin{4};
end

d=sqrt((Y-y).^2+(X-x).^2);
d=d(:);
% sort according to the distance
[dsort, ind]=sort(d, 'ascend');
% indp is the index of the closest point
ind1=1;
indp=ind(ind1);

if exist('mask', 'var')
    while (mask(indp)==0)
        ind1=ind1+1;
        indp=ind(ind1);
    end
end
if ~( strcmpi(varargin(end),'auto'))
    try
        pcolor(X,Y, mask), axis equal, shading flat
        hold(aa,'on')
    catch
        fprintf('Cannot plot the mask.\n')
    end
    
    
    try
        hold(aa,'on')
    catch
        try
            aa= gca;
            hold(aa,'on')
        catch
            disp('no active figure');
        end
    end
    
    
    plot(X(indp),Y(indp), 'xk',X(indp),Y(indp), '+r',X(indp),Y(indp), 'ow')
end

lx=size(X,2);
ly=size(Y,1);

%% remapping the indices for the 2D array from a column vector
ix=ceil(indp/ly);
if(mod(indp,ly)==0)
    iy=ly;
else
    iy=mod(indp,ly);
end
[etap xip]=ind2sub([ly lx], indp);
if ~(etap == iy && xip ==ix)
    disp('Indices do not match. Check formulation')
end
if ~(strcmpi(varargin(end),'auto'))
    fprintf('Point index (etap, xip)= (%d,%d)\n', etap, xip);
end
end


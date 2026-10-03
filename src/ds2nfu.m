function varargout=ds2nfu(varargin)
% DS2NFU Convert data coordinates to normalized figure coordinates.
%
% Usage:
%
%   [Xf,Yf]=ds2nfu(X,Y)
%   [Xf,Yf]=ds2nfu(hAx,X,Y)
%   POSf=ds2nfu(POS)
%   POSf=ds2nfu(hAx,POS)
%
% The returned coordinates are normalized to the parent figure.

narginchk(1,3);

% Determine whether the first input is an axes handle.
if isgraphics(varargin{1},'axes')
    hAx=varargin{1};
    args=varargin(2:end);
else
    hAx=gca;
    args=varargin;
end

if isempty(args)
    error('ds2nfu:InvalidInput',...
        'Coordinate input is missing.');
end

if numel(args)==1
    pos=args{1};

    if numel(pos)~=4
        error('ds2nfu:InvalidPosition',...
            'Position must contain four elements.');
    end

    pos=double(pos(:).');
    isPosition=true;

elseif numel(args)==2
    x=args{1};
    y=args{2};

    if ~isequal(size(x),size(y))
        error('ds2nfu:SizeMismatch',...
            'X and Y must have the same size.');
    end

    x=double(x);
    y=double(y);
    isPosition=false;

else
    error('ds2nfu:InvalidInput',...
        'Use X,Y or a four-element position vector.');
end

% Store original units.
oldAxesUnits=get(hAx,'Units');
oldFigureUnits=get(ancestor(hAx,'figure'),'Units');

% Use normalized units for the axes position.
set(hAx,'Units','normalized');
axPos=get(hAx,'Position');

% Get current axis limits.
axLim=axis(hAx);

xLim=axLim(1:2);
yLim=axLim(3:4);

xWidth=diff(xLim);
yHeight=diff(yLim);

if xWidth==0 || yHeight==0
    error('ds2nfu:InvalidAxesLimits',...
        'Axes limits must have nonzero width and height.');
end

if isPosition
    pos(1)=axPos(1)+(pos(1)-xLim(1))*axPos(3)/xWidth;
    pos(2)=axPos(2)+(pos(2)-yLim(1))*axPos(4)/yHeight;
    pos(3)=pos(3)*axPos(3)/xWidth;
    pos(4)=pos(4)*axPos(4)/yHeight;

    varargout{1}=pos;
else
    Xf=axPos(1)+(x-xLim(1))*axPos(3)/xWidth;
    Yf=axPos(2)+(y-yLim(1))*axPos(4)/yHeight;

    varargout{1}=Xf;
    varargout{2}=Yf;
end

% Restore units.
set(hAx,'Units',oldAxesUnits);
set(ancestor(hAx,'figure'),'Units',oldFigureUnits);

end


% function varargout = ds2nfu(varargin)
% % DS2NFU  Convert data space units into normalized figure units. 
% %
% % [Xf, Yf] = DS2NFU(X, Y) converts X,Y coordinates from
% % data space to normalized figure units, using the current axes.  This is
% % useful as input for ANNOTATION.  
% %
% % POSf = DS2NFU(POS) converts 4-element position vector, POS from
% % data space to normalized figure units, using the current axes.  The
% % position vector has the form [Xo Yo Width Height], as defined here:
% %
% %      web(['jar:file:D:/Applications/MATLAB/R2006a/help/techdoc/' ...
% %           'help.jar!/creating_plots/axes_pr4.html'], '-helpbrowser')
% %
% % [Xf, Yf] = DS2NFU(HAX, X, Y) converts X,Y coordinates from
% % data space to normalized figure units, on specified axes HAX.  
% %
% % POSf = DS2NFU(HAX, POS) converts 4-element position vector, POS from
% % data space to normalized figure units, using the current axes. 
% % Michelle Hirsch
% % mhirsch@mathworks.com
% % Copyright 2006-2014 The MathWorks, Inc
% 
% %% Process inputs
% narginchk(1, 3)
% 
% % Determine if axes handle is specified
% if length(varargin{1})== 1 && ishandle(varargin{1}) && strcmp(get(varargin{1},'type'),'axes')	
% 	hAx = varargin{1};
% 	varargin = varargin(2:end);
% else
% 	hAx = gca;
% end
% errmsg = ['Invalid input.  Coordinates must be specified as 1 four-element \n' ...
% 	'position vector or 2 equal length (x,y) vectors.'];
% % Proceed with remaining inputs
% if length(varargin)==1	% Must be 4 elt POS vector
% 	pos = varargin{1};
% 	if length(pos) ~=4
% 		error(errmsg);
%     end
% else
% 	[x,y] = deal(varargin{:});
% 	if length(x) ~= length(y)
% 		error(errmsg)
% 	end
% end
% 
% 	
% %% Get limits
% axun = get(hAx,'Units');
% set(hAx,'Units','normalized');
% axpos = get(hAx,'Position');
% axlim = axis(hAx);
% axwidth = diff(axlim(1:2));
% axheight = diff(axlim(3:4));
% 
% 
% %% Transform data
% if exist('x','var')
% 	varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
% 	varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
% else
% 	pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
% 	pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
% 	pos(3) = pos(3)*axpos(3)/axwidth;
% 	pos(4) = pos(4)*axpos(4)/axheight;
% 	varargout{1} = pos;
% end
% 
% 
% %% Restore axes units
% set(hAx,'Units',axun)
% 

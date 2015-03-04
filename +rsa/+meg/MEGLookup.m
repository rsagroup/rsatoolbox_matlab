function [varargout] = MEGLookup(varargin)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

% [g1 g2 m] = MEGLookup(c, n)
%
% Given a channel index, c, this will return the labels of the first and second grads and of the mag.  Here, n is the number of sensor types being used.
%
% 	-~- OR -~-
%
% c = MEGLookup(descrip, i)
%
% Given a channel, i, and a channel description, descrip, (which is a string; either 'g1', 'g2' or 'm'), this will return c, a channel location index

if nargin == 2

	if ~ischar(varargin{1})
	
		c = varargin{1};
		n = varargin{2};
	
		if c > 102 || c < 1, error('c must be in [1, 102] \cap N'); end%if

		for i = 1:n

			varargout{i} = n * c - (n - i);

		end%for:i
	
	else
	
		descrip = varargin{1};
		cIn = varargin{2};
	
		switch descrip
			case 'g1'
				cOut = (cIn + 2) / 3;
			case 'g2'
				cOut = (cIn + 1) / 3;
			case 'm'
				cOut = cIn / 3;
			otherwise
				error('Only ''g1'', ''g2'' and ''m'' are accepted as channel names.')
		end%switch:descrip
	
		if cOut ~= floor(cOut)
			error([num2str(cIn) ' is not a channel of type ' descrip '.']);
		end%if
	
		varargout{1} = cOut;

	end%if

else

	error(['Wrong number of arguments: ' num2str(nargin) ' is not 2.']);

end%if

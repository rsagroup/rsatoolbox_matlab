function RDMs = averageRDMs_subjectSession(varargin)
% RDMs = averageRDMs_subjectSession(RDMs): returns the input unchanged (no
% operation). RDMs = averageRDMs_subjectSession(RDMs, 'subject'): computes
% the subject-averaged RDMs. RDMs = averageRDMs_subjectSession(RDMs,
% 'session') RDMs = averageRDMs_subjectSession(RDMs, 'subject', 'session')
% RDMs = averageRDMs_subjectSession(RDMs, 'session', 'subject')
%
% Averages a struct of RDMs which is [nMasks, nSubjects, nSessions] based
% on arguments and returns the averaged struct. It must be noted that RDMs
% should be given in a particular format so that this function can be
% aplied for averaging. The first dimension would be brain region (e.g.
% hIT), the second dimension would be the subjects (e.g. subject5) and the
% third would be the recording session (e.g. session1).  
% The naming convention is the
% following:
% RDMs(1).name = 'RoiName | subjectName | sessionName/number';
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

	if nargin == 1 || nargin == 2 || nargin == 3

		RDMs = varargin{1};

		if nargin == 2

			if strcmpi(varargin{2}, 'subject')
				RDMs = aSub(RDMs);
			elseif strcmpi(varargin{2}, 'session')
				RDMs = aSes(RDMs);
			else
				wrongStringsError();
			end%if

		elseif nargin == 3

			if strcmpi(varargin{2}, 'subject')
				RDMs = aSub(RDMs);
			elseif strcmpi(varargin{2}, 'session')
				RDMs = aSes(RDMs);
			else
				wrongStringsError();
			end%if

			if strcmpi(varargin{3}, 'subject') && ~strcmpi(varargin{2}, 'subject')
				RDMs = aSub(RDMs);
			elseif strcmpi(varargin{3}, 'session') && ~strcmpi(varargin{2}, 'session')
				RDMs = aSes(RDMs);
			else
				wrongStringsError();
			end%if

		end%if

	else

		wrongNargsError();

	end%if

end%function

%% Subfunctions %%

function wrongStingsError
	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
	error('Please put ''subject'' or ''session'' in for the string arguments.');
end%function

function wrongNargsError
	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
	error('Only accepts 1, 2 or 3 arguments.');
end%function

function aveRDMs = aSub(RDMs)
	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
	nMa = size(RDMs, 1);
	nSu = size(RDMs, 2);
	nSe = size(RDMs, 3);
	for ma = 1:nMa
		for se = 1:nSe
			for su = 1:nSu
				if su == 1
					suSum = RDMs(ma, su, se).RDM;
				else
					suSum = suSum + RDMs(ma, su, se).RDM;
				end%if:su==1
			end%for:su
			aveRDMs(ma, 1, se).RDM = suSum ./ nSu;
			aveRDMs(ma, 1, se).color = RDMs(ma, 1, se).color;
			oldName = RDMs(ma, 1, se).name;
			bothBits = findstr(oldName, ' | ');
			if numel(bothBits) == 2
				firstBit = bothBits(1);
				secondBit = bothBits(2);
				newName = [oldName(1:firstBit-1) oldName(secondBit:length(oldName))];
			elseif numel(bothBits) == 1
				firstBit = bothBits(1);
				newName = oldName(1:firstBit-1);
			end%if
			aveRDMs(ma, 1, se).name = newName;
		end%for:se
	end%for:ma
end%function

function aveRDMs = aSes(RDMs)
	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
	nMa = size(RDMs, 1);
	nSu = size(RDMs, 2);
	nSe = size(RDMs, 3);
	for ma = 1:nMa
		for su = 1:nSu
			for se = 1:nSe
				if se == 1
					seSum = RDMs(ma, su, se).RDM;
				else
					seSum = seSum + RDMs(ma, su, se).RDM;
				end%if:se==1
			end%for:se
			aveRDMs(ma, su, 1).RDM = seSum ./ nSe;
			aveRDMs(ma, su, 1).color = RDMs(ma, 1, se).color;
			oldName = RDMs(ma, su, 1).name;
			bothBits = findstr(oldName, ' | ');
			if numel(bothBits) == 2
				firstBit = bothBits(1);
				secondBit = bothBits(2);
				newName = oldName(1:secondBit-1);
			elseif numel(bothBits) == 1
				firstBit = bothBits(1);
				newName = oldName(1:firstBit-1);
			end%if
			aveRDMs(ma,su, 1).name = newName;
		end%for:su
	end%for:ma
end%function
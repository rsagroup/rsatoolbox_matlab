function makeGifti(file,varargin)
% function rsa_makeGifti(file,varargin)
% Convert file into the GIFTI file format and write it into
% the same directory. Aligns surface with the volume anatomical for the
% subject.
%
% INPUTS
%   file: full path and file name for the file to be converted 
%
% USEROPTIONS / VARARGIN:
%	type: surface (default)
%
% EXAMPLE:
%   % Convert a surface file to gifti
%   opt.type = 'surf';
%   rsa_makeGifti('lh.white','type','surf');
%
% Naveed Ejaz
% n.ejaz@ucl.ac.uk
% 2/2015

import rsa.util.*

Opt.type = 'surf';
Opt      = rsa.getUserOptions(varargin,Opt); 

% Get file path and name
[fpath,fname,ext] = fileparts(file);

switch(Opt.type)
    case 'surf'
        [A1,A2] = readAffineTransMatrix(fpath);
        A       = A2/A1;
        [Vw,F]  = read_surf(file);
        [data(:,1),data(:,2),data(:,3)] = rsa.affine_transform(Vw(:,1),Vw(:,2),Vw(:,3),A);        

        g.vertices  = data;
        g.faces     = F+1;          % 1 based indexing
        g           = gifti(g);
        g.mat(1:3,4)  = zeros(3,1);

        rsa_save_gii(g,fullfile(fpath,[fname ext '.surf.gii']));    % File name must have 
                                                                    % .surf postfix to work in caret     
    otherwise
        disp('Incorrect file type ...');
end;

function [A1,A2] = readAffineTransMatrix(fpath)
% Reads the affine transformation matrices
%   - vox2ras-tkr
%   - vox2ras
anaFile = fullfile(fpath,'../mri/brain.mgz');

[status,result] = system(['mri_info ' anaFile ' --vox2ras-tkr']);
A               = sscanf(result,'%f');
A1              = reshape(A,4,4)';

[status,result] = system(['mri_info ' anaFile ' --vox2ras']);
A               = sscanf(result,'%f');
A2              = reshape(A,4,4)';

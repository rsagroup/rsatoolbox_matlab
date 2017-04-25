function T = readSurf(white,pial)
% function rsa_readSurf(file,varargin)
% Reads surfaces for each subject and returns an array of structures with
% the corresponding surfaces. Only GIFTI input files formats are supported
% for now.
%
% INPUTS
%   white: Nx1 cell array of white surfaces (left or right sided)
%   pial:  Nx1 cell array of corresponding pial surfaces
% OUTPUT:
%   T: Nx1 cell array of structures containing the following fields
%           - white: vertices belonging to the white surface
%           - pial: vertices belonging to the corresponding pial surface
%           - topo: topology for the c1/c2 surface
% EXAMPLE:
%   % Convert a surface file to gifti
%   white = {'lh.white.surf.gii','rh.white.surf.gii'};
%   pial  = {'lh.pial.surf.gii','rh.pial.surf.gii'};
%   T     = rsa_readSurf(white,pial);
%
% Naveed Ejaz
% n.ejaz@ucl.ac.uk
% 2/2015

import rsa.util.*

nSurf = length(white);
for i=1:nSurf
    [fpath,fname,ext] = fileparts(white{i});
    
    switch(ext)
        case '.gii'
            g1 = gifti(white{i});
            g2 = gifti(pial{i});
            
            T{i}.topo  = g1.faces';
            T{i}.white = g1.vertices';
            T{i}.pial  = g2.vertices';
        otherwise
            disp('Unrecognized file extension');
    end;    
end;

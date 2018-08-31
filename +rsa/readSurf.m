function T = readSurf(white,pial,topo)
% function rsa_readSurf(file,varargin)
% Reads surfaces for each subject and returns an array of structures with
% the corresponding surfaces. Only GIFTI input files formats are supported
% for now.
%
% INPUTS
%   white: Nx1 cell array of white surfaces (left or right sided)
%   pial:  Nx1 cell array of corresponding pial surfaces
%   topo:  Nx1 cell array of corresponding toplogies (if not included in
%           the white/pial file. 
% OUTPUT:
%   T: Nx1 cell array of structures containing the following fields
%           - white: vertices belonging to the white surface
%           - pial: vertices belonging to the corresponding pial surface
%           - topo: topology for the c1/c2 surface
% EXAMPLE 1:
%   % Read a series of gifti files as surfaces 
%   white = {'lh.white.surf.gii','rh.white.surf.gii'};
%   pial  = {'lh.pial.surf.gii','rh.pial.surf.gii'};
%   T     = rsa_readSurf(white,pial);
%
% EXAMPLE 2:
%   % Read a series of gifti files as surfaces 
%   white = {'lh.white.coord','rh.white.coord'};
%   pial  = {'lh.pial.coord','rh.pial.coord'};
%   topo  = {'lh.CLOSED.topo','rh.CLOSED.topo'};
%   T     = rsa_readSurf(white,pial);
%
% Naveed Ejaz
% n.ejaz@ucl.ac.uk
% 2/2015
% V2 03/2016- added '.coord' case
% V3 10/2017- moved from +fmri to base rsa folder

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
        case '.coord' 
            g1 = caret_load(white{i});
            g2 = caret_load(pial{i});
            g3 = caret_load(topo{i});
            
            T(i).topo  = g3.data';
            T(i).white = g1.data';
            T(i).pial  = g2.data';
        otherwise
            disp('Unrecognized file extension');
    end;    
end;

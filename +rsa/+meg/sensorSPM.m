% based on http://imaging.mrc-cbu.cam.ac.uk/meg/SensorSpm
% to be run in SPM 5 EEG
% runs group statistics on files in "imgfiles"
% output written into directory "outdir"
% OH, March 2009

clear all

addpath /imaging/local/meg_misc;    % for meg_batch_anova

% output directory
outdir = [pwd, '/stats'];
% root directory common to all input fiff-files (not required if full pathnames specified in imgfiles{})
root_dir = [pwd, '/spm/'];

% Here: One group of subjects, 2 conditions
% You can also specify more groups of subjects, or test one condition against zero:
% More groups of subjects: Vary first cell index of imgfiles (imgfiles{1}{cc}...imgfiles{2}{cc}...)
% One-sample t-test against zero mean: just specify one file per line,
% e.g. imgfiles{1}{cc} = ['/thispath/subj1_con1.img']; cc=cc+1; etc.
if exist('imgfiles')~=1,
    cc = 1;
    
    imgfiles{1}{cc} = 'meg08_0319_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0320_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0323_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0324_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0327_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0348_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0350_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0363_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0366_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0371_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0372_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0377_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0380_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0397_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0400_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0401_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
    imgfiles{1}{cc} = 'meg08_0402_RDM-stem_i_SPM-mags/trialtype1/saverage.img'; 
%     
%     imgfiles{1}{cc} = 'meg08_0319_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0320_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0323_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0324_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0327_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0348_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0350_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0363_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0366_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0371_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0372_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0377_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0380_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0397_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0400_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0401_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; cc=cc+1;
%     imgfiles{1}{cc} = 'meg08_0402_RDM-suffix_i_SPM-mags/trialtype1/saverage.img'; 
end;
nr_subjects = cc;

% Attach root directory
clear Panova;
for i=1:nr_subjects, % create full image file names, including root directory etc.
    [m,n] = size(imgfiles{1}{i});
    for j=1:m,
        Panova{1}{i}(j,:) = fullfile( root_dir, imgfiles{1}{i}(j,:) );
    end;    %..j
end;    %..i

meg_batch_anova(Panova,outdir);  % Run SPM stats

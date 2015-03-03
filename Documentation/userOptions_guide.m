%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.maskPath
%                        A string describing the path to the location of the
%                        files for the definition of RoI masks.
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in fullBrainVols.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Defaults to the fieldnames of the first subject
%                        of binaryMasks_nS. If working in MEG sensor space,
%                        (right now) you should either have a single entry
%                        called 'masked', or leave this unset.
%                userOptions.betaPath
%                        A string which contains the absolute path to the
%                        location of the beta images. It can contain the
%                        following wildcards which would be replaced as
%                        indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%                                [[betaIdentifier]]
%                                        To be replaced by filenames as provided
%                                        by betaCorrespondence.
%                userOptions.conditionLabels
%                        A cell array containing the names of the conditions in
%                        this experiment.
%                userOptions.voxelSize
%                        A tripple consisting of the [x y z] dimensions of each
%                        voxel in mm.
%                userOptions.structuralsPath
%                        A string which contains the absolute path to the
%                        location of the structural images and the normalisation
%                        warp definition file. It can contain the following
%                        wildcards which would be replaced as indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%                userOptions.distanceMeasure
%                        A string descriptive of the distance measure to be used
%                        to compare two RDMs. Defaults to 'Spearman'.
%                userOptions.saveFigurePDF
%                        A boolean value. If true, the figure is saved as a PDF.
%                        Defaults to false.
%                userOptions.saveFigurePS
%                        A boolean value. If true, the figure is saved as a PS.
%                        Defaults to false.
%                userOptions.saveFigureFig
%                        A boolean value. If true, the figure is saved as a
%                        MATLAB .fig file. Defaults to false.
%                userOptions.saveFiguresJpg
%                        A boolean value. If true, the figure is saved as a jpeg
%                        image.
%                userOptions.displayFigures
%                        A boolean value. If true, the figure remains open after
%                        it is created. Defaults to true.
%                userOptions.imagelables
%                        Defaults to empty (no image labels).
%                userOptions.rankTransform
%                        Boolean value. If true, values in each RDM are
%                        separately transformed to lie uniformly in [0,1] with
%                        ranks preserved. If false, true values are represented.
%                        Defaults to true. This is for display only.
%                userOptions.criterion
%                        The criterion which will be minimised to optimise the
%                        MDS arrangement. Defaults to metric stress.
%                userOptions.rubberbands
%                        Boolean value. If true, rubberbands indicating MDS
%                        distortion are drawn on the MDS plot. Defaults to true.
%                userOptions.distance
%                        A string indicating the distance measure with which to
%                        calculate the RDMs. Defaults to 'Correlation'.
%                userOptions.RoIColor
%                        A triple indicating the [R G B] value of the colour
%                        which should be used to indicated data RDMs on various
%                        diagrams. Defaults to black ([0 0 0]).
%                userOptions.ModelColor
%                        A triple indicating the [R G B] value of the colour
%                        which should be used to indicated model RDMs on various
%                        diagrams. Defaults to black ([0 0 0]).
%                userOptions.nResamplings
%                        How many bootstrap resamplings should be performed?
%                        Defaults to 1000.
%                userOptions.resampleSubjects
%                        Boolean. If true, subjects will be bootstrap resampled.
%                        Defaults to false.
%                userOptions.resampleConditions
%                        Boolean. If true, conditions will be resampled.
%                        Defaults to true.             
%                userOptions.dotSize
%                        A number specifying the size of the circles used
%                        to indicate points in MDS plots. Defaults to 8,
%                        probably.
%                userOptions.significanceTestPermutations
%                        An integer which describes the number of random
%                        permutations to be used to calculate significance.
%                        Defaults to 10,000.
%                userOptions.conditionColours
%                        A [nConditions 3]-sized matrix indicating an [R G B]
%                        tripple colour for each condition. Defaults to all
%                        black.
%                userOptions.convexHulls
%                        A vector of length equal to the number of conditions.
%                        Each entry in the vector corresponds to the same-index
%                        condition, and the number for this entry represents a
%                        category for this condition. Convex hulls are drawn
%                        around all conditions of the same category on the MDS
%                        plots, coloured by the first colour in
%                        userOptions.conditionColours for the points in this
%                        ategory. For example:
%                                [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]
%                        Would represent four categories, each with four
%                        conditions. If unset, convex hulls will not be drawn.
%                userOptions.colourScheme
%                        A colour scheme for the RDMs. Defualts to jet(64).
%                userOptions.distanceMeasure
%                        A string descriptive of the distance measure to be used
%                        to compare two RDMs. Defaults to 'Spearman'.
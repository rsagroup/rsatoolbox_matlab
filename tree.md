.
├── Demos
│   ├── 92imageData
│   │   ├── 92_behavRDMs.mat
│   │   ├── 92_brainRDMs.mat
│   │   ├── 92_modelRDMs.mat
│   │   ├── faceAnimateInaniClustersRDM.mat
│   │   ├── Kriegeskorte_Neuron2008_supplementalData.mat
│   │   ├── rdm92_HMAXnatImPatch.mat
│   │   ├── rdm92_V1model.mat
│   │   └── simTruePatterns.mat
│   ├── anatomy.mat
│   ├── compareRefRDM2candRDMs_barGraph.pdf
│   ├── compareRefRDM2candRDMs_RDMcomparisonPvalues.pdf
│   ├── defineUserOptions2.m
│   ├── defineUserOptions.m
│   ├── DEMO1_RSA_ROI_simulatedAndRealData.m
│   ├── DEMO2_RSA_ROI_sim.m
│   ├── DEMO3_LDt_sim.m
│   ├── DEMO4_RSAsearchlight_sim.m
│   ├── modelRDMs_demo2.m
│   ├── modelRDMs_SL_sim.m
│   ├── projectOptions_DEMO1.m
│   ├── projectOptions_demo.m
│   ├── sampleMask_org.mat
│   ├── simulationOptions_demo_LDt.m
│   ├── simulationOptions_demo.m
│   └── simulationOptions_demo_SL.m
├── Documentation
│   ├── license.txt
│   ├── toolbox documentation.pdf
│   └── userOptions_guide.m
├── README.md
├── Recipes
│   ├── betaCorrespondence.m
│   ├── defineUserOptions.m
│   ├── Recipe_fMRI.m
│   └── Recipe_fMRI_searchlight.m
├── +rsa
│   ├── compareRefRDM2candRDMs.m
│   ├── constructModelRDMs.m
│   ├── constructRDMs.m
│   ├── crossvalIPM.m
│   ├── defineSearchlight.m
│   ├── defineSearchlight_surface.m
│   ├── defineSearchlight_volume.m
│   ├── dendrogramConditions.m
│   ├── distanceLDC.m
│   ├── +fig
│   │   ├── addComparisonBars.m
│   │   ├── addHeading.m
│   │   ├── addImageSequenceToAxes.m
│   │   ├── colorScale.m
│   │   ├── exportCurrentFigAsPDF.m
│   │   ├── exportCurrentFigAsPostscript.m
│   │   ├── exportfig.m
│   │   ├── figureDendrogram.m
│   │   ├── figureMDSArrangement.m
│   │   ├── handleCurrentFigure.m
│   │   ├── imageRDMs.m
│   │   ├── image_thr.m
│   │   ├── mat2RGBimage.m
│   │   ├── pageFigure.m
│   │   ├── paneling.m
│   │   ├── plotDotsWithTextLabels.m
│   │   ├── plotTextLabels.m
│   │   ├── RDMcolormap.m
│   │   ├── rubberbandGraphPlot.m
│   │   ├── selectPlot.m
│   │   ├── setPapertoFigPos.m
│   │   ├── shepardPlot.m
│   │   ├── showRDMs.m
│   │   ├── showVol.m
│   │   ├── showVoxObj.m
│   │   └── xticklabel_rotate.m
│   ├── figureRDMs.m
│   ├── +fmri
│   │   ├── addBinaryMapToVol.m
│   │   ├── addRoiToVol.m
│   │   ├── boyntonModel.m
│   │   ├── fMRIDataMasking.m
│   │   ├── fMRIDataPreparation.m
│   │   ├── fMRIMaskPreparation.m
│   │   ├── fMRISearchlight.m
│   │   ├── generateCognitiveModel_fastButTrialsNeedToStartOnVols.m
│   │   ├── readMask.m
│   │   ├── readSurf.m
│   │   ├── searchlightMapping_fMRI.m
│   │   ├── spatiallySmooth4DfMRI_mm.m
│   │   ├── sphericalRelativeRoi.m
│   │   └── temporallySmoothTimeSpaceMatrix.m
│   ├── getUserOptions.m
│   ├── @gifti
│   │   ├── Contents.m
│   │   ├── display.m
│   │   ├── export.m
│   │   ├── fieldnames.m
│   │   ├── gifti.m
│   │   ├── isfield.m
│   │   ├── plot.m
│   │   ├── private
│   │   │   ├── base64decode.m
│   │   │   ├── base64encode.m
│   │   │   ├── dunzip.m
│   │   │   ├── dzip.m
│   │   │   ├── getdict.m
│   │   │   ├── isintent.m
│   │   │   └── read_gifti_file.m
│   │   ├── save.m
│   │   ├── struct.m
│   │   ├── subsasgn.m
│   │   └── subsref.m
│   ├── +gifti
│   │   ├── makeGifti.m
│   │   └── save_gii.m
│   ├── MDSConditions.m
│   ├── MDSRDMs.m
│   ├── pairwiseCorrelateRDMs.m
│   ├── +rdm
│   │   ├── averageRDMs_subjectSession.m
│   │   ├── categoricalRDM.m
│   │   ├── concatenateRDMs.m
│   │   ├── concatRDMs.m
│   │   ├── concatRDMs_unwrapped.m
│   │   ├── interleaveRDMs.m
│   │   ├── squareRDM.m
│   │   ├── squareRDMs.m
│   │   ├── stripNsquareRDMs.m
│   │   ├── unwrapRDMs.m
│   │   ├── vectorizeRDM.m
│   │   ├── vectorizeRDMs.m
│   │   ├── wrapAndNameRDMs.m
│   │   └── wrapRDMs.m
│   ├── runSearchlightLDC.m
│   ├── runSearchlight.m
│   ├── +sim
│   │   ├── generateBetaPatterns.m
│   │   ├── simulateClusteredfMRIData_fullBrain.m
│   │   ├── simulateClusteredfMRIData.m
│   │   └── simulateDataFiles.m
│   ├── +spm
│   │   ├── crossvalIPMraw.m
│   │   ├── distanceLDCraw.m
│   │   ├── getDataFromSPM.m
│   │   ├── noiseNormalizeBeta.m
│   │   ├── spm_create_vol.m
│   │   ├── spm_type.m
│   │   ├── spm_write_plane.m
│   │   └── spm_write_vol.m
│   ├── +stat
│   │   ├── bootstrapRDMs.m
│   │   ├── ceilingAvgRDMcorr.m
│   │   ├── compactPvalueString.m
│   │   ├── covdiag.m
│   │   ├── FDRthreshold.m
│   │   ├── features2ModelRDMs.m
│   │   ├── fisherDiscrTRDM.m
│   │   ├── fisherTransform.m
│   │   ├── fitModelOLS.m
│   │   ├── gaussian_nk.m
│   │   ├── modelMANOVA.m
│   │   ├── pairCACt.c
│   │   ├── pairCACt.m
│   │   ├── pairCACt.mexmaci64
│   │   ├── raeSpearmanCorr.m
│   │   ├── randomlyPermute.m
│   │   ├── randomPermutation.m
│   │   ├── rankCorr_Kendall_taua.m
│   │   ├── RDMCorrMat.m
│   │   ├── signrank_onesided.m
│   │   └── varianceLDC.m
│   └── +util
│       ├── affine_transform.m
│       ├── brightness.m
│       ├── dataframe2struct.m
│       ├── deunderscore.m
│       ├── gotoDir.m
│       ├── indicatorMatrix.m
│       ├── map2vec.m
│       ├── map2vol.m
│       ├── mask2roi.m
│       ├── overwritePrompt.m
│       ├── pairMatrix.m
│       ├── rankTransform_equalsStayEqual.m
│       ├── rankTransform.m
│       ├── relRankIn_includeValue_lowerBound.m
│       ├── removeSpaces.m
│       ├── replaceWildcards.m
│       ├── scale01.m
│       ├── setIfUnset.m
│       ├── struct2dataframe.m
│       ├── underscoresToSpaces.m
│       ├── upSample.m
│       └── vec2map.m
└── tree.md

15 directories, 179 files

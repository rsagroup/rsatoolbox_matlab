%  modelRDMs is a user-editable function which specifies the models which
%  brain-region RDMs should be compared to, and which specifies which kinds of
%  analysis should be performed.
%
%  Models should be stored in the "Models" struct as a single field labeled
%  with the model's name (use underscores in stead of spaces).
%  
%  Cai Wingfield 11-2009

function Models = modelRDMs()

nconditions=72;
% define animate and inanimate index vectors
OwnBodyParts = 1:3; OwnFaces = 4:8; OwnPet = 9; OwnPlaces = 10:15; OwnObjects=16:18;
OtherBodyParts = 19:21; OtherFaces = 22:26; OtherPet = 27; OtherPlaces = 28:33; OtherObjects=34:36;
GeneralBodyParts = 37:44; GeneralFaces = 45:52; GeneralPets = [53 54]; GeneralPlaces=55:66; GeneralObjects=67:72;

animates = [OwnBodyParts OwnFaces OwnPet OtherBodyParts OtherFaces OtherPet GeneralBodyParts GeneralFaces GeneralPets];
inanimates = [OwnPlaces OwnObjects OtherPlaces OtherObjects GeneralPlaces GeneralObjects];%setdiff(1:nconditions,animates)

faces = [OwnFaces OtherFaces GeneralFaces];
nonfaces = setdiff(1:nconditions,faces);

Models.facenonface = ones(nconditions,nconditions);
Models.facenonface(faces,faces)=0;
Models.facenonface(nonfaces,nonfaces)=0;
Models.facenonface(logical(eye(nconditions)))=0; % fix the zero-diagonal


Models.animateInanimate = ones(nconditions,nconditions);
Models.animateInanimate(animates,animates)=0;
Models.animateInanimate(inanimates,inanimates)=0;
Models.animateInanimate(logical(eye(nconditions)))=0; % fix the zero-diagonal

% familiar = [OwnBodyParts OwnFaces OwnPet OwnPlaces OwnObjects];
% unfamiliar = [OtherBodyParts OtherFaces OtherPet OtherPlaces OtherObjects];
% 
% Models.familiarUnfamiliar = nan(nconditions,nconditions);
% Models.familiarUnfamiliar(familiar,familiar)=0;
% Models.familiarUnfamiliar(unfamiliar,unfamiliar)=0;
% Models.familiarUnfamiliar(familiar,unfamiliar)=1;
% Models.familiarUnfamiliar(unfamiliar,familiar)=1;
% Models.familiarUnfamiliar(logical(eye(nconditions)))=0; % fix the zero-diagonal

Models.random = squareform(pdist(rand(nconditions,nconditions)));

% Models.noStructure = ones(64,64);
% Models.noStructure(logical(eye(64)))=0;

end%function

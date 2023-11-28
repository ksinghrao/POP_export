% Takes image plots non axial slices and plots them

% Go to path with dicom image and rtst 
pathIn = '/Users/kamal_singhrao/Desktop/Pt files/Non axial slices/2017-01__Studies';
cd(pathIn);

% Load image and reinterpolate
dicomIMdir = dir([pathIn filesep 'IM/' filesep '*.dcm']);

for iter1 = 1:length(dicomIMdir)
    IMs = dicomread([dicomIMdir(iter1).folder filesep dicomIMdir(iter1).name]);
    IMsInfo = dicominfo([dicomIMdir(iter1).folder filesep dicomIMdir(iter1).name]);
    IM(:,:,IMsInfo.InstanceNumber) = IMs;
    instanceNumber(iter1) = IMsInfo.InstanceNumber;
    SOPClassUID{iter1} = IMsInfo.SOPInstanceUID;
end

sIM = size(IM);

% Overlay contours in axial plan
dicomInfoDir = dir([pathIn filesep 'rtst/' filesep '*.dcm']);

rtst = dicominfo([dicomInfoDir.folder filesep dicomInfoDir.name]);

contours = dicomContours(rtst);

imagePositionPat = IMsInfo.ImagePositionPatient;

%% Match contour instance ID to image set
% Pick a contour, start with 8 because its the PTV
roiSet = rtst.ROIContourSequence;
roiSetCell = struct2cell(roiSet);

contourSet = cell(length(roiSetCell),3);
%contourItem = rtst.ROIContourSequence.Item_8;

for cIter = 1:length(roiSetCell)
    contourItem = roiSetCell{cIter};
    contourItemCell = struct2cell(contourItem);
    
    
    contourN = zeros(sIM);
    
    % Pick a slice of the contour
    contourSequence = struct2cell(contourItem.ContourSequence);
    %contourSlice = contourItem.ContourSequence.Item_1;
    sizeContourSequence = length(contourSequence);
    
    contourCoordinates = cell(sIM(3),1);
    for iter2 = 1:sizeContourSequence
        contourSlice = contourSequence{iter2};
        
        % Find Associated CT slice with contour
        refSOPUID = contourSlice.ContourImageSequence.Item_1.ReferencedSOPInstanceUID;
        
        % Find corresponding image slice
        index = find(contains(SOPClassUID,refSOPUID));
        sliceN = instanceNumber(index);
        
        % Extract contour data
        contourData = contourSlice.ContourData;
        nRows = reshape(contourData,3,round(length(contourData)/3));
        imOffset = nRows - IMsInfo.ImagePositionPatient;
        
        % voxel size
        imVoxSize = [IMsInfo.PixelSpacing(1) IMsInfo.PixelSpacing(1) IMsInfo.SliceThickness]';
        contourPixelCoordinates = imOffset./imVoxSize;
        
        xPos = round(contourPixelCoordinates(1,:));
        yPos = round(contourPixelCoordinates(2,:));
        
        zPos = sliceN(1)*ones(size(yPos));
        % Write xy coordinates
        contourCoordinates{sliceN(1)} = [xPos; yPos; zPos]';

        % Write mask
        contourN(:,:,sliceN(1)) = poly2mask(xPos,yPos,512,512);
    end
    
    % Extract ROI name
    sSetROISequence = struct2cell(rtst.StructureSetROISequence);
    roiName = sSetROISequence{cIter}.ROIName;
    
    contourSet{cIter,3} = contourCoordinates;
    contourSet{cIter,2} = contourN;
    contourSet{cIter,1} = roiName;
end


%% For each structure set extract 3D binary mask
% List strcutures of interst
sPlot = [3,4,8,13,14,15,16];

% Load the binary mask of each contour

for sNum = 1:length(sPlot)
    bVol{sNum} = contourSet{sPlot(sNum),2};
end

%% Visualize images
% Find centroid of PTV
structID = find(contains({contourSet{:,1}},'PTV'));
centerSliceLocSt = regionprops(contourSet{8,2},'centroid');
centerSliceLoc = round(centerSliceLocSt.Centroid);

f = figure;
f.Position = [50 50 1600 400];tiledlayout(1,3)

% Axial
nexttile
imagesc(IM(:,:,centerSliceLoc(3)),[-600 3000])
colormap('gray')
axis off;


%bCor = bVol{1}(:,:,centerSliceLoc(3));
%[B,L] = bwboundaries(bCor,'noholes');

%hold on;
%plot(B{1}(:,2),B{1}(:,1))

for iterAx = 1:length(bVol)
    bAx = bVol{iterAx}(:,:,centerSliceLoc(3));
    [B,L] = bwboundaries(bAx,'noholes');
    
    if isempty(B)
        continue;
    end
    hold on;
    p = plot(B{1}(:,2),B{1}(:,1));
    p.LineWidth = 2;
end


% Sagittal
nexttile
imagesc(rot90(squeeze(IM(:,centerSliceLoc(1),:)),3),[-600 3000])
colormap('gray')
axis off;

for iterSag = 1:length(bVol)
    bCor = rot90(squeeze(bVol{iterSag}(:,centerSliceLoc(1),:)),3);
    [B,L] = bwboundaries(bCor,'noholes');
    
    if isempty(B)
        continue;
    end
    hold on;
    p = plot(B{1}(:,2),B{1}(:,1));
    p.LineWidth = 2;
end
% Coronal
nexttile
imagesc(rot90(squeeze(IM(centerSliceLoc(2),:,:)),3),[-600 3000])
colormap('gray')
axis off;

for iterCor = 1:length(bVol)
    bCor = rot90(squeeze(bVol{iterCor}(centerSliceLoc(2),:,:)),3);
    [B,L] = bwboundaries(bCor,'noholes');
    
    if isempty(B)
        continue;
    end
    hold on;
    p = plot(B{1}(:,2),B{1}(:,1));
    p.LineWidth = 2;
end


%{
nexttile
imagesc(rot90(squeeze(IM(centerSliceLoc(1),:,:)),3))
colormap('gray')
axis off;

%}

%%
%{
%% Visualize images
% List strcutures of interst
sPlot = [3,4,8,13,14,15,16];

% Find centroid of PTV
structID = find(contains({contourSet{:,1}},'PTV'));
centerSliceLocSt = regionprops(contourSet{8,2},'centroid');
centerSliceLoc = round(centerSliceLocSt.Centroid);

figure;

imagesc(IM(:,:,centerSliceLoc(3)))
colormap('gray')
axis off;

for iter3 = 1:length(sPlot)
    hold on;
    xPoints = contourSet{sPlot(iter3),3}{centerSliceLoc(3)};
    yPoints =contourSet{sPlot(iter3),3}{centerSliceLoc(3)} ;

    if isempty(xPoints)
        continue;
    end

    h = plot(contourSet{sPlot(iter3),3}{centerSliceLoc(3)}(:,1), ...
        contourSet{sPlot(iter3),3}{centerSliceLoc(3)}(:,2));

    set(h(1),'linewidth',2);
end

%% For each contour find the ones at centerslice location
imagesc(squeeze(IM(:,centerSliceLoc(2),:)))
colormap('gray')
axis off;

%for iter4 = 1:length(sPlot)
    iter4=3;
    hold on;
    cPoints = contourSet{sPlot(iter4),3};
    
    xPointsR = [];
    yPointsR = [];
    zPointsR = [];

    for iter5 = 1:length(cPoints)
        if isempty(cPoints{iter5})
            continue;
        end
            xPoints = cPoints{iter5}(:,1);
            yPoints = cPoints{iter5}(:,2);
            zPoints = cPoints{iter5}(:,3);

            % Find points in slice
            sLocs = find(xPoints==centerSliceLoc(2));
            yPointsR = [yPoints(sLocs); yPointsR];
            zPointsR = [zPoints(sLocs); zPointsR];


  %  end

    h = plot(xPointsR,yPointsR);
    set(h(1),'linewidth',2);
end

%%
%{
%%

contCoordinates = nan(1,3);

for i1 = 1:length(contourSet{sPlot(iter3),3})
    if isempty(contourSet{sPlot(iter3),3}{i1})
        continue;
    end
    
    contCoordinates = [contourSet{sPlot(iter3),3}{i1}; contCoordinates];
end

contCoordinates(end,:) = [];

[k1,av1] = convhull(contCoordinates(:,1),contCoordinates(:,2),contCoordinates(:,3));

x = contCoordinates(:,1);
y = contCoordinates(:,2);
z = contCoordinates(:,3);

trisurf(k1,x,y,z,'FaceColor','cyan')

indexCV = find(z==70);

%%
%}
%}
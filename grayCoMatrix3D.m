function [ coMat, harMat ] = grayCoMatrix3D( I,distance,numLevels,numHarFeature,...
    offSet )
%Gray level cooccurance 3D
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%I = the 3D image matrix
%distance = a vector of the distances to analyze in
%numLevels = the number of graylevels to be used
%numHarFeature = the number of haralick features to compute
%offSet = a matrix of the directions to analyze in
%
%harMat = a matrix of the haralick features in the format harMat(direction,
%feature)
%distance: a nx1 array of distances that will be used when analyzing the
%image. Default is [1,2,4,8];
%coMat the Co-Occurrence matrices produced
%direction: a nx3 array of direction offsets in [row, column, vertical]
%format. The vertical value increases from top to bottom 
%   [0 1 -1]   0 degrees, 45 degrees
%   [0 0 -1]   straight up
%   [0 -1 -1]  0 degrees, 135 degrees
%   [-1 0 -1]  90 degrees, 45 degrees
%   [1 0 -1]   90 degrees, 135 degrees
%   [-1 1 -1]  45 degrees, 45 degrees
%   [1 -1 -1]  45 degrees, 135 degrees
%   [-1 -1 -1] 135 degrees, 45 degrees
%   [1 1 -1]   135 degrees, 135 degrees
%   Default is all 13 directions.
%**************Variable initialization/Declaration**********************
harMat =0;


noDirections = size(offSet,1); %number of directions, currently 13
coMat = zeros(numLevels,numLevels,noDirections);


%************************graylevel resizing*******************************
numLevels = numLevels-1; %don't touch. Logical adding issue.
minImage = min(min(min(I)));
display (minImage);
display (size(I));
I = I-(minImage);
min(min(min(I)));
maxImage = max(max(max(I)));
tempShift = double(maxImage)/double(numLevels);
I = floor(double(I)/double(tempShift));
I=I+1;
numLevels = numLevels+1; %don't touch. Logical adding issue.
if max(max(max(I))) > numLevels
    disp('Error is graylevel resizing.')
    disp('cooc3d.m');
    return
end
%matlabpool open;

%**************************Beginning analysis*************************
%Order of loops: Direction, slice, graylevel, graylevel locations
parfor direction =1:noDirections %currently 13 (for the 3d image)

    tempMat = zeros(numLevels,numLevels,size(I,3));
    for slicej =1:size(I,3)
         for j=1:numLevels %graylevel
             
             %finds all the instances of that graylevel
            [rowj,colj] = find(I(:,:,slicej)==j);  

            %populating the Cooc matrix.
            for tempCount = 1:size(rowj,1) 
                rowT = rowj(tempCount) + distance*offSet(direction,1);
                colT = colj(tempCount) + distance*offSet(direction,2);
                sliceT = slicej+ distance*offSet(direction,3);
                [I1, I2, I3] = size(I);
                if rowT <= I1 && colT <= I2 && sliceT <= I3
                    if rowT > 0 && colT > 0 && sliceT > 0
                        
                        %Error checking for NANs and Infinite numbers
                        IIntensity = I(rowT,colT,sliceT);
                        if ~isnan(IIntensity)
                            if ~isinf(IIntensity)
                                %Matlab doesn't have a ++ operator.
                                tempMat(j,IIntensity,slicej)= tempMat...
                                    (j,IIntensity,slicej)+1;
                            end
                        end
                    end
                end
            end
            
        end

    end
    for slicej =1:size(I,3)
        coMat(:,:,direction)= coMat(:,:,direction)+tempMat(:,:,slicej);
    end
    %vectorized version
    %coMat(:,:,direction)= coMat(:,:,direction)+tempMat(:,:,1:size(I,3));
end
%matlabpool close;
%extracting the Haralick features from the Co-Occurrence matrices
harMat = harFeatures(coMat,numHarFeature);
return


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%coMat = Co-occurrence matrices 2D stack upon eachother (i j k) k is number
%of directions analyzed. For 3d that's 13. Created in cooc3d.m
%numHarFeature is the number of variables you will be extracting.
%Haralick order
%Energy, Entropy, Correlation, Contrast, Variance, SumMean, Inertia, 
%Cluster Shade, Cluster tendendy, Homogeneity,MaxProbability, 
%Inverse Variance.
%harMat = matrix of the haralick features in the format harMat(direction,
%feature)
function [harMat]= harFeatures(coMat, numHarFeature)

%numHarFeature=12;
%numPosFeature=12; %If you add any more features bump this up.
numLevels = size(coMat,1); %number of graylevels
harMat = zeros(numHarFeature,size(coMat,3));
%%%%%%tempHarMat = zeros(numPosFeature,1);  %continue working here....
%tempCoMat=zeros(size(coMat,1),size(coMat,2));


for iteration = 1:size(coMat,3) %directions

    
%%%%%%%%%%%%%%%%%%%%Preparation

%%%%%%%determining various p values

    pij = sum(sum(coMat(:,:,iteration))); %already normalized
    coMat(:,:,iteration)=coMat(:,:,iteration)./pij;

    tempmux=0;
    tempmuy=0;
    for j=1:numLevels
        for i=1:numLevels
            tempmux =  tempmux+(i*(coMat(j,i,iteration)));
            tempmuy =  tempmuy+(j*(coMat(j,i,iteration)));
        end
    end
    mux=tempmux; %mux
    muy=tempmuy;

    tempx=0;
    tempy=0;
    for j=1:numLevels
        for i=1:numLevels
            tempx = tempx+ (i-mux)^2*coMat(j,i,iteration);
            tempy = tempy+ (j-muy)^2*coMat(j,i,iteration);
        end
    end
    sigx=tempx; %sigx
    sigy=tempy;
    



%Calculations
    tempEnergy =0;
    tempEntropy=0;
    tempCorr=0;
    tempCont=0;
    tempGen=0;
    tempVar=0;
    tempMean=0;
    tempInert=0;
    tempShade=0;
    tempTen=0;
    tempInVar=0;
    for j=1:numLevels
        for i=1:numLevels
            value = coMat(j,i,iteration);
            
            tempEnergy = tempEnergy+ value^2;
            if(value~=0) 
                tempEntropy = tempEntropy + (value * log10(value));
            end
            tempCorr = tempCorr+ ((i-mux)*(j-muy)*(value/(sigy*sigx)));
            n=(abs(i-j))^2;
            tempCont = tempCont+ value*n;
            tempGen = tempGen+ value/(1+abs(1-j));
            tempVar = tempVar + ((i - mux)^2)*value+((j-muy)^2)*value;
            tempMean = tempMean + (i+j)*(value);
            tempInert = tempInert+ (i-j)^2*(value);
            tempShade=tempShade+ ((i+j-mux-muy)^3)*(value);
            tempTen = tempTen+ (((i + j - mux - muy)^4) .* (value));
            if i~=j
                tempInVar=tempInVar+ value/(i-j)^2;
            end
        end
    end
    harMat(1,iteration)=tempEnergy;         %Energy
    harMat(2,iteration) = -tempEntropy;     %Entropy
    harMat(3,iteration)=tempCorr;           %Correlation
    harMat(4,iteration)=tempCont;           %Contrast
    harMat(5,iteration) = tempGen;          %Homogeneity
    harMat(6,iteration) = tempVar/2;        %Variance
    harMat(7,iteration)=tempMean/2;         %Sum Mean
    harMat(8,iteration)=tempInert;          %Inertia
    harMat(9,iteration)=tempShade;          %Cluster Shade
    harMat(10,iteration) = tempTen;         %Cluster Tendency
    harMat(11,iteration) = max(max(coMat(:,:,iteration))); %Max Probability
    harMat(12,iteration) = tempInVar;       %Inverse Variance
    
    clear 'tempEnergy' 'tempEntropy' 'tempCorr' 'tempCont' 'tempGen';
    clear 'tempVar' 'tempMean' 'tempInert' 'tempShade';
    clear 'tempTen' 'tempInVar';

end
%makes it so that rows are cases
harMat = harMat';
return


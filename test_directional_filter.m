% ========================
% DTM Generation (dsm2dtm)
% ========================
% Here you find a naive and thus simple Matlab implementation of the Digital Terrain Model (DTM) generation algorithm
% proposed in our papers [1,2]. From a given input Digital Surface Model (DSM) a DTM is extracted according to
% Algorithm 1 [2]. An optimized version can be found within the commercial software Remote Sensing Software Graz (RSG) [3].
% There the version using two passes (cf. Figure 16, [2]) with only one label image is implemented.
%
% If you use the code or any results in your work then refer to our papers [1,2].
%
% ==========
% References
% ==========
% [1] Roland Perko, Hannes Raggam, Karlheinz Gutjahr, and Mathias Schardt. Advanced DTM generation from very high resolution
%     satellite stereo images. In ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences,
%     volume II-3/W4, pages 165-172, Munich, Germany, March 2015.
%
%     @InProceedings{perko2015advanced,  
%       title                    = {Advanced {DTM} generation from very high resolution satellite stereo images},  
%       author                   = {Roland Perko and Hannes Raggam and Karlheinz Gutjahr and Mathias Schardt},  
%       booktitle                = {ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences},  
%       year                     = {2015},  
%       address                  = {Munich, Germany},  
%       month                    = {March},  
%       pages                    = {165-172},  
%       volume                   = {II-3/W4},  
%       doi                      = {10.5194/isprsannals-II-3-W4-165-2015},  
%       timestamp                = {2015-03-26}}  
% 
% [2] Roland Perko, Hannes Raggam, and Peter M. Roth. Mapping with Pleiades - End-to-end workflow. Remote Sensing,
%     11(17):2052, September 2019.
% 
%     @Article{perko2019mapping,  
%       author    = {Roland Perko and Hannes Raggam and Peter M. Roth},  
%       title     = {Mapping with {P}l{\'e}iades -- {E}nd-to-End Workflow},  
%       journal   = {Remote Sensing},  
%       year      = {2019},  
%       volume    = {11},  
%       number    = {17},  
%       pages     = {2052},  
%       issn      = {2072-4292},  
%       doi       = {10.3390/rs11172052},  
%       url       = {https://www.mdpi.com/2072-4292/11/17/2052},   
%       month     = {September},  
%       timestamp = {2019-09-01}}  
% 
% [3] https://www.remotesensing.at/remote-sensing-software/
% 
%
% Roland Perko, 2019
%
function [] = test_directional_filter(sName,iVerbose)

%close all; clear;

disp([' --- DSM2DTM conversion ' sName ' ---']);

bNoCorrection  = false; % also extract DTM without slope correction
bOnlyLeftRight = false; % take only left and right direction
bSmooth        = true;

dSpacing       =    1.0; % one meter spacing (i.e. GSD of input DSM)

% parameter setup
iDistance      =     91; % distance to search for minimum
dSigma         =     25; % sigma for gaussian bluring
dThrDeltaMin   =    3.0; % height difference to minimal value in meters (height_treshold in RSG)
                         % a large value yield more ground pixels
dThrDeltaSlope =     30; % max slope in degrees
iMinGround     =      3; % in RSG this number is always larger by 1. e.g. Matlab 3 is 4 in RSG
dOutputNodata  = -99999; % the output nodata value

if ( bOnlyLeftRight )
    iMinGround = 0;
end

%             X  Y
iDir(1,:) = [ 1  0];
iDir(2,:) = [ 1  1];
iDir(3,:) = [ 0  1];
iDir(4,:) = [-1  1];
iDir(5,:) = [-1  0];
iDir(6,:) = [-1 -1];
iDir(7,:) = [ 0 -1];
iDir(8,:) = [ 1 -1];

%     7
%   6 | 8
% 5 - o - 1
%   4 | 2
%     3

oDSM = imread(sName);
%oDSM = oDSM(1:1024,1:1024);
oDSM = double(oDSM);

if ( bSmooth )
    tic;
    disp('DSM smoothing');
    
    iSize = max(1,fix(6*dSigma));
    iSize = iSize + (1 - mod(iSize,2)); % must be odd
    iSizeH = (iSize-1)/2;
    
    h = fspecial('gaussian',[iSize iSize],dSigma);
    oDSMs = imfilter(oDSM,h,'symmetric');
    fprintf(' %6.2f secs\n',toc);
    if ( iVerbose > 0 )
        sNameGauss = [sName(1:end-4) '_test_gauss.tif'];
        imwrite2tif(oDSMs,[],sNameGauss,'single');
    end
end


iDim = size(oDSM);

% get slope in all 8 directions
sNameDelta       = [sName(1:end-4) '_test_delta.tif'];
sNameDeltaSmooth = [sName(1:end-4) '_test_delta_smooth.tif'];

%if ( true )
if ( exist(sNameDelta,'file') == 0 )
    
    disp('DSM slope extraction');
    
    % get slope in all 8 directions
    oDD  = zeros(iDim(1),iDim(2),8,'single');
    oDDs = zeros(iDim(1),iDim(2),8,'single');
    iBorder = 10;
    tic;
    for iY = 1+iBorder:iDim(1)-iBorder
        if ( mod(iY,ceil(iDim(1)/64)) == 0 )
            fprintf('*');
        end
        for iX = 1+iBorder:iDim(2)-iBorder
            for iD = 1:8
                oDD(iY,iX,iD)  = oDSM(iY,iX)  - oDSM(iY+iDir(iD,2),iX+iDir(iD,1));
                oDDs(iY,iX,iD) = oDSMs(iY,iX) - oDSMs(iY+iDir(iD,2),iX+iDir(iD,1));
            end
        end
    end
    fprintf(' %6.2f secs\n',toc);
    
    
    %     oDMin  =  min(abs(oDelta),[],3);
    %     oDMax  =  max(abs(oDelta),[],3);
    %     oDMean = mean(abs(oDelta),3);
    %     oDStd  =  std(abs(oDelta),[],3);
    %
    %     figure; imshow(oDMean,[]); title('mean');
    %     figure; imshow(oDStd,[]);  title('std');
    %     figure; imshow(oDMin,[]);  title('min');
    %     figure; imshow(oDMax,[]);  title('max');
    
    %    imwrite2tif(oDMin ,[],[sName(1:end-4) '_test_min.tif'],'single');
    %    imwrite2tif(oDMax ,[],[sName(1:end-4) '_test_max.tif'],'single');
    %    imwrite2tif(oDMean,[],[sName(1:end-4) '_test_mean.tif'],'single');
    %    imwrite2tif(oDStd ,[],[sName(1:end-4) '_test_std.tif'],'single');
    if ( iVerbose > 0 )
        disp('  store images');
        tic;
        imwrite2tif(oDD,[],sNameDelta,'single');
        imwrite2tif(oDDs,[],sNameDeltaSmooth,'single');
        fprintf('\n');
        fprintf(' %6.2f secs\n',toc);
    end
    
else
    
    %    oDMin  = imread([sName(1:end-4) '_test_min.tif']);
    %    oDMax  = imread([sName(1:end-4) '_test_max.tif']);
    %    oDMean = imread([sName(1:end-4) '_test_mean.tif']);
    %    oDStd  = imread([sName(1:end-4) '_test_std.tif']);
    oDD    = imread(sNameDelta);
    oDDs   = imread(sNameDeltaSmooth);
    
end

% get minimum filter

%fun1 = @(x) min(x(11,:));
%fun2 = @(x) min(x(:,11));
%fun3 = @(x) min(diag(x));
%fun4 = @(x) min(diag(fliplr(x)));
%oMin1=nlfilter(oDSM,[21 21],fun1);
%oMin2=nlfilter(oDSM,[21 21],fun2);
%oMin3=nlfilter(oDSM,[21 21],fun3);
%oMin4=nlfilter(oDSM,[21 21],fun4);


% get minimal value in 4 directions
sNameMin     = [sName(1:end-4) '_test_minimum.tif'];
sNameMinCorr = [sName(1:end-4) '_test_minimum_corrected.tif'];

%if ( exist(sNameMin,'file') == 0 )
if ( true )
    disp('DSM slope dependent minima extraction');
    
    nhood(1) = iDistance;
    nhood(2) = nhood(1);
    nhoodH   = (nhood(1)-1)/2;
    padval   = 0;
    
    % Expand A
    [ma,na] = size(oDSM);
    % aa <- DSM
    aa = zeros(size(oDSM)+nhood-1,'single');
    aa(floor((nhood(1)-1)/2)+(1:ma),floor((nhood(2)-1)/2)+(1:na)) = oDSM;
    
    % bb <- DSM smoothed
    bb = zeros(size(oDSM)+nhood-1,'single');
    bb(floor((nhood(1)-1)/2)+(1:ma),floor((nhood(2)-1)/2)+(1:na)) = oDSMs;
    
    % Find out what output type to make.
    rows = 0:(nhood(1)-1);
    cols = 0:(nhood(2)-1);
    
    if ( bNoCorrection )
        oMin = zeros([size(oDSM) 8],'single');
    end
    oMinCorr = zeros([size(oDSM) 8],'single');
    
    % calculate slope corrected minimum
    X = -nhoodH:nhoodH;
    
    
    %     7
    %   6 | 8
    % 5 - o - 1
    %   4 | 2
    %     3
    
    % Apply fun to each neighborhood of a
    tic;
    for i=1:ma
        if ( mod(i,ceil(iDim(1)/64)) == 0 )
            fprintf('*');
        end
        for j=1:na
            x  = aa(i+rows,j+cols); % DSM
            %xs = bb(i+rows,j+cols); % DSMs
            
            % simple minimum (here two directions share the same minimum)
            x1 = x(nhoodH+1,:);
            x2 = diag(x)';
            x3 = x(:,nhoodH+1)';
            x4 = diag(fliplr(x))';
            
%            if ( i == 879 && j == 3439 ) % x==j, y==i
%                                         % just for illustration of the
%                                         % robust fit
%
%                %xs = bb(i+rows,j+cols); % DSMs
%                xs = -oDDs(i+rows,j+cols,3);
%                figure;
%                subplot(1,2,1); imshow(x ,[]); hold on;
%                line([nhoodH nhoodH],[1 nhood(1)],'Color','g');
%                plot(nhoodH(1),nhoodH(1),'y+');
%                plot(nhoodH(1),nhood(1)-4,'gv');
%                
%                subplot(1,2,2); imshow(xs,[]); hold on;
%                line([nhoodH nhoodH],[1 nhood(1)],'Color','g');
%                plot(nhoodH(1),nhoodH(1),'y+');
%                plot(nhoodH(1),nhood(1)-4,'gv');
%                
%                set(findobj('type','line'),'linewidth',3);
%                set(gcf,'Color','white');
%                savefigure(['D:\per\documents\Projekt-2013-HighSensII\DSM2DTM\DSM2DTM_smoothing.png'],true,0);
%                         
%                
%                % visulization for direction 3 (=down)
%                Yc = x3 + X.*oDDs(i,j,3); % slope corrected DSM
%                
%                [b,stats] = robustfit(X,x3);
%                YrLine = b(1)+b(2)*X; % robust line fit
%                
%                
%                YsLine = b(1) - X.*oDDs(i,j,3); % gaussian based line fit
%                
%                iFigurePosition = [100 250 800 500];
%                figure('Position',iFigurePosition);
%                grid on; hold on;
%                plot(X,x3,'k');
%                plot(X,YrLine,'g--')
%                plot(X,YsLine,'c:');                
%                plot(X,x3-YrLine+YrLine(nhoodH+1),'r--');
%                plot(X,Yc,'b:');
%                                
%                legend('input data','robust fit','simple fit','robustly corrected data','simply corrected data','Location','SouthWest');
%                set(findobj('type','line'),'linewidth',3);
%                set(gcf,'Color','white');
%                savefigure(['D:\per\documents\Projekt-2013-HighSensII\DSM2DTM\DSM2DTM_1D.png'],true,0);
%               
%            end

if ( bNoCorrection )
    oMin(i,j,1) = min(x1);
    oMin(i,j,5) = oMin(i,j,1);
    oMin(i,j,3) = min(x3);
    oMin(i,j,7) = oMin(i,j,3);
    oMin(i,j,2) = min(x2);
    oMin(i,j,6) = oMin(i,j,2);
    oMin(i,j,4) = min(x4);
    oMin(i,j,8) = oMin(i,j,4);
end
            
%             if ( iVerbose > 0 )
%                 %if ( i==325 && j==777)
%                 %if ( i==2097 && j==1142 )
%                 if ( i==1150 && j==360 )
%                     figure;
%                     subplot(1,2,1); imshow(x ,[]);
%                     subplot(1,2,2); imshow(xs,[]);
%                     impixelinfo;
%                     % robust linear fit  and minimum
%                     X = -nhoodH:nhoodH;
%                     Y = x(nhoodH,:);
%                     [b,stats] = robustfit(X,Y);
%                     
%                     %dGradient = oDDs(i+nhoodH,j+nhoodH,1);
%                     dGradient = oDDs(i,j,1);
%                     Yd = X.*dGradient;
%                     Ycorr1 = Y + Yd;
%                     
%                     dGradient = -oDDs(i,j,5);
%                     Yd = X.*dGradient;
%                     Ycorr5 = Y + Yd;
%                     
%                     Yr = b(1)+b(2)*X;
%                     figure;
%                     grid on; hold on;
%                     plot(X,Y);
%                     plot(X,Yr,'g','LineWidth',2)
%                     plot(X,Y-Yr+Yr(nhoodH+1),'r--');
%                     plot(X,Ycorr1,'m:');
%                     plot(X,Ycorr5,'k:');
%                     legend('input data','robust fit','correct (fit)','corrected (slope)','corrected (slope)');
%                     
%                     1
%                     
%                 end
%             end
            
            % direction 1 & 5
            oMinCorr(i,j,1) = min(x1 + X.*oDDs(i,j,1));
            oMinCorr(i,j,5) = min(x1 + X.*(-oDDs(i,j,5)));
            
            % direction 2 & 6
            oMinCorr(i,j,2) = min(x2 + X.*oDDs(i,j,2));
            oMinCorr(i,j,6) = min(x2 + X.*(-oDDs(i,j,6)));
            
            % direction 3 & 7
            oMinCorr(i,j,3) = min(x3 + X.*oDDs(i,j,3));
            oMinCorr(i,j,7) = min(x3 + X.*(-oDDs(i,j,7)));
            
            % direction 4 & 8
            oMinCorr(i,j,4) = min(x4 + X.*oDDs(i,j,4));
            oMinCorr(i,j,8) = min(x4 + X.*(-oDDs(i,j,8)));
        end
    end
    fprintf(' %6.2f secs\n',toc);
    
    clear aa bb;
    %clear oDDs;
    
    if ( iVerbose > 0 )
        if ( bNoCorrection )
            imwrite2tif(oMin ,[],sNameMin,'single');
        end
        imwrite2tif(oMinCorr ,[],sNameMinCorr,'single');
    end
    
else
    
    %oMin     = imread(sNameMin);
    oMinCorr = imread(sNameMinCorr);
    
end


% get simple labeling in x-direction

disp('DSM classification');

iGround    = 1;
iNonGround = 2;

%iD       = 1; % direction 1 --> right
%iIDDelta = 1; % direction 1 --> right

iBorder = 10;
tic;
% loop over two types of minima extraction

if ( bNoCorrection )
    iTypeStart = 1;
else
    iTypeStart = 2;
end

for iType = iTypeStart:2 %1:2
    if ( iType == 1 )
        fprintf('[1] ');
    else
        fprintf('[2] ');
    end
    
    oLabel = zeros([iDim(1) iDim(2) 8],'uint8');
    %oSlope = zeros(iDim(1),iDim(2),8,'single');
    
    % loop over first 4 directions (PASS 1)
    fprintf('(P1)');
    for iD = 1:4
        
        for iY = 1+iBorder:iDim(1)-iBorder
            if ( mod(iY,ceil(iDim(1)/4)) == 0 )
                fprintf('*');
            end
            for iX = 1+iBorder:iDim(2)-iBorder
                
                % height difference to lowest point
                if ( iType == 1 )
                    dHDelta = ( oDSM(iY,iX) - oMin(iY,iX,iD) );
                else % use slope corrected minimum
                    dHDelta = ( oDSM(iY,iX) - oMinCorr(iY,iX,iD) );
                end
                
                if ( dHDelta > dThrDeltaMin )
                    oLabel(iY,iX,iD) = iNonGround;
                else
                    
                 if ( iType == 1 )
                    dDelta = oDD(iY,iX,iD); % + means down, - means up
                    dSignDelta = -sign(dDelta);
                    dSlopeLocal = atan2(abs(dDelta),dSpacing)*180/pi;
                    dSlope = dSlopeLocal * dSignDelta;
                else % use slope corrected local slope
                    dDelta = oDD(iY,iX,iD)-oDDs(iY,iX,iD); % + means down, - means up
                    dSignDelta = -sign(dDelta);
                    dSlopeLocal = atan2(abs(dDelta),dSpacing)*180/pi;
                    dSlope = dSlopeLocal * dSignDelta;                    
                end                   
                    
                    if ( dSlope > dThrDeltaSlope )
                        oLabel(iY,iX,iD) = iNonGround;
                    else
                        oLabel(iY,iX,iD) = oLabel(iY-iDir(iD,2),iX-iDir(iD,1),iD); % label as previous point
                    end
                    
                    if ( dSlope < 0 ) % negative slope
                        % get nearest ground point
                        oLabel(iY,iX,iD) = iGround;
                    end
                end
                
            end
        end
        
    end
    
    
    % loop over second 4 directions (PASS 2)
    fprintf('(P2)');
    for iD = 5:8
        
        for iY = iDim(1)-iBorder:-1:1+iBorder % 1+iBorder:iDim(1)-iBorder
            if ( mod(iY,ceil(iDim(1)/4)) == 0 )
                fprintf('*');
            end
            for iX = iDim(2)-iBorder:-1:1+iBorder % 1+iBorder:iDim(2)-iBorder
                
                % height difference to lowest point
                if ( iType == 1 )
                    dHDelta = ( oDSM(iY,iX) - oMin(iY,iX,iD) );
                else
                    %        dHDelta = abs( oDSM(iY,iX) - oMinCorr(iY,iX,iD) );
                    dHDelta = ( oDSM(iY,iX) - oMinCorr(iY,iX,iD) );
                    
                end
                
                if ( dHDelta > dThrDeltaMin )
                    oLabel(iY,iX,iD) = iNonGround;
                else
                    
                 if ( iType == 1 )
                    dDelta = oDD(iY,iX,iD); % + means down, - means up
                    dSignDelta = -sign(dDelta);
                    dSlopeLocal = atan2(abs(dDelta),dSpacing)*180/pi;
                    dSlope = dSlopeLocal * dSignDelta;
                else % use slope corrected local slope
                    dDelta = oDD(iY,iX,iD)-oDDs(iY,iX,iD); % + means down, - means up
                    dSignDelta = -sign(dDelta);
                    dSlopeLocal = atan2(abs(dDelta),dSpacing)*180/pi;
                    dSlope = dSlopeLocal * dSignDelta;
                end                   
                    if ( dSlope > dThrDeltaSlope )
                        oLabel(iY,iX,iD) = iNonGround;
                    else
                        oLabel(iY,iX,iD) = oLabel(iY-iDir(iD,2),iX-iDir(iD,1),iD); % label as previous point
                    end
                    
                    if ( dSlope < 0 ) % negative slope
                        % get nearest ground point
                        oLabel(iY,iX,iD) = iGround;
                    end
                end
                
            end
        end
        
    end
    
    %fprintf('\n');
    
    oL0 = (oLabel == 0);
    oL1 = (oLabel == 1);
    oL2 = (oLabel == 2);
    
    if ( bOnlyLeftRight )
        oL0d = oL0(:,:,1) + oL0(:,:,5); % only left and right
        oL1d = oL1(:,:,1) + oL1(:,:,5); % only left and right
        oL2d = oL2(:,:,1) + oL2(:,:,5); % only left and right
    else
        oL0d = sum(oL0,3);
        oL1d = sum(oL1,3);
        oL2d = sum(oL2,3);
    end
    
    oG = uint8(oL1d > iMinGround );
    
    fprintf(' %6.2f secs\n',toc);

    if ( iType == 1 )
        sType = '';
    else
        sType = '_corr';
    end
    
    oGE = imerode(oG,strel('disk',1));
    
    oDTM = oDSM;
    oDTM(oG==0) = dOutputNodata;
    
    if ( iVerbose > 0 )
%        imwrite2tif(oLabel,[],[sName(1:end-4) '_test_label' sType '.tif'],'single');
        imwrite(oG,[sName(1:end-4) '_test_ground' sType '.tif']);
        imwrite(oGE,[sName(1:end-4) '_test_ground_erode' sType '.tif']);
%        imwrite2tif(oSlope,[],[sName(1:end-4) '_test_slope' sType '.tif'],'single');
%        imwrite(uint8(oL1d),[sName(1:end-4) '_test_label_sum' sType '.tif']);
    end

    imwrite2tif(oDTM,[],[sName(1:end-4) '_test_DTM' sType '.tif'],'single');
    
    % interpolate the DTM
    tic;
   
    disp('interpolate resulting DTM (fill holes)');    
    [X,Y] = meshgrid(1:iDim(2),1:iDim(1));
    X = X(oG~=0);
    Y = Y(oG~=0);
    oDTM = oDTM(oG~=0);
    [Xq,Yq] = meshgrid(1:iDim(2),1:iDim(1));
    oDTMfill = griddata(X,Y,oDTM,Xq,Yq,'cubic');
    oDTMfill(isnan(oDTMfill))=dOutputNodata; % set nan to nodata value
    fprintf(' %6.2f secs\n',toc);

    imwrite2tif(oDTMfill,[],[sName(1:end-4) '_test_DTM_fill' sType '.tif'],'single');

end
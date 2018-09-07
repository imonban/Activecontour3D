function [phi,MEAN, STD, T] = active_contour3Dnew(im,phi,mu,nu, lambda1,lambda2, dt, NarrowBandSize,print,normalize_force_flag)

% The code is chan and vese 3D segmentation with the improvements of narrow-band and normalization

%  big mu - less curve
% shift operations
%shiftD = @(M) M([1 1:end-1],:,:);
%shiftL = @(M) M(:,[2:end end],:);
%shiftR = @(M) M(:,[1 1:end-1],:);
%shiftU = @(M) M([2:end end],:,:);
% derivatives
%Dx = @(M) (shiftL(M) - shiftR(M))/2;
%Dy = @(M) (shiftU(M) - shiftD(M))/2;

ii=1;                       % init iterations
[dimy, dimx, dimz] = size(im);    % 3D imsize
ObjectSizeX=dimx;
ObjectSizeY=dimy;
ObjectSizeZ=dimy;
% compute Global Contrast and Correlation of 3D volume
PxDist = 1;
offSet = [0 PxDist 0; -PxDist PxDist 0; -PxDist 0 0; -PxDist -PxDist 0]; %2D Co-Occurrence directions
%o,45,90,135 degrees

%the additional 9 directions that make 3D Co-Occurrence from 2D
dimension3 = [0 PxDist -PxDist; 0 0 -PxDist; 0 -PxDist -PxDist; -PxDist 0 -PxDist; PxDist 0 -PxDist; -PxDist PxDist -PxDist; PxDist -PxDist -PxDist;
           -PxDist -PxDist -PxDist; PxDist PxDist -PxDist];
offSet = cat(1,offSet,dimension3);

%[coocMat] = cooc3dTest (im, offSet);
[GTSDM] = GrayCoM3D(im, offSet, 8);
stats = graycoprops(GTSDM);

WholeCont=mean(stats.Contrast); % vector -> one value per angle
WholeHomo=mean(stats.Homogeneity); % vector -> one value per angle % homogenity
clear GTSDM stats dimension3 offSet;

% parfor init
force = 0;
if strcmp(print,'on')
   figure
end

%while (ii<=numIter)
while (1)
    
    % init narrow band size
    idx = find(phi <= NarrowBandSize & phi >= -NarrowBandSize);  % get the curve's narrow band indexes
    %idx = find(phi == NarrowBandSize );
    while (isempty(idx))
        NarrowBandSize=NarrowBandSize+0.05;
        idx = find(phi <= NarrowBandSize & phi >= -NarrowBandSize);
    end
    NarrowBandSize = 1;
    
    [y, x, z] = ind2sub(size(phi),idx); % points in narrow band
  
    xOrig=x;
    yOrig=y;
    zOrig=z;
    idxOrig=idx;
    
    for m=1:length(x)
        DelId = find(x>=max((x(m)-1),1) & x<=min((x(m)+1),size(im,1))& y>=max((y(m)-1),1) & y<=min((y(m)+1),size(im,2))& z>=max((z(m)-1),1) & z<=min((z(m)+1),size(im,3)));
        x(DelId(DelId~=m))=0;
        y(DelId(DelId~=m))=0;
        z(DelId(DelId~=m))=0;
    end
    idnew=find(x~=0);
    xnew=x(idnew);
    ynew=y(idnew);
    znew=y(idnew);
    idxnew=idx(idnew);

    clear x y z idx
    
    x=xnew;
    y=ynew;
    z=znew;
    idx=idxnew;
    
    % myRads is an array for window size per pixel in contour
    % X window size - myRads(1,:)
    % Y window size - myRads(2,:)
    % Z window size - myRads(3,:)
    if (ii == 1)
        myRads(1, 1:length(idx)) = 5;   % local window X
        myRads(2, 1:length(idx)) = 5;   % local window Y
        myRads(3, 1:length(idx)) = 5;   % local window Z
        % for each point in the narrowband, create a cube -> check stats
        % in the square -> decice the force -> update
        if strcmp(print,'on')
            %pause(0.2);
            title(['iteration#',num2str(ii)])
            showCurveAndPhi(im, phi, 'w');
           
            %                 saveas(gcf,['./img/',num2str(ii),'.png'],'png')
        end
 
    else
        % re-init ryRads
        myRads = [];
        TotCont=[];
        TotHomo=[];
        % init i_array
        i_array = 1:length(x);
        %visualize the iso-surface
 %--------------------------------------- 
 % may be it is better to visualize the evolution of the iso-surface
       
        if strcmp(print,'on')
            %pause(0.2);
            title(['iteration#',num2str(ii)])
            showCurveAndPhi(im, phi, 'w');
           
            %                 saveas(gcf,['./img/',num2str(ii),'.png'],'png')
        end
 
%-------------------------------------------------------------------        
        % parfor init
        i_array_indexes = {};
        if matlabpool('size')
            %                 display('use_parfor: YES');
            i=1;
            c = 1;
            while i<=length(i_array)
                step = round(length(i_array)/matlabpool('size'));
                
                i_step = i+step;
                if i_step>length(i_array)
                    i_step = length(i_array);
                end
                
                i_array_indexes{c} = i_array(i:i_step);
                c = c+1;
                i = i+step+1;
            end
        else
            %                 display('use_parfor: NO');
            i_array_indexes{1} = i_array;
        end
        
        parfor i_index=1:length(i_array_indexes)
            for j_index = 1:length(i_array_indexes{i_index});
                i = i_array_indexes{i_index}(j_index);
                mone=1;
                flag = 0;
                myHalfX = 3;
                myHalfY = 3;
                myHalfZ = 3;
                
                while ~flag
                    % mean value of the force controlls the size of the local window (myRads) on each iteration
                    % -> control the convergance of the local window size
                    PrevForce=force;
                    convergence=force;
                    NormConvergence = (convergence + abs(min(convergence)));
                    NormConvergence = NormConvergence / max(NormConvergence(:));
                    meanConv = mean(NormConvergence);
                    
                    % compute the local stats (Contrast) for each pixel in the narrowband
                    % by tacking in account each pixel in the local area (in
                    % and outside the narrow band)
                    SmallArea=(im(round(max(y(i)-myHalfY,1)):round(min(y(i)+myHalfY,dimy)),...
                        round(max(x(i)-myHalfX,1)):round(min(x(i)+myHalfX,dimx)),...
                        round(max(z(i)-myHalfZ,1)):round(min(z(i)+myHalfZ,dimz)))); % local area
     
                    offSet = [0 PxDist 0; -PxDist 0 0; 0 0 -PxDist]; % 0, 90, and stright up
                    %[coocMat] = cooc3dTest (SmallArea, offSet);
                    [GTSDM] = GrayCoM3D(SmallArea, offSet, 8);
                    stats = graycoprops(GTSDM);
                    
                    myRads_x = ObjectSizeX/(log10(ObjectSizeX)*((WholeHomo+stats.Contrast(1)+1/WholeCont)+(sqrt(1/meanConv))));
                    myRads_y = ObjectSizeY/(log10(ObjectSizeY)*((WholeHomo+stats.Contrast(2)+1/WholeCont)+(sqrt(1/meanConv))));
                    myRads_z = ObjectSizeZ/(log10(ObjectSizeZ)*((WholeHomo+stats.Contrast(3)+1/WholeCont)+(sqrt(1/meanConv))));
                    
                    myRads_x=max(2,round(myRads_x));
                    myRads_y=max(2,round(myRads_y));
                    myRads_z=max(2,round(myRads_z));
                    % if windows size changed -> recalculate the window size
                    % until convergance
                    if (abs(myRads_x-myHalfX) == 0 && abs(myRads_y-myHalfY) == 0 && abs(myRads_z-myHalfZ)== 0)
                        flag = 1;
                    else
                        myHalfX = myRads_x;
                        myHalfY = myRads_y;
                        myHalfZ = myRads_z;
                    end
                    
                    % number of iterations for the rectangle to converge
                    mone=mone+1;
                    if (mone==5)
                        flag=1;
                    end
                end
                myRads = [myRads [myRads_x; myRads_y; myRads_z]];
                TotCont=[TotCont WholeCont]; 
                TotHomo=[TotHomo WholeHomo];
            end
        end
    end
    
    display(['ii:',num2str(ii)]);
    
    % constructing the window
    xneg = x - myRads(1,:)';  xpos = x+myRads(1,:)';  %get subscripts for local regions
    yneg = y-myRads(2,:)';  ypos = y+myRads(2,:)';
    zneg = z-myRads(3,:)';  zpos = z+myRads(3,:)';
    
    
    
    xneg(xneg<1)=1; yneg(yneg<1)=1;  zneg(zneg<1)=1;               %check bounds
    xpos(xpos>dimx)=round(dimx); ypos(ypos>dimy)=round(dimy); zpos(zpos>dimz)=round(dimz);
    
    %-- re-initialize u,v,Ain,Aout
    u=zeros(size(idx)); v=zeros(size(idx));
    Ain=zeros(size(idx)); Aout=zeros(size(idx));
    
    %-- compute local stats
    for i = 1:numel(idx)  % for every point in the narrow band
        img = im(yneg(i):ypos(i),xneg(i):xpos(i), zneg(i):zpos(i)); %sub image
        P = phi(yneg(i):ypos(i),xneg(i):xpos(i), zneg(i):zpos(i)); %sub phi
        upts = find(P<=0);            %local interior
        Ain(i) = length(upts);
        if (Ain(i)==0)
            Ain(i)=1;
        end
        u(i) = sum(img(upts))/Ain(i);
        
        vpts = find(P>0);             %local exterior
        Aout(i) = length(vpts);
        if (Aout(i)==0)
            Aout(i)=1;
        end
        v(i) = sum(img(vpts))/Aout(i);
    end
    
    [curv,phi_x2,phi_y2, phi_z2] = get_curvature3D(phi,idx); %get the curvature of the 3D point
    %     curv=(fx2.*fyy + fy2.*fxx -2.*fx.*fy.*fxy)./den;
    grad=phi_x2+phi_y2+phi_z2;
    grad_m=grad.^0.5;
    kappa=curv.*grad_m;
    
    %% CV
    force=  -nu + lambda1*(im(idx)-u).^2 - lambda2*(im(idx)-v).^2; % change the sign of the lambdas to change the evol. direction
    %% Yezzi
    %         force = -nu+lambda1*((im(idx)-u)).^2./Ain-lambda2*((im(idx)-v)).^2./Aout;
    
    % normalize the force
    
    
    if normalize_force_flag
        force = force./max(force(:));
        dphidt = force + mu* kappa;
        dt = .45/(max(dphidt(:))+eps);
    else
        dphidt=force+ mu* kappa;
    end
    ii=ii+1;
    PrevPhi = phi;
    phi(idx)=phi(idx)+dt*dphidt;
    phi = sussman(phi, .5);
    Data.Energy(ii)=sum(force);
    Data.AbsEnergy(ii)=sum(abs(force));
         m = zeros(size(im,1));
         mm = zeros(size(im,1));
         o=zeros(size(im,1));
    %% Comptute the stopping creteria: mean and std inside the lesion
    InT = im(find(phi<=0)); %inside lesion
    MEAN(ii) = mean2(InT);
    STD(ii) = std(InT);
    T(ii) = MEAN(ii) - 2 * STD(ii);
    disp(T(ii));
 if(ii > 100 && abs(T(ii) - T(ii -10)) <= 0.03) % minimum 50 iterations is needed
    break;
end
    %% Input - Bounding box
 %    CI=20;
 %    MinX=round(min(data(:,1)))-CI;
  %   MinY=round(min(data(:,2)))-CI;
   %  MaxX=round(max(data(:,1)))+CI;
    % MaxY=round(max(data(:,2)))+CI;
   %       mm(MinY:MaxY,MinX:MaxX)=phi;
   %       m(mm<0)=1; % Automatic Vs. Manaul
    %      o = roipoly(WholeImg,(data(:,1)),(data(:,2)));  %for Original manual marking
   %       AreaM(i)=sum(sum(m)); 
    %      AreaA(i)=sum(sum(o));
     %    [DiceMine]=sevaluate(m,double(o));
    %     Data.Dice(ii)=DiceMine*100;
end
%          close all
% %         imshow(im,[])
% %         hold on
% %         plot(x,y,'.g')
% %         plot(x(rad<round(alpha)),y(rad<round(alpha)),'.m')
% % %         for i=1:length(x)
% % %             DrawCircle(x(i), y(i), rad(i), 20, 'r')
% % %         end
% end

% if ii>2
%     Rad.Std=mean(StdRad);
%     Rad.Mean=mean(MeanRad);
%     Rad.Max=mean(MaxRad);
%     Rad.Min=mean(MinRad);
%     Rad.Two=sum(RadTwo);
%     Rad.Distribution=MeanRad;
%     % Rad.Homogeneity=mean(RadHomo);
%     Rad.ROISize=mean([dimx dimy]);
%     Rad.Box=[x, y, myRads'];
%     %     Rad.Cont=RadCont;
%     %     Rad.Homo=RadHomo;
%     Data.GlobCont=GlobCont;
     % Data.LocalCont=mean(TotCont);
     % Data.GlobHomo=mean(TotHomo);

%
%     %     Rad.lambda1=mean(SegmentationParams.lambda1);
%     %     Rad.lambda2=mean(SegmentationParams.lambda2);
%     %     Rad.mu=mean(SegmentationParams.mu);
%     %     Rad.nu=mean(SegmentationParams.nu);
%     % Rad.Whole=mean(SegmentationParams.Whole);
% else
%     Rad=[];
%     Data=[];
% end
if strcmp(print,'on')
    pause(0.2);
    showCurveAndPhi(im, phi, 'w');
   title(['iteration#',num2str(ii)])
end
return


%% internal functions
%% internal functions
function D = sussman(D, dt)
% forward/backward differences
a = D - shiftR(D); % backward
b = shiftL(D) - D; % forward
c = D - shiftD(D); % backward
d = shiftU(D) - D; % forward
e = D - shiftF(D); % backward
f = shiftB(D) - D; % forward

a_p = a;  a_n = a; % a+ and a-
b_p = b;  b_n = b;
c_p = c;  c_n = c;
d_p = d;  d_n = d;
e_p = e;  e_n = e;
f_p = f;  f_n = f;  

a_p(a < 0) = 0;
a_n(a > 0) = 0;
b_p(b < 0) = 0;
b_n(b > 0) = 0;
c_p(c < 0) = 0;
c_n(c > 0) = 0;
d_p(d < 0) = 0;
d_n(d > 0) = 0;
e_p(e < 0) = 0; 
e_n(e > 0) = 0; 
f_p(f < 0) = 0;
f_n(f > 0) = 0;

dD = zeros(size(D));
D_neg_ind = find(D < 0);
D_pos_ind = find(D > 0);
dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
    + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2) ...
    + max(e_p(D_pos_ind).^2, f_n(D_pos_ind).^2))- 1;
dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
    + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2) ...
    + max( e_n(D_neg_ind).^2, f_p(D_neg_ind).^2))- 1;

D = D - dt .* sussman_sign(D) .* dD;

%-- whole matrix derivatives
function shift = shiftD(M)
  shift = [ M(1,:,:); M(1:size(M,1)-1,:,:) ];

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2),:) M(:,size(M,2),:) ];


function shift = shiftU(M)
  shift = [ M(2:size(M,1),:,:); M(size(M,1),:,:) ];


function shift = shiftR(M)
  shift = [ M(:,1,:) M(:,1:size(M,2)-1,:) ];
  

function shift = shiftF(M)
 % assumes 3D
  shift = cat(3, M(:,:,1), M(:,:,1:size(M,3)-1));
  
function shift = shiftB(M)
  % assumes 3D
  shift = cat(3, M(:,:,2:size(M,3)), M(:,:,size(M,3)));


function S = sussman_sign(D)
S = D ./ sqrt(D.^2 + 1);

% Calculate curvature
function [curvature,phi_x2,phi_y2, phi_z2] = get_curvature3D(phi,idx)
[dimy, dimx, dimz] = size(phi);
[y,x,z] = ind2sub([dimy,dimx,dimz],idx);  % get subscripts

%-- get subscripts of neighbors
ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1; zm1 = z-1; zp1 = z+1;

%-- bounds checking
ym1(ym1<1) = 1; xm1(xm1<1) = 1; zm1(zm1<1) = 1; 
yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx; zp1(zp1>dimz)=dimz;

%-- get indexes for 18 neighbors
idup = sub2ind(size(phi),yp1,x,z);
iddn = sub2ind(size(phi),ym1,x, z);
idlt = sub2ind(size(phi),y,xm1,z);
idrt = sub2ind(size(phi),y,xp1,z);
idul = sub2ind(size(phi),yp1,xm1,z);
idur = sub2ind(size(phi),yp1,xp1,z);
iddl = sub2ind(size(phi),ym1,xm1,z);
iddr = sub2ind(size(phi),ym1,xp1,z);
idzup = sub2ind(size(phi),y,x,zp1);
idzdown = sub2ind(size(phi),y,x,zm1);
idzlt = sub2ind(size(phi),y,xm1,zp1);
idzrt = sub2ind(size(phi),y, xp1,zm1);
idzul = sub2ind(size(phi),y, xp1,zp1);
idzur = sub2ind(size(phi),y, xm1,zm1);
idzdl = sub2ind(size(phi),ym1,x,zp1);
idzdr = sub2ind(size(phi),yp1,x,zm1);
idz2dl = sub2ind(size(phi),yp1, x, zp1);
idz2dr = sub2ind(size(phi), ym1,x, zm1);

%-- get central derivatives of SDF at x,y,z
    phi_x  = (phi(idlt)-phi(idrt))./2; % (l-r)/2
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt); % l-2c+r
    phi_x2 = phi_x.^2;

    phi_y  = -phi(iddn)+phi(idup); % (u-d)/2
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup); % u-2c+d
    phi_y2 = phi_y.^2;

    phi_z  = - phi(idzdown)+ phi(idzup); % (b-f)/2
    phi_zz = phi(idzdown)-2*phi(idx)+phi(idzup); % b-2c+f
    phi_z2 = phi_z.^2;

 
%(ul+dr-ur-dl)/4
    phi_xy = -0.25*phi(iddl)-0.25*phi(idur)+0.25*phi(iddr)+0.25*phi(idul);

%(lf+rb-rf-lb)/4
    phi_xz = -0.25 * phi(idzul)- 0.25 * phi(idzur) + 0.25* phi(idzlt)+ 0.25*phi(idzrt);

%(uf+db-df-ub)/4
    phi_yz = -0.25 * phi(idz2dl)- 0.25 * phi(idz2dr)+ phi(idzdl)+phi(idzdr);

%-- compute curvature (Kappa)
   curvature = (((phi_x2+phi_z2).*phi_yy + (phi_y2+phi_z2).*phi_xx + (phi_x2+phi_y2).*phi_zz - 2*phi_x.*phi_y.*phi_xy - 2*phi_x.*phi_z.*phi_xz - 2*phi_y.*phi_z.*phi_yz)./...
    (phi_x2 + phi_y2 ++phi_z2 +eps).^(3/2)).*(phi_x2 + phi_y2 + +phi_z2).^(1/2); 

    
    
%-- Displays the image with curve superimposed
function phi = showCurveAndPhi(I, phi, color)
    % subplot(numRows, numCols, plotNum)
    subplot(3, 2, 1);
    imagesc(I(:,:,59)); colormap(gray);
    hold on,
    contour(phi(:,:,59), [0 0], color);
    hold off;

    
    subplot(3, 2, 2);
    imshow(phi(:,:,59), []);
    
    subplot(3, 2, 3);
    imagesc(I(:,:,61)); colormap(gray);
    hold on,
    contour(phi(:,:,61), [0 0], color);
    hold off;

    subplot(3, 2, 4);
    imshow(phi(:,:,61), []); 

    subplot(3, 2, 5);
    imagesc(I(:,:,63));colormap(gray);
    hold on,
    contour(phi(:,:,63), [0 0], color);
    hold off;

    subplot(3, 2, 6)
    imshow(phi(:,:,63), []);     
  


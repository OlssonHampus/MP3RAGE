%MP3RAGE main - final version
%Define sequence parameter settings here
ti = [0 960 3821];   %Inversion times (in ms) for S0, S1 and S2 (1st is always 0)
finv = 0.87;         %Global inversion efficiency
tr = 7.45;           %FLASH TR in ms
alow = 3;            %Low flip angle (in degrees)
ahigh = 16;          %High flip angle (in degrees)
def_path = 'C:\Users\'; %Location of MP3RAGE data

%Colors for parameter maps
map_color{1} = 'jet';    %rsquare
map_color{2} = 'gray';   %SPD
map_color{3} = 'gray';   %T1*
map_color{4} = 'gray';   %T1app
map_color{5} = 'gray';   %T1
map_color{6} = 'hot';    %B1+ (fT)

%Nothing to define beyond this point

addpath('C:\Users\Hampus\Documents\MATLAB\NIfTI_20140122\'); %NIfTI package needed
[s1_file,s_path] = uigetfile('*.nii.gz','Select S1 data',...
    def_path, 'MultiSelect', 'on');
[s2_file,~] = uigetfile('*.nii.gz','Select S2 data',...
    s_path, 'MultiSelect', 'on');
[s3_file,~] = uigetfile('*.nii.gz','Select S3 data',...
    s_path, 'MultiSelect', 'on');
file{1} = s1_file;
file{2} = s2_file;
file{3} = s3_file;

[mask_file,mask_path] = uigetfile('*.nii.gz','Select mask data',...
    s_path, 'MultiSelect', 'on');
nii_mask = load_untouch_nii(strcat(mask_path,'\',mask_file));
mask = double(nii_mask.img);
mask(mask~=0) = 1;

x = input('2D or 3D evaluation? [0/1]');
y = input('Show image input? [0/1]');

for ix = 1:3
    nii_s{ix} = load_untouch_nii(strcat(s_path,'\',file{ix}));
    s(:,:,:,ix+1) = double(nii_s{ix}.img);
end
s = s.*mask;
xdim = length(s(:,1,1,1));
ydim = length(s(1,:,1,1));
zdim = length(s(1,1,:,1));

slice = ceil(length(s(1,1,:,1))/1.6);

alow = alow*pi/180;
ahigh = ahigh*pi/180;
s(:,:,:,1) = -finv*s(:,:,:,4)*sin(alow)/sin(ahigh); %Eq. (4)
if x==1
    st1 = s(:,:,:,4);
else
    st1 = squeeze(s(:,:,slice,4));
end
s0 = s(:,:,:,1);
s1 = s(:,:,:,2);
s2 = s(:,:,:,3);

c = [0 1000];
c_s0 = [-100 0];
if y==1
    figure
    rows = 1;
    cols = 4;
    for ix = 1:4
        subplot(rows,cols,ix)
        if ix == 1
            imshow(flip(rot90(s(:,:,slice,ix)),2),c_s0)
        else
            imshow(flip(rot90(s(:,:,slice,ix)),2),c)
        end
        title_string = strcat('{\itS}_',num2str(ix-1));
        title(title_string)
        colorbar
    end
end

fit_model = fittype('spd+(s0-spd)*exp(-ti/t1s)','independent',{'ti'},...
    'coefficients',{'spd','s0','t1s'});
if x==1
    t1s = zeros(xdim,ydim,zdim);
    spd = zeros(xdim,ydim,zdim);
    r2 = zeros(xdim,ydim,zdim);
    flag = zeros(xdim,ydim,zdim);
else
    t1s = zeros(xdim,ydim);
    spd = zeros(xdim,ydim);
    r2 = zeros(xdim,ydim);
    flag = zeros(xdim,ydim);
end

t1start = 1500;
tic
if x==1
    for xix = 1:xdim
        disp(strcat(num2str(xix),'/',num2str(xdim)))
        for yix = 1:ydim
            parfor zix = 1:zdim
                if sum(mask(xix,yix,zix))~=0
                    sv = squeeze(s(xix,yix,zix,1:3));
                    fo = fitoptions('method','NonlinearLeastSquares',...
                        'StartPoint',[s2(xix,yix,zix),s0(xix,yix,zix),t1start],...
                        'Lower',[s2(xix,yix,zix),s0(xix,yix,zix),0],...
                        'Upper',[inf,s0(xix,yix,zix),inf],...
                        'MaxFunEvals',1e3,'MaxIter',1e3);
                    [cf,gof,output] = fit(ti',sv,fit_model,fo);
                    t1s(xix,yix,zix)=cf.t1s;
                    spd(xix,yix,zix) = cf.spd;
                    r2(xix,yix,zix) = gof.rsquare;
                    flag(xix,yix,zix) = output.exitflag;
                end
            end
        end
    end
else
    for xix = 1:xdim
        disp(strcat(num2str(xix),'/',num2str(xdim)))
        parfor yix = 1:ydim
            if sum(mask(xix,yix,slice))~=0
                sv = squeeze(s(xix,yix,slice,1:3));
                fo = fitoptions('method','NonlinearLeastSquares',...
                    'StartPoint',[s2(xix,yix,slice),s0(xix,yix,slice),t1start],...
                    'Lower',[s2(xix,yix,slice),s0(xix,yix,slice),0],...
                    'Upper',[inf,s0(xix,yix,slice),inf],...
                    'MaxFunEvals',1e3,'MaxIter',1e3);
                [cf,gof,output] = fit(ti',sv,fit_model,fo);
                t1s(xix,yix)=cf.t1s;
                spd(xix,yix) = cf.spd;
                r2(xix,yix) = gof.rsquare;
                flag(xix,yix) = output.exitflag;
            end
        end
    end
end
toc

t1app = 2*tr*(spd/alow-st1/ahigh)./(st1*ahigh-spd*alow); %Eq. (5)
t1 = t1s.*(1+t1app*alow^2/(2*tr));                       %Eq. (7)
ft = sqrt(t1app./t1);                                    %Eq. (8)
if x==1
    param_map(:,:,:,1) = r2;
    param_map(:,:,:,2) = spd;
    param_map(:,:,:,3) = t1s;
    param_map(:,:,:,4) = t1app;
    param_map(:,:,:,5) = t1;
    param_map(:,:,:,6) = ft;
else
    param_map(:,:,1) = r2;
    param_map(:,:,2) = spd;
    param_map(:,:,3) = t1s;
    param_map(:,:,4) = t1app;
    param_map(:,:,5) = t1;
    param_map(:,:,6) = ft;
end
map_name{1} = '{\itr}^2'; map_name{2} = '{\itS}_P_D';...
map_name{3} = '{\itT}_1^*'; map_name{4} = '{\itT}_1_,_a_p_p';...
map_name{5} = '{\itT}_1'; map_name{6} = '{\itB}_1^+';
r2_range = [0.95 1.05];
ft_range = [0.4 1.7];
t1_range = [0 4000];
map_range{1} = r2_range; map_range{2} = c; map_range{3} = t1_range;...
map_range{4} = t1_range; map_range{5} = t1_range; map_range{6} = ft_range;

figure('units','normalized','outerposition',[0 0 1 1])
rows = 2;
cols = 5;
for ix = 1:4
    ax(ix) = subplot(rows,cols,ix);
    if ix==1
        imshow(flip(rot90(s(:,:,slice,ix)),2),c_s0)
    else
        imshow(flip(rot90(s(:,:,slice,ix)),2),c)
    end
    title_string = strcat('{\itS}_',num2str(ix-1));
    title(title_string)
    colormap(ax(ix),'gray')
    colorbar
end
for ix = 1:6
    ax(4+ix) = subplot(rows,cols,4+ix);
    if x==1
        imshow(flip(rot90(param_map(:,:,slice,ix)),2),map_range{ix})
    else
        imshow(flip(rot90(param_map(:,:,ix)),2),map_range{ix})
    end
    title(map_name{ix})
    colormap(ax(4+ix),map_color{ix})
    colorbar
end

%Save files
nii = nii_s{1};
[~,s1_file_name] = fileparts(s1_file);
[~,s1_file_name] = fileparts(s1_file_name);
[~,s2_file_name] = fileparts(s2_file);
[~,s2_file_name] = fileparts(s2_file_name);
[~,s3_file_name] = fileparts(s3_file);
[~,s3_file_name] = fileparts(s3_file_name);

if x==0
    nii.hdr.dime.dim = [2 xdim ydim 1 1 1 1 1];
end

nii.img = t1;
save_untouch_nii(nii,strcat(s_path,'maps\',s1_file_name,...
    s2_file_name,s3_file_name,'_T1.nii.gz'));
nii.img = spd;
save_untouch_nii(nii,strcat(s_path,'maps\',s1_file_name,...
    s2_file_name,s3_file_name,'_SPD.nii.gz'));
nii.img = t1s;
save_untouch_nii(nii,strcat(s_path,'maps\',s1_file_name,...
    s2_file_name,s3_file_name,'_T1s.nii.gz'));
nii.img = t1app;
save_untouch_nii(nii,strcat(s_path,'maps\',s1_file_name,...
    s2_file_name,s3_file_name,'_T1app.nii.gz'));
nii.img = ft;
save_untouch_nii(nii,strcat(s_path,'maps\',s1_file_name,...
    s2_file_name,s3_file_name,'_fT.nii.gz'));
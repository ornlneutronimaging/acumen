%% Create Test Data %%

clear all
close all
clc

% wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/TURBINECT_0180';
% wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/TURBINECT_0010';
wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/CT_Aluminum_ironcenter_0090';
% wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/Derek_injec_0040';
fid = fopen([wdir,'/reconstruction/fileinfo_projections.dat'],'r');
nfiles= fread(fid,1,'single');
nimage = fread(fid,1,'single');
n_files_sino = fread(fid,1,'single');
n_files_rec  = fread(fid,1,'single');
nsub_t = 16;
%fseek(fid,4*9,-1);
fclose(fid);

mkdir(wdir,'fbp_sub_images')
mkdir(wdir,'fbp_sub_tv_images')

add_Nx = 0;
if mod(nimage,2)==0
    add_Nx = 1;
end

% get clim %
rtmp2 = nimage/n_files_sino;
rtmp = nimage/n_files_rec;

figure,
for j_slice = 1:nimage

    js = ceil(j_slice/rtmp2);
    offset = (j_slice-floor((js-1)*rtmp)-1)*(nimage+add_Nx)*(nimage+add_Nx);

    fid = fopen([wdir,'/reconstruction/fbp_sub',num2str(js),'.dat'],'r');
    fseek(fid,4*offset,-1);
    fbp_rec = reshape(fread(fid,(nimage+add_Nx)*(nimage+add_Nx),'single'),nimage+add_Nx,nimage+add_Nx);
    fclose(fid);

    js = ceil(j_slice/rtmp2);
    offset = (j_slice-floor((js-1)*rtmp2)-1)*(nimage+add_Nx)*(nimage+add_Nx);

    fid = fopen([wdir,'/reconstruction/fbp_sub_tv',num2str(js),'.dat'],'r');
    fseek(fid,4*offset,-1);
    fbp_tv = reshape(fread(fid,(nimage+add_Nx)*(nimage+add_Nx),'single'),nimage+add_Nx,nimage+add_Nx);
    fclose(fid);

    
    subplot(1,2,1)
    imagesc(fbp_rec)
    title(['Slice: ',num2str(j_slice)],'FontSize',16)
    axis off equal
    colormap gray
    subplot(1,2,2)
    imagesc(fbp_tv)
    title(['Slice: ',num2str(j_slice)],'FontSize',16)
    axis off equal
    colormap gray
    drawnow
    imwrite(fbp_rec,[wdir,'/fbp_sub_images/fbp_sub',num2str(j_slice),'.tiff'],'tiff');
    imwrite(fbp_tv,[wdir,'/fbp_sub_tv_images/fbp_sub_tv',num2str(j_slice),'.tiff'],'tiff');
end
close all

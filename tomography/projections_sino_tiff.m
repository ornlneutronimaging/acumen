%% Create Test Data %%

clear all
close all
clc

% wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/TURBINECT_0180';
% wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/TURBINECT_0010';
% wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/CT_Aluminum_ironcenter_0090';
wdir = '/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data/Derek_injec_0040';

mkdir(wdir,'projections_images')
mkdir(wdir,'sinogram_images')
fid = fopen([wdir,'/reconstruction/fileinfo_projections.dat'],'r');
nfiles= fread(fid,1,'single');
nimage = fread(fid,1,'single');
n_files_sino = fread(fid,1,'single');
n_files_rec  = fread(fid,1,'single');
fclose(fid);

fid = fopen([wdir,'/reconstruction/theta.dat'],'r');
theta= fread(fid,nfiles,'single');
fclose(fid);

add_Nx = 0;
if mod(nimage,2)==0
    add_Nx = 1;
end


rtmp2 = nimage/n_files_sino;
rtmp = nfiles/n_files_sino;

figure,
for j_slice = 1:nfiles
    js = ceil(j_slice/rtmp);
    offset = (j_slice-floor((js-1)*rtmp)-1)*nimage*nimage;
    
    fid = fopen([wdir,'/reconstruction/projections',num2str(js),'.dat'],'r');
    fseek(fid,4*offset,-1);
    proj = reshape(fread(fid,nimage*nimage,'single'),nimage,nimage);
    fclose(fid);
    imagesc(proj)
    title(['Slice: ',num2str(j_slice)],'FontSize',16)
    axis off 
    colormap gray
    drawnow
    imwrite(proj,[wdir,'/projections_images/projections',num2str(j_slice),'.tiff'],'tiff');
end
for j_slice = 1:nimage
    js = ceil(j_slice/rtmp2);
    offset = (j_slice-floor((js-1)*rtmp2)-1)*(nimage+add_Nx)*nfiles;
    
    fid = fopen([wdir,'/reconstruction/sinograms',num2str(js),'.dat'],'r');
    fseek(fid,4*offset,-1);
    sino = reshape(fread(fid,(nimage+add_Nx)*nfiles,'single'),nimage+add_Nx,nfiles);
    fclose(fid);
    imagesc(sino)
    title(['Slice: ',num2str(j_slice)],'FontSize',16)
    axis off
    colormap gray
    drawnow
    imwrite(sino,[wdir,'/sinogram_images/sinogram',num2str(j_slice),'.tiff'],'tiff');
end
close all

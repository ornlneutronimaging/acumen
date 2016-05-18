program main
	
	!*****************************************************************************
	!
	!  Purpose:
	!
	!    Create corrected HPC projection file from raw projections
	!	 Everything written in single precision since fits files are single
	!
	!    Designed to run on one node on labtop
	!*****************************************************************************	

	USE MPI
	
	implicit none
	

    integer ::  rank, ierr, nprocs,readwrite,blocksize,naxes(2),nfound,group,firstpix
	integer ::  ctheta1,ctheta2,nfiles,nfiles_bg,j,js,tmpint,stat,unit,nimage
	integer ::  jtheta_s,jtheta_e,n_files_sino,n_files_rec
	integer ::  fh_theta,fh_pro,fh_info,fh_bg,add_N,add_Nx,nfiles_df
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer (kind = MPI_OFFSET_KIND) :: offset, empty = 0
	real :: r,nullval,start_time, finish_time,lambda,mu,eps=1e-16
	real, dimension(:), allocatable :: f_theta,f_theta_tmp 
	real, dimension(:,:), allocatable :: image_buffer,bg_image,projection_data,df_image 
	character(LEN=100), dimension(:), allocatable :: projectionFileNames,bgFileNames,dfFileNames
	character(len=150) :: c_string,root_dir,exp_name
	character(len=1) :: d_string
	logical anynull
	
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
	
	

	
	if (nprocs.GT.1) then
		write(*,*) " Designed to run for one processor"
	else

		! find pwd and set up files for sinograms and information !
		call cpu_time(start_time)
		
		
		call get_environment_variable("ROOT_CT_DIR", root_dir)
		call get_environment_variable("EXP_CT_NAME", exp_name)
		
		! Create directorys for data !

		
		
		ctheta1 = len_trim(exp_name)+11
		ctheta2 = len_trim(exp_name)+11+4
		
		call system("mkdir -p " // trim(root_dir) // "/" // trim(exp_name) //  "/reconstruction")
		call system("ls " // trim(root_dir) // "/" // trim(exp_name) // "/projections > fileContents.txt")
		
		open(31,FILE='fileContents.txt',action="read")
		!  Determine number of projections !
		j = 0
		do
			read(31,FMT='(a)',iostat=tmpint) r
			if (tmpint/=0) EXIT
			j = j+1
		end do

		nfiles = j
		write(*,*) "Number of raw projections: " , nfiles
		allocate(projectionFileNames(nfiles))
		allocate(f_theta(nfiles))
		rewind(31)
		do j=1,nfiles
			read(31,'(a)') projectionFileNames(j)
			c_string = trim(projectionFileNames(j))
			f_theta(j) = 0.0
			d_string = c_string(ctheta1:ctheta1)
			read(d_string,*,iostat=stat) tmpint
			f_theta(j) = f_theta(j)+100.0*tmpint
			d_string = c_string(ctheta1+1:ctheta1+1)
			read(d_string,*,iostat=stat) tmpint
			f_theta(j) = f_theta(j)+10.0*tmpint
			d_string = c_string(ctheta1+2:ctheta1+2)
			read(d_string,*,iostat=stat) tmpint
			f_theta(j) = f_theta(j)+1.0*tmpint

			d_string = c_string(ctheta2:ctheta2)
			read(d_string,*,iostat=stat) tmpint
			f_theta(j) = f_theta(j)+0.1*tmpint
			d_string = c_string(ctheta2+1:ctheta2+1)
			read(d_string,*,iostat=stat) tmpint
			f_theta(j) = f_theta(j)+0.01*tmpint
			d_string = c_string(ctheta2+2:ctheta2+2)
			read(d_string,*,iostat=stat) tmpint
			f_theta(j) = f_theta(j)+0.001*tmpint
		enddo
		close(31)
	
		call system('rm fileContents.txt')
		
		! Write HPC theta file !
		fh_theta = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //  & 
			"/reconstruction/" //"theta.dat",MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh_theta, ierr)
		call MPI_FILE_SET_VIEW(fh_theta, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_AT(fh_theta, offset, f_theta, nfiles, MPI_REAL4, status, ierr)
		call MPI_FILE_CLOSE(fh_theta, ierr)
	
		

	

		call system("ls " // trim(root_dir) // "/" // trim(exp_name) // "/OB > fileContents.txt")
		open(41,FILE='fileContents.txt',action="read")

		!  Determine number of background projections !
		j = 0
		do
			read(41,FMT='(a)',iostat=tmpint) r
			if (tmpint/=0) EXIT
			j = j+1
		end do

		nfiles_bg = j
		write(*,*) "Number of background projections: " , nfiles_bg
		allocate(bgFileNames(nfiles_bg))
		rewind(41)
		do j=1,nfiles_bg
			read(41,'(a)') bgFileNames(j)
		enddo
		close(41)
		call system('rm fileContents.txt')
		
		call system("ls " // trim(root_dir) // "/" // trim(exp_name) // "/DF > fileContents.txt")
		open(41,FILE='fileContents.txt',action="read")

		!  Determine number of background projections !
		j = 0
		do
			read(41,FMT='(a)',iostat=tmpint) r
			if (tmpint/=0) EXIT
			j = j+1
		end do

		nfiles_df = j
		write(*,*) "Number of dark field projections: " , nfiles_df
		allocate(dfFileNames(nfiles_df))
		rewind(41)
		do j=1,nfiles_df
			read(41,'(a)') dfFileNames(j)
		enddo
		close(41)
		call system('rm fileContents.txt')
	

			
	    group=1
	    firstpix=1
	    nullval=-999
		call ftgiou(unit,stat)

		readwrite = 0
		stat = 0
		call ftgiou(unit,stat)
		call ftopen(unit,trim(root_dir) // "/" // trim(exp_name) // "/OB/" &
			//trim(bgFileNames(1)),readwrite,blocksize,stat)
		call ftgknj(unit,'NAXIS',1,2,naxes,nfound,stat)
		call ftclos(unit, stat)
		call ftfiou(unit, stat)
	
		nimage = naxes(1)
		allocate(image_buffer(nimage,nimage))
		allocate(bg_image(nimage,nimage))
		allocate(df_image(nimage,nimage))
		allocate(projection_data(nimage,nimage))

	
		bg_image(:,:) = 0
		
		write(*,*) "Determining Average Background" 
		do j=1,nfiles_bg
			stat = 0
			call ftgiou(unit,stat)
			call ftopen(unit,trim(root_dir) // "/" // trim(exp_name) // "/OB/" & 
				//trim(bgFileNames(j)),readwrite,blocksize,stat)
			call ftgpve(unit,group,firstpix,nimage*nimage,nullval,image_buffer,anynull,stat)
			call ftclos(unit, stat)
			call ftfiou(unit, stat)
			bg_image(:,:) = bg_image(:,:)+image_buffer(:,:)
		end do
	
		bg_image(:,:) = bg_image(:,:)/nfiles_bg
		
		df_image(:,:) = 0
		
		write(*,*) "Determining Average Dark Field" 
		do j=1,nfiles_df
			stat = 0
			call ftgiou(unit,stat)
			call ftopen(unit,trim(root_dir) // "/" // trim(exp_name) // "/DF/" & 
				//trim(dfFileNames(j)),readwrite,blocksize,stat)
			call ftgpve(unit,group,firstpix,nimage*nimage,nullval,image_buffer,anynull,stat)
			call ftclos(unit, stat)
			call ftfiou(unit, stat)
			df_image(:,:) = df_image(:,:)+image_buffer(:,:)
		end do
	
		df_image(:,:) = df_image(:,:)/nfiles_df
		
		add_Nx= 0
		if (mod(nimage,2) .EQ. 0) then
		    add_Nx= 1
		end if
		
		n_files_sino = ceiling(real(4.0*(nimage+add_Nx)*nfiles*nimage)/8.0e9)
		n_files_rec  = ceiling(real(4.0*(nimage+add_Nx)*(nimage+add_Nx)*nimage)/8.0e9)
		r = real(nfiles)/real(n_files_sino)
		! If nfiles .GT. 512 this hits a MPI_ERR_QUOTA on my Mac.  Rescale to 512 !

		do j=1,n_files_sino
			
			write(c_string,*) j
			
			fh_pro = 0
			offset = 0
			call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) // &
				"/reconstruction/" //"projections" // trim(adjustl(c_string)) //".dat", & 
				MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh_pro, ierr)
	
		    call MPI_FILE_SET_VIEW(fh_pro, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
			
			do js=1+floor((j-1)*r),min(nfiles,floor(j*r))
			
				stat = 0
				call ftgiou(unit,stat)
				call ftopen(unit,trim(root_dir) // "/" // trim(exp_name) // "/projections/" &
					//trim(projectionFileNames(js)),readwrite,blocksize,stat)
				call ftgpve(unit,group,firstpix,nimage*nimage,nullval,image_buffer,anynull,stat)
				call ftclos(unit, stat)
				call ftfiou(unit, stat)
				projection_data(:,:) = -log(image_buffer(:,:)/max(bg_image(:,:),1.0))
				!projection_data(:,:) = -log(max(image_buffer(:,:)-df_image(:,:),eps)/max(bg_image(:,:)-df_image(:,:),1.0))
				offset = (js-ceiling((j-1)*r)-1)*nimage*nimage
				call MPI_FILE_WRITE_AT(fh_pro, offset,projection_data,nimage*nimage, MPI_REAL4, status, ierr)
				write(*,*) "Processing data ",trim(projectionFileNames(js))
			end do
		
			call MPI_FILE_CLOSE(fh_pro, ierr)
			
		end do
		! Write HPC fileinfo !
		!write(*,*),nfiles,nimage,nsub,nsub_t,f_theta
		fh_info = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) // &
			"/reconstruction/" //"fileinfo_projections.dat",MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh_info, ierr)
		call MPI_FILE_SET_VIEW(fh_info, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_AT(fh_info, offset, real(nfiles), 1, MPI_REAL4, status, ierr)
		offset = 1
		call MPI_FILE_WRITE_AT(fh_info, offset, real(nimage), 1, MPI_REAL4, status, ierr)
		offset = 2
		call MPI_FILE_WRITE_AT(fh_info, offset, real(n_files_sino), 1, MPI_REAL4, status, ierr)
		offset = 3
		call MPI_FILE_WRITE_AT(fh_info, offset, real(n_files_rec), 1, MPI_REAL4, status, ierr)
		call MPI_FILE_CLOSE(fh_info, ierr)
		
		call cpu_time(finish_time)
		write(*,*) "Total write HPC projection data: ",finish_time-start_time
		
		deallocate(image_buffer,bg_image,projection_data,projectionFileNames,f_theta)

	endif
	call MPI_FINALIZE(ierr)
end






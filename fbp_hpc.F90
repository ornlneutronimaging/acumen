module Tomography
	
	integer, public :: ncg,ninnner,nouter,nimage,ntheta,nfilt,nstart,wpu,add_Nx
	real, public :: pi,midpoint,lambda,mu
    real, public, dimension(:), allocatable :: theta,filter
	real, public, dimension(:,:), allocatable :: sinogramCS,cmask,img_rec
	real, public, dimension(:,:), allocatable :: uk,ukp,dk_x,dk_y,bk_x,bk_y,ifkt,rhs,iukt,r,p,irpt,Ap,f,fk,ruk,rp

contains
	
	subroutine int_radon()
		integer :: jx,halfFilterSize
		real :: rad_n,rtmp
		
	    midpoint = (nimage+1)/2
		pi = 4.0*atan(1.0)
	
	    nfilt = 2*(nimage)+1
	    nstart = nimage
	    halfFilterSize = nimage

	
		allocate(filter(nfilt))
	    allocate(theta(ntheta))
		allocate(cmask(nimage,nimage))
		allocate(img_rec(nimage,nimage))
		allocate(sinogramCS(nimage,ntheta))
		allocate(uk(nimage,nimage))
		allocate(ukp(nimage,nimage))
		allocate(dk_x(nimage,nimage))
		allocate(dk_y(nimage,nimage))
		allocate(bk_x(nimage,nimage))
		allocate(bk_y(nimage,nimage))
		allocate(ifkt(nimage,nimage))
		allocate(rhs(nimage,nimage))
		allocate(iukt(nimage,nimage))
		allocate(r(nimage,nimage))
		allocate(p(nimage,nimage))
		allocate(irpt(nimage,nimage))
		allocate(Ap(nimage,nimage))
		allocate(f(nimage,ntheta))
		allocate(fk(nimage,ntheta))
		allocate(ruk(nimage,ntheta))
		allocate(rp(nimage,ntheta))
		
		
	    !$omp barrier
	    !$omp master
		!$acc enter data create(filter,cmask,img_rec,sinogramCS,uk,ukp,dk_x,dk_y,bk_x,bk_y,ifkt,rhs,iukt,r,p,irpt,Ap,f,fk,ruk,rp)
	    !$omp end master
	    !$omp barrier
		
		!$acc parallel loop gang vector collapse(1) present(filter)
		do jx=1,nfilt,1
			filter(jx) = -2.0*(nimage)/((pi**2)*((4*(jx-1-halfFilterSize)**2) -1))
		enddo
		
		rad_n = real(nimage-1)/2.0-wpu
		
		!$acc parallel loop gang vector collapse(2) present(cmask)
		do jx = 1,nimage,1
		    do jy = 1,nimage,1
		        rtmp = 1-2.0*(sqrt(real(jx-(nimage-1.0)/2.0-1.0)**2.0+real(jy-(nimage-1.0)/2.0-1.0)**2.0)-rad_n)/real(wpu)
		        if (rtmp.le.0) then
		            cmask(jx,jy) = 0.0
		        else
					if (rtmp.ge.1) then
		            	cmask(jx,jy) = pi/real(ntheta)
		        	else
		            	cmask(jx,jy) = pi*exp(1-1/(1-exp(1-1/(1-rtmp))))/real(ntheta)
					endif	
				endif
			enddo
		enddo
		

 	   return
 	end subroutine int_radon
	
	subroutine fin_radon()
		!deallocate(filter,cmask,img_rec,sinogramCS,uk,ukp,dk_x,dk_y,bk_x,bk_y,ifkt,rhs,iukt,r,p,irpt,Ap,f,fk,ruk,rp)
	end subroutine fin_radon

	subroutine iradon(img_rec,sinogram)
	
		real, intent(in) :: sinogram(nimage,ntheta)
		real, intent(out) :: img_rec(nimage,nimage)
		real :: xr,yr,ctheta,sum
		integer :: jx,jy,jtheta,y1,y2

	    !! Use Bilinear interpolation to preform Radon transform !!
		
		!$acc parallel loop gang vector collapse(2) present(sinogram,sinogramCS,filter)
	    do jtheta=1,ntheta,1
			do jy=1,nimage,1
				sum = 0
	        	do jx=1,nimage,1	
					jk = jy-jx+1+nstart
	                if ((jk>1) .AND. (jk < nfilt+1)) then
	                    sum = sum+sinogram(jx,jtheta)*filter(jk)
	                end if
	            enddo
				sinogramCS(jy,jtheta) = sum
	        enddo
	    enddo
		
		!$acc parallel loop gang vector collapse(2) present(sinogramCS,theta,cmask,img_rec)
	    
		do jy=1,nimage,1
        	do jx=1,nimage,1
				sum = 0	
				do jtheta=1,ntheta,1
					ctheta = pi*(90.0-theta(jtheta))/180.0
	                yr = sin(ctheta)*(jx-midpoint) + cos(ctheta)*(jy-midpoint)+midpoint
	                y1 = floor(yr)
	                y2 = y1+1
	                if ((y1>1) .AND. (y2 < nimage+1)) then
	                    sum  = sum+sinogramCS(y1,jtheta)*(y2-yr)+sinogramCS(y2,jtheta)*(yr-y1)
	                end if
	            enddo
				img_rec(jx,jy) = cmask(jx,jy)*sum
	        enddo
	    enddo
		
	   return
	end subroutine iradon
	
	
	subroutine get_slice(img_rec_fbp)
	
		real, intent(out) :: img_rec_fbp(nimage,nimage)
		
	    !$omp barrier
	    !$omp master
		!$acc update device(f) 
	    !$omp end master
	    !$omp barrier
	
		call iradon(uk,f)
		
	    !$omp barrier
	    !$omp master
		!$acc update host(uk) 
	    !$omp end master
	    !$omp barrier
		
		img_rec_fbp = uk
		
	    return
	end subroutine get_slice
	


end module Tomography

program main 
	
	!*****************************************************************************
	!
	!  Purpose:
	!
	!    Do Filtered Back Projection from corrected HPC sinogram file 
	!	 Everything written in single precision since fits files are single
	!
	!*****************************************************************************	
	
	use MPI
	use Tomography
	
	implicit none

	integer ierr, rank, nprocs,request_node_id,head_node_id,job_id,tag,jtheta_s,jtheta_e
	integer ::  nsub,nsub_t,nfiles,stat,fh_theta,fh_info,n_files_sino,n_files_rec
	integer ::  jstart,jend,jx,jy,js,jxs,jys,ind,cjx,cjy,x1,y1,x2,y2,jlist,nimagef
	integer, dimension(:), allocatable :: fh_rec,fh_sino
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer (kind = MPI_OFFSET_KIND) :: offset, empty = 0
	real :: rtmp2,rtmp,tm,cmx_0,cmy_0,cmx_180,cmy_180,best_shift,best_theta,midpointx,midpointy
	real :: register(27),xr,yr,start_time, finish_time
	real, dimension(:), allocatable ::f_theta
	real, dimension(:,:), allocatable :: sino,img_rec_fbp
	character(len=150) :: root_dir,exp_name,c_string

	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) 
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	head_node_id = 0
	
	if (nprocs .EQ. 1) then
		write(*,*),"Doing nothing, nproc must be greater than 2 for simple master-slave load balancing"
	else

		! Read HPC fileinfo !
		call get_environment_variable("ROOT_CT_DIR", root_dir)
		call get_environment_variable("EXP_CT_NAME", exp_name)
	
		fh_info = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"fileinfo_projections.dat", &
			MPI_MODE_RDWR, MPI_INFO_NULL, fh_info, ierr)
	
		call MPI_FILE_SET_VIEW(fh_info, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
		
		if (rank.EQ.head_node_id) then
			call cpu_time(start_time)
			write(*,*) 'Updating from file: ', trim(root_dir) // "/" // trim(exp_name) // "/input_parm.txt"
			open(99, file=trim(root_dir) // "/" // trim(exp_name) // "/input_parm.txt",&
					status='old', action='read', position='rewind')
			read(99,*) nsub
			read(99,*) nsub_t
			read(99,*) ninnner
			read(99,*) nouter
			read(99,*) ncg
			read(99,*) lambda
			read(99,*) mu
			read(99,*) wpu
			close(99) 
			
			offset = 2
			call MPI_FILE_WRITE_AT(fh_info, offset, real(nsub), 1, MPI_REAL4, status, ierr)
			offset = 3
			call MPI_FILE_WRITE_AT(fh_info, offset, real(nsub_t), 1, MPI_REAL4, status, ierr)
			offset = 4
			call MPI_FILE_WRITE_AT(fh_info, offset, real(ninnner), 1, MPI_REAL4, status, ierr)
			offset = 5
			call MPI_FILE_WRITE_AT(fh_info, offset, real(nouter), 1, MPI_REAL4, status, ierr)
			offset = 6
			call MPI_FILE_WRITE_AT(fh_info, offset, real(ncg), 1, MPI_REAL4, status, ierr)
			offset = 7
			call MPI_FILE_WRITE_AT(fh_info, offset, lambda, 1, MPI_REAL4, status, ierr)
			offset = 8
			call MPI_FILE_WRITE_AT(fh_info, offset, mu, 1, MPI_REAL4, status, ierr)
			offset = 9
			call MPI_FILE_WRITE_AT(fh_info, offset, real(wpu), 1, MPI_REAL4, status, ierr)
		endif
		
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		
		offset = 0
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		nfiles = int(rtmp)
		offset = 1
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		nimagef = int(rtmp)
		offset = 2
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		nsub = int(rtmp)
		offset = 3
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		nsub_t = int(rtmp)
		offset = 4
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		ninnner = int(rtmp)
		offset = 5
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		nouter = int(rtmp)
		offset = 6
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		ncg = int(rtmp)
		offset = 7
		call MPI_FILE_READ_AT(fh_info, offset, lambda, 1, MPI_REAL4, status, ierr)

		offset = 8
		call MPI_FILE_READ_AT(fh_info, offset, mu, 1, MPI_REAL4, status, ierr)

		offset = 9
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		wpu = int(rtmp)
		offset = 10
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		n_files_sino = int(rtmp)
		offset = 11
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		n_files_rec = int(rtmp)
		
		
		
		!  Do at full resolution for just FPB !
		nsub = 1
		nsub_t = 1
		
		allocate(fh_rec(n_files_rec))
		allocate(fh_sino(n_files_sino))
		
		!! Create HPC Filtered Back Projection File !!
		offset = 0
		!! Open HPC Projection File !!
		do js=1,n_files_rec
			write(c_string,*) js
			fh_rec(js) = 0
		
			call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"fbp" & 
					// trim(adjustl(c_string)) // ".dat", MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh_rec(js), ierr)
			call MPI_FILE_SET_VIEW(fh_rec(js), empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)	
		end do
	
		!! Open HPC Sinogram File !!
		do js=1,n_files_sino
			write(c_string,*) js
			fh_sino(js) = 0
			
			call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"sinograms" & 
					// trim(adjustl(c_string)) // ".dat", MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh_sino(js), ierr)
			call MPI_FILE_SET_VIEW(fh_sino(js), empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)	
			
		end do	
		
		! Read HPC theta file !
		allocate(f_theta(nfiles))
		
		!! Open HPC Theta File !!
		fh_theta = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"theta.dat", &
					MPI_MODE_RDONLY, MPI_INFO_NULL, fh_theta, ierr)
		call MPI_FILE_SET_VIEW(fh_theta, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_READ_AT(fh_theta, offset, f_theta, nfiles, MPI_REAL4, status, ierr)

		jtheta_s = minloc(abs(f_theta), DIM=1)
		jtheta_e = minloc(abs(f_theta - 180.0), DIM=1)
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)

		ntheta = jtheta_e-jtheta_s+1
		add_Nx= 0
		if (mod(nimagef,2) .EQ. 0) then
		    add_Nx= 1
		endif
		nimage = nimagef+add_Nx
	
		allocate(sino(nimagef+add_Nx,nfiles))
	
		call int_radon()
		do jy=1,ntheta,1
		    theta(jy) = f_theta((jy-1)+jtheta_s)
		enddo
		
		if (rank.NE.head_node_id) then
		    !$omp barrier
		    !$omp master
			!$acc enter data copyin(theta)
		    !$omp end master
		    !$omp barrier
		else
			if ((nimage-1)*nsub .NE. nimagef+add_Nx-1) then
				write(*,*), "Warning:  nsub must evenly sub-divide"
			endif
		    write(*,*),"Reconstruction with Resolution: ",nimage,"with",ntheta,"projections"
		endif

		call MPI_FILE_CLOSE(fh_theta, ierr)
	
	
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		! Begin Master Slave Load Balance Sinogram Calculations !
	
		if (rank .EQ. head_node_id) then
			do jlist = 1,nimagef
				call MPI_Recv ( request_node_id, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, status, ierr )
				write(*,*),"Node",request_node_id,"Working on projection",jlist
				tag = 2
				call MPI_Send (jlist, 1, MPI_INTEGER, request_node_id, tag, MPI_COMM_WORLD, ierr )
			enddo
			do jlist = 1,nprocs-1
				call MPI_Recv ( request_node_id, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, status, ierr )
				write(*,*),"Shutting down node",request_node_id
				tag = 2
				call MPI_Send ( -1, 1, MPI_INTEGER, request_node_id, tag, MPI_COMM_WORLD, ierr )
			enddo
		else
			job_id = 1
			allocate(img_rec_fbp(nimage,nimage))
			
			rtmp2 = real(nimagef)/real(n_files_sino)
			rtmp = real(nimagef)/real(n_files_rec)
			do while (job_id .GT. 0)
				tag = 1
				call MPI_Send ( rank, 1, MPI_INTEGER, head_node_id, tag, MPI_COMM_WORLD, ierr )
				tag = 2
				call MPI_Recv ( job_id, 1, MPI_INTEGER, head_node_id, tag,MPI_COMM_WORLD, status, ierr )
				if ( job_id .GT. 0) then
					
					js = ceiling(job_id/rtmp2)
					offset = (job_id-floor((js-1)*rtmp2)-1)*(nimagef+add_Nx)*nfiles
					
					call MPI_FILE_READ_AT(fh_sino(js), offset, sino,(nimagef+add_Nx)*nfiles, MPI_REAL4, status, ierr)
					
					do jy=1,ntheta,1
						do jx=1,nimage,1
						    f(jx,jy) = sino(jx,(jy-1)+jtheta_s)
						enddo
					enddo
					
					call get_slice(img_rec_fbp)

					js = ceiling(job_id/rtmp)
					offset = (job_id-floor((js-1)*rtmp)-1)*nimage*nimage
					call MPI_FILE_WRITE_AT(fh_rec(js), offset, img_rec_fbp,(nimage)*(nimage), MPI_REAL4, status, ierr)
					
				endif
			enddo
			deallocate(img_rec_fbp)
		endif
		
		! Finalize Program !
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		do js=1,n_files_sino
			call MPI_FILE_CLOSE(fh_sino(js), ierr)
		end do
		do js=1,n_files_rec
			call MPI_FILE_CLOSE(fh_rec(js), ierr)
		end do
		
		if (rank.EQ.head_node_id) then
			call cpu_time(finish_time)
			write(*,*) "Total Compute time: ",finish_time-start_time
		else
			call fin_radon()
		endif
	
	endif

	call MPI_FINALIZE(ierr)
	
end	
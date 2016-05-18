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
	
	subroutine radon(sinogram,img)
	
		real, intent(in) :: img(nimage,nimage)
		real, intent(out) :: sinogram(nimage,ntheta)
		real :: xr,yr,ctheta,sum
		integer :: jx,jy,jtheta,x1,y1,x2,y2

		
		!$acc parallel loop gang vector collapse(2) present(img,theta,sinogram)
		do jtheta=1,ntheta,1
			do jy=1,nimage,1
				sum = 0.0
				do jx=1,nimage,1	
					ctheta = pi*(theta(jtheta)-90.0)/180.0
					xr = cos(ctheta)*(jx-midpoint) - sin(ctheta)*(jy-midpoint)+midpoint
					yr = sin(ctheta)*(jx-midpoint) + cos(ctheta)*(jy-midpoint)+midpoint
					y1 = floor(yr)
					y2 = y1+1
					x1 = floor(xr)
					x2 = x1+1
					if ((min(x1,y1)>1) .AND. (max(x2,y2) < nimage+1)) then
						sum = sum+img(x1,y1)*(x2-xr)*(y2-yr)+img(x2,y1)*(xr-x1)*(y2-yr) &
										+img(x1,y2)*(x2-xr)*(yr-y1)+img(x2,y2)*(xr-x1)*(yr-y1)
					end if
				enddo
				sinogram(jy,jtheta) = sum/nimage
			enddo
		enddo

	   return
	end subroutine radon
	
	
	subroutine tvopt(img_rec_tv,img_rec_fbp,jslice)
	
		real, intent(out) :: img_rec_tv(nimage,nimage),img_rec_fbp(nimage,nimage)
		integer, intent(in) :: jslice
		real :: Err_diff,rsold,num,alpha,rsnew,tmpv
		integer :: jncg,jouter,jinner,jx,jxl,jxr,jy,jyl,jyr
		character(len=64) :: int_string,file_string
		
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
		
		
		!$acc parallel loop gang vector collapse(2) present(dk_x,dk_y,bk_x,bk_y)
		do jy=1,nimage,1
			do jx=1,nimage,1
				dk_x(jx,jy) = 0.0
				dk_y(jx,jy) = 0.0
				bk_x(jx,jy) = 0.0
				bk_y(jx,jy) = 0.0
			enddo
		enddo
		
		!$acc parallel loop gang vector collapse(2) present(fk,f)
		do jy=1,ntheta,1
			do jx=1,nimage,1
				fk(jx,jy) = f(jx,jy)
			enddo
		enddo
		
		
		do jouter = 1,nouter,1
		    do jinner = 1,ninnner,1
				!$acc parallel loop gang vector collapse(2) present(ukp,uk)
				do jy=1,nimage,1
					do jx=1,nimage,1
						ukp(jx,jy) = uk(jx,jy)
					enddo
				enddo
				

				call iradon(ifkt,fk)
				
				!$acc parallel loop gang vector collapse(2) present(rhs,ifkt,dk_x,dk_y,bk_x,bk_y)
				do jy=1,nimage,1
					do jx=1,nimage,1
						jxl = jx+1-floor(real(jx+1)/real(nimage+1))* nimage
						jyl = jy+1-floor(real(jy+1)/real(nimage+1))* nimage
						rhs(jx,jy) = mu*ifkt(jx,jy)+lambda*(dk_x(jxl,jy)-bk_x(jxl,jy)-dk_x(jx,jy)+bk_x(jx,jy) & 
										+dk_y(jx,jyl)-bk_y(jx,jyl)-dk_y(jx,jy)+bk_y(jx,jy))
					enddo
				enddo
				

				call radon(ruk,uk)
				call iradon(iukt,ruk)
				
				!$acc parallel loop gang vector collapse(2) present(p,r,rhs,iukt,uk)
				do jy=1,nimage,1
					do jx=1,nimage,1
						jxr = jx+1-floor(real(jx+1)/real(nimage+1))* nimage
		                jyr = jy+1-floor(real(jy+1)/real(nimage+1))* nimage
		                jxl = jx-1+nimage-floor(real(jx-1+nimage)/real(nimage+1))* nimage
		                jyl = jy-1+nimage-floor(real(jy-1+nimage)/real(nimage+1))* nimage
		                r(jx,jy) = rhs(jx,jy)-mu*iukt(jx,jy)-lambda*(4*uk(jx,jy)-uk(jxl,jy)-uk(jxr,jy)-uk(jx,jyl)-uk(jx,jyr))
						p(jx,jy) = r(jx,jy)		
					enddo
				enddo
				
				rsold=0
				!$acc parallel loop gang vector collapse(2) present(r) reduction(+:rsold)
				do jy=1,nimage,1
					do jx=1,nimage,1
						rsold = rsold+r(jx,jy)**2			
					enddo
				enddo
			    
				
		        do jncg=1,ncg,1
					call radon(rp,p)
					call iradon(irpt,rp)

			
					!$acc parallel loop gang vector collapse(2) present(irpt,Ap,p)
					do jy=1,nimage,1
						do jx=1,nimage,1
							jxr = jx+1-floor(real(jx+1)/real(nimage+1))* nimage
		                    jyr = jy+1-floor(real(jy+1)/real(nimage+1))* nimage
		                    jxl = jx-1+nimage-floor(real(jx-1+nimage)/real(nimage+1))* nimage
		                    jyl = jy-1+nimage-floor(real(jy-1+nimage)/real(nimage+1))* nimage
		                    Ap(jx,jy) = mu*irpt(jx,jy)+lambda*(4*p(jx,jy)-p(jxl,jy)-p(jxr,jy)-p(jx,jyl)-p(jx,jyr))
						enddo
					enddo
					
					num = 0
					!$acc parallel loop gang vector collapse(2) present(Ap,p) reduction(+:num)
					do jy=1,nimage,1
						do jx=1,nimage,1
							num = num+p(jx,jy)*Ap(jx,jy)
						enddo
					enddo
					
			
					alpha=rsold/num
					!$acc parallel loop gang vector collapse(2) present(uk,p,r,Ap)
					do jy=1,nimage,1
						do jx=1,nimage,1
							uk(jx,jy) = uk(jx,jy)+alpha*p(jx,jy)
		                    r(jx,jy) =r(jx,jy)-alpha*Ap(jx,jy)
						enddo
					enddo
					
					rsnew = 0
					!$acc parallel loop gang vector collapse(2) present(r) reduction(+:rsnew)
					do jy=1,nimage,1
						do jx=1,nimage,1
		                    rsnew = rsnew+r(jx,jy)**2
						enddo
					enddo
	
					!$acc parallel loop gang vector collapse(2) present(p,r)
					do jy=1,nimage,1
						do jx=1,nimage,1
							p(jx,jy)=r(jx,jy)+(rsnew/rsold)*p(jx,jy)
						enddo
					enddo
					 
					
					rsold=rsnew
		        enddo
				
				!$acc parallel loop gang vector collapse(2) present(uk,dk_x,dk_y,bk_x,bk_y)
				do jy=1,nimage,1
					do jx=1,nimage,1
						jxl = jx-1+nimage-floor(real(jx-1+nimage)/real(nimage+1))* nimage
		                tmpv = uk(jxl,jy)-uk(jx,jy)+bk_x(jx,jy)
						if (tmpv .EQ. 0.0) then
							dk_x(jx,jy)=0.0
						else
		                	dk_x(jx,jy)=max(abs(tmpv)-1/lambda,0.0)*tmpv/abs(tmpv)
						endif
		                bk_x(jx,jy) = tmpv-dk_x(jx,jy)
		                jyl = jy-1+nimage-floor(real(jy-1+nimage)/real(nimage+1))* nimage
		                tmpv = bk_y(jx,jy)+uk(jx,jyl)-uk(jx,jy)
						if (tmpv .EQ. 0.0) then
							dk_y(jx,jy)=0.0
						else
		                	dk_y(jx,jy)=max(abs(tmpv)-1/lambda,0.0)*tmpv/abs(tmpv)
						endif  
		                bk_y(jx,jy) = tmpv-dk_y(jx,jy)
					enddo
				enddo	
				
				Err_diff = 0
				!$acc parallel loop gang vector collapse(2) present(uk,ukp) reduction(+:Err_diff)
				do jy=1,nimage,1
					do jx=1,nimage,1
						Err_diff = Err_diff+abs(uk(jx,jy)-ukp(jx,jy))
					enddo
				enddo
				
			    
			enddo
			

			call radon(ruk,uk) 

			!$acc parallel loop gang vector collapse(2) present(fk,f,ruk)
			do jy=1,ntheta,1
				do jx=1,nimage,1
					fk(jx,jy) = fk(jx,jy)+f(jx,jy)-ruk(jx,jy)
				enddo
			enddo
			
			
			
			if (mod(jouter,10) .EQ. 0) then
				write (*,*) "Slice:",jslice,"Loop: ", jouter, "delta:", Err_diff
			end if
		enddo
		
	    !$omp barrier
	    !$omp master
		!$acc update host(uk) 
	    !$omp end master
	    !$omp barrier

	 	img_rec_tv = uk
	    return
	end subroutine tvopt
	


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
	integer ::  nsub,nsub_t,nfiles,stat,fh_theta,fh_fbp_sub,fh_fbp_tv,fh_info,n_files_sino,n_files_rec
	integer ::  jstart,jend,jx,jy,jxs,jys,ind,cjx,cjy,x1,y1,x2,y2,jlist,nimagef,js
	integer, dimension(:), allocatable :: fh_sino
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer (kind = MPI_OFFSET_KIND) :: offset, empty = 0
	real :: rtmp,tm,cmx_0,cmy_0,cmx_180,cmy_180,best_shift,best_theta,midpointx,midpointy
	real :: register(27),xr,yr,start_time, finish_time
	real, dimension(:), allocatable ::f_theta
	real, dimension(:,:), allocatable :: sino,img_rec_tv,img_rec_fbp
	character(len=150) :: root_dir,exp_name,c_string

	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) 
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	head_node_id = 0
	
	if (nprocs .NE. 1) then
		write(*,*),"Doing nothing, nproc must only one"
	else

		! Read HPC fileinfo !
		call get_environment_variable("ROOT_CT_DIR", root_dir)
		call get_environment_variable("EXP_CT_NAME", exp_name)
	
		fh_info = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"fileinfo_projections.dat", &
			MPI_MODE_RDWR, MPI_INFO_NULL, fh_info, ierr)
	
		call MPI_FILE_SET_VIEW(fh_info, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
		
		call cpu_time(start_time)
			
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
		
		call MPI_FILE_CLOSE(fh_info, ierr)
		
		
		allocate(fh_sino(n_files_sino))
		
		!! Create HPC TV Filtered Back Projection File !!
		fh_fbp_tv = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"fbp_tv_slice.dat", &
					MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh_fbp_tv, ierr)
		call MPI_FILE_SET_VIEW(fh_fbp_tv, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)	
		
		!! Create HPC Filtered Back Projection File for sub resolution !!
		fh_fbp_sub = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"fbp_sub_slice.dat", &
					MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh_fbp_sub, ierr)
		call MPI_FILE_SET_VIEW(fh_fbp_sub, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)	
	
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

		ntheta = floor(real(jtheta_e-jtheta_s+1)/real(nsub_t))
		add_Nx= 0
		if (mod(nimagef,2) .EQ. 0) then
		    add_Nx= 1
		endif
		nimage = floor(real(nimagef+add_Nx-1)/real(nsub))+1
	
		allocate(sino(nimagef+add_Nx,nfiles))
	
		call int_radon()
		do jy=1,ntheta,1
		    theta(jy) = f_theta((jy-1)*nsub_t+jtheta_s)
		enddo
		

	    !$omp barrier
	    !$omp master
		!$acc enter data copyin(theta)
	    !$omp end master
	    !$omp barrier


		job_id = nimage/2
		allocate(img_rec_tv(nimage,nimage))
		allocate(img_rec_fbp(nimage,nimage))

		offset = (job_id-1)*(nimagef+add_Nx)*nfiles
		call MPI_FILE_READ_AT(fh_sino, offset, sino,(nimagef+add_Nx)*nfiles, MPI_REAL4, status, ierr)
		
		do jy=1,ntheta,1
			do jx=1,nimage,1
			    f(jx,jy) = sino((jx-1)*nsub+1,(jy-1)*nsub_t+jtheta_s)
			enddo
		enddo
		
		call tvopt(img_rec_tv,img_rec_fbp,job_id)

		
		offset = 0
		call MPI_FILE_WRITE_AT(fh_fbp_tv, offset, img_rec_tv,(nimage)*(nimage), MPI_REAL4, status, ierr)
		call MPI_FILE_WRITE_AT(fh_fbp_sub, offset, img_rec_fbp,(nimage)*(nimage), MPI_REAL4, status, ierr)
				
		deallocate(img_rec_tv,img_rec_fbp)
	
		
		! Finalize Program !
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		do js=1,n_files_sino
			call MPI_FILE_CLOSE(fh_sino(js), ierr)
		end do
		call MPI_FILE_CLOSE(fh_fbp_tv, ierr)
		call MPI_FILE_CLOSE(fh_fbp_sub, ierr)
		
	endif

	call MPI_FINALIZE(ierr)
	
end	
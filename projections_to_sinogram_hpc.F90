program main 
	
	!*****************************************************************************
	!
	!  Purpose:
	!
	!    Create corrected HPC sinogram file from raw projections
	!	 Everything written in single precision since fits files are single
	!
	!*****************************************************************************	
	
	use MPI
	
	implicit none

	integer ierr, rank, nprocs,request_node_id,head_node_id,job_id,tag,n_sp,n_files_sino,n_files_rec
	integer ::  nsub,nsub_t,nfiles,stat,fh_theta,fh_info,nimage,jtheta_s,jtheta_e
	integer, dimension(:), allocatable :: fh_pro,fh_sino
	integer ::  add_Nx,jstart,jend,jx,jy,jxs,jys,js,ind,cjx,cjy,x1,y1,x2,y2,jlist,jsp,j,ji
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer (kind = MPI_OFFSET_KIND) :: offset, empty = 0
	real :: ra,rb,rtmp,rtmp2,best_shift,best_theta,dshift,dtheta,shift_scale
	real :: best_error,perror_diff,error_diff,dF0_dtheta,dF0_dshift,grd_theta,grd_shift
	real :: reg_max,reg_min,xr,yr,start_time, finish_time,midpointx,midpointy
	real :: reg_tmp1,reg_tmp2,tm,cmx_0,cmx_180,n_mag,n_min
	real, dimension(:), allocatable :: f_theta 
	real, dimension(:,:), allocatable :: image_buffer,sino_0,sino_180,sino_a,sino_b 
	real, dimension(:,:), allocatable :: datal,datac,datar,sino,V180,V0
	real, dimension(:,:,:), allocatable :: register
	character(len=150) :: root_dir,exp_name,c_string

	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) 
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	head_node_id = 0
	
	if (nprocs .EQ. 1) then
		write(*,*),"Doing nothing, nproc must be greater than 2 for simple master-slave load balancing"
	else

		if (rank.EQ.head_node_id) then
			call cpu_time(start_time)
		endif
	
		!call get_environment_variable('PWD',root_dir)
		
		call get_environment_variable("ROOT_CT_DIR", root_dir)
		call get_environment_variable("EXP_CT_NAME", exp_name)
		
		
		! Read HPC fileinfo !
		fh_info = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"fileinfo_projections.dat", &
			MPI_MODE_RDONLY, MPI_INFO_NULL, fh_info, ierr)
		call MPI_FILE_SET_VIEW(fh_info, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		nfiles = int(rtmp)
		offset = 1
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		nimage = int(rtmp)
		offset = 2
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		n_files_sino = int(rtmp)
		offset = 3
		call MPI_FILE_READ_AT(fh_info, offset, rtmp, 1, MPI_REAL4, status, ierr)
		n_files_rec = int(rtmp)
		call MPI_FILE_CLOSE(fh_info, ierr)
		
		if (rank.EQ.head_node_id) then
			call cpu_time(start_time)
			write(*,*) "Working on: ", trim(exp_name), " With number of projections: ",nfiles
		endif
		
		allocate(fh_pro(n_files_sino))
		allocate(fh_sino(n_files_sino))
		
		offset = 0
		
		!! Create HPC Sinogram File !!
		!! Open HPC Projection File !!
		do js=1,n_files_sino
			write(c_string,*) js
			fh_sino(js) = 0
			
			call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"sinograms" & 
					// trim(adjustl(c_string)) // ".dat", MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh_sino(js), ierr)
			call MPI_FILE_SET_VIEW(fh_sino(js), empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)	
			
			fh_pro(js) = 0
			call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"projections" & 
					// trim(adjustl(c_string)) // ".dat",MPI_MODE_RDONLY, MPI_INFO_NULL, fh_pro(js), ierr)
			call MPI_FILE_SET_VIEW(fh_pro(js), empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)	
			
		end do
		

	
		!! Open HPC Theta File !!
		fh_theta = 0
		offset = 0
		call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(root_dir) // "/" // trim(exp_name) //"/reconstruction/" //"theta.dat", &
					MPI_MODE_RDONLY, MPI_INFO_NULL, fh_theta, ierr)
		call MPI_FILE_SET_VIEW(fh_theta, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
	
		if (rank .NE. head_node_id) then
			! Read HPC theta file !
			allocate(f_theta(nfiles))
			call MPI_FILE_READ_AT(fh_theta, offset, f_theta, nfiles, MPI_REAL4, status, ierr)
		endif
		call MPI_FILE_CLOSE(fh_theta, ierr)
		
	
		if (rank .NE. head_node_id) then
			
			add_Nx= 0
			if (mod(nimage,2) .EQ. 0) then
			    add_Nx= 1
			endif
			
			! Allocate Working Arrays !
			allocate(image_buffer(nimage,nimage))
			allocate(sino_0(nimage,nimage))
			allocate(sino_180(nimage,nimage))
			allocate(sino_a(nimage,nimage))
			allocate(sino_b(nimage,nimage))
			allocate(datal(nimage,nimage))
			allocate(datac(nimage,nimage))
			allocate(datar(nimage,nimage))
			allocate(sino(nimage+add_Nx,nimage))
			allocate(V0(nimage+add_Nx,nimage))
			allocate(V180(nimage+add_Nx,nimage))
			allocate(register(nimage,nimage,27))
		    !$omp barrier
		    !$omp master
			!$acc enter data create(V0,V180,sino,sino_0,sino_180,datal,datac,datar,image_buffer,register)
		    !$omp end master
		    !$omp barrier

			jtheta_s = minloc(abs(f_theta), DIM=1)
			jtheta_e = minloc(abs(f_theta-180), DIM=1)
			rtmp = real(nfiles)/real(n_files_sino)
			
			if (f_theta(jtheta_s) .EQ. 0.0) then
				j = ceiling(jtheta_s/rtmp)
				offset = (jtheta_s-floor((j-1)*rtmp)-1)*nimage*nimage
				
				call MPI_FILE_READ_AT(fh_pro(j), offset,sino_0,nimage*nimage, MPI_REAL4, status, ierr)
			else
				if (f_theta(jtheta_s) .LE. 0 ) then
					if (jtheta_s .LT. nfiles) then
						j = ceiling(jtheta_s/rtmp)
						offset = (jtheta_s-floor((j-1)*rtmp)-1)*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_a,nimage*nimage, MPI_REAL4, status, ierr)
						
						j = ceiling((jtheta_s+1)/rtmp)
						offset = (jtheta_s-floor((j-1)*rtmp))*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_b,nimage*nimage, MPI_REAL4, status, ierr)
						sino_0(:,:) = (0.0-f_theta(jtheta_s))*(sino_b(:,:)-sino_a(:,:))/(f_theta(jtheta_s+1)-f_theta(jtheta_s))+sino_a(:,:)
					else
						j = ceiling(jtheta_s/rtmp)
						offset = (jtheta_s-floor((j-1)*rtmp)-1)*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_0,nimage*nimage, MPI_REAL4, status, ierr)
					endif
				else
					if (jtheta_s .GT. 1) then
						j = ceiling((jtheta_s-1)/rtmp)
						offset = (jtheta_s-floor((j-1)*rtmp)-2)*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_a,nimage*nimage, MPI_REAL4, status, ierr)
						
						j = ceiling(jtheta_s/rtmp)
						offset = (jtheta_s-floor((j-1)*rtmp)-1)*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_b,nimage*nimage, MPI_REAL4, status, ierr)
						sino_0(:,:) = (0.0-f_theta(jtheta_s-1))*(sino_b(:,:)-sino_a(:,:))/(f_theta(jtheta_s)-f_theta(jtheta_s-1))+sino_a(:,:)
					else
						j = ceiling(jtheta_s/rtmp)
						offset = (jtheta_s-floor((j-1)*rtmp)-1)*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_0,nimage*nimage, MPI_REAL4, status, ierr)
					endif
				endif
			endif

			if (f_theta(jtheta_e) .EQ. 180.0) then
				j = ceiling(jtheta_e/rtmp)
				offset = (jtheta_e-floor((j-1)*rtmp)-1)*nimage*nimage
				call MPI_FILE_READ_AT(fh_pro(j), offset,sino_180,nimage*nimage, MPI_REAL4, status, ierr)
			else
				if (f_theta(jtheta_e) .LE. 180 ) then
					if (jtheta_e .LT. nfiles) then
						j = ceiling(jtheta_e/rtmp)
						offset = (jtheta_e-floor((j-1)*rtmp)-1)*nimage*nimage
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_a,nimage*nimage, MPI_REAL4, status, ierr)
						j = ceiling((jtheta_e+1)/rtmp)
						offset = (jtheta_e-floor((j-1)*rtmp))*nimage*nimage
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_b,nimage*nimage, MPI_REAL4, status, ierr)
						sino_180(:,:) = (180.0-f_theta(jtheta_e))*(sino_b(:,:)-sino_a(:,:))/(f_theta(jtheta_e+1)-f_theta(jtheta_e))+sino_a(:,:)
					else
						j = ceiling(jtheta_e/rtmp)
						offset = (jtheta_e-floor((j-1)*rtmp)-1)*nimage*nimage
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_180,nimage*nimage, MPI_REAL4, status, ierr)
					endif
				else
					if (jtheta_e .GT. 1) then
						j = ceiling((jtheta_e-1)/rtmp)
						offset = (jtheta_e - floor((j-1)*rtmp)-2)*nimage*nimage
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_a,nimage*nimage, MPI_REAL4, status, ierr)
						j = ceiling(jtheta_e/rtmp)
						offset = (jtheta_e - floor((j-1)*rtmp)-1)*nimage*nimage
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_b,nimage*nimage, MPI_REAL4, status, ierr)
						sino_180(:,:) = (180.0-f_theta(jtheta_e-1))*(sino_b(:,:)-sino_a(:,:))/(f_theta(jtheta_e)-f_theta(jtheta_e-1))+sino_a(:,:)
					else
						j = ceiling(jtheta_e/rtmp)
						offset = (jtheta_e -floor((j-1)*rtmp)-1)*nimage*nimage
						call MPI_FILE_READ_AT(fh_pro(j), offset,sino_180,nimage*nimage, MPI_REAL4, status, ierr)
					endif
				endif
			endif
			
		    !$omp barrier
		    !$omp master
			!$acc update device(sino_0,sino_180) 
		    !$omp end master
		    !$omp barrier
			tm = 0.0
			cmx_0 = 0.0
			cmx_180 = 0.0
			!$acc parallel loop gang vector collapse(2) present(sino_0,sino_180) reduction(+:tm,cmx_0,cmx_180)
			do jx=1,nimage
				do jy=1,nimage
					
					tm = tm+sino_0(jx,jy)+sino_180(jx,jy)
					cmx_0 = cmx_0+sino_0(jx,jy)*(2*real(nimage+add_Nx)*real(jx-1)/real(nimage+add_Nx-1)-nimage-add_Nx)
					cmx_180 = cmx_180+sino_180(jx,jy)*(2*real(nimage+add_Nx)*real(jx-1)/real(nimage+add_Nx-1)-nimage-add_Nx)
					
				enddo
			enddo
			
			shift_scale = real(nimage)
			
			!! Find center of rotation !!
			best_shift = -(cmx_0/tm+cmx_180/tm)/(2.0*shift_scale)
			best_theta = 0.0

			midpointy = real(nimage-1)/2.0
			midpointx = real(nimage+add_Nx-1)/2.0
			
			dshift = 1e-4
			dtheta = 1e-4
			n_sp = 6
			n_mag = 1.05
			n_min = 2.0/real(n_sp)
			
			!$acc parallel loop gang vector collapse(2) present(V0,V180,sino_0,sino_180)	
			do jx=1,nimage+add_Nx
				do jy=1,nimage
					xr = min(max(cos(best_theta)*(jx-midpointx-shift_scale*best_shift) - &
									sin(best_theta)*(jy-midpointy)+midpointx,1.0),real(nimage))
		            yr = min(max(sin(best_theta)*(jx-midpointx-shift_scale*best_shift) + &
									cos(best_theta)*(jy-midpointy)+midpointy,1.0),real(nimage))
		            y1 = min(max(floor(yr),1),nimage)
		            y2 = min(max(y1+1,1),nimage)
		            x1 = min(max(floor(xr),1),nimage)
		            x2 = min(max(x1+1,1),nimage)
		           	V0(jx,jy) = (sino_0(x1,y1)*(x2-xr)*(y2-yr)+sino_0(x2,y1)*(xr-x1)*(y2-yr)+ &
									sino_0(x1,y2)*(x2-xr)*(yr-y1)+sino_0(x2,y2)*(xr-x1)*(yr-y1))
		           	V180(jx,jy) = (sino_180(x1,y1)*(x2-xr)*(y2-yr)+sino_180(x2,y1)*(xr-x1)*(y2-yr)+ &
									sino_180(x1,y2)*(x2-xr)*(yr-y1)+sino_180(x2,y2)*(xr-x1)*(yr-y1))
				enddo
			enddo
			
			best_error = 0
			!$acc parallel loop gang vector collapse(2) present(V0,V180) reduction(+:best_error)
			do jx=1,nimage+add_Nx
				do jy=1,nimage
					best_error = best_error +abs(V0(jx,jy)-V180(nimage+add_Nx-jx+1,jy))
				enddo
			enddo

			
			do while(max(dtheta,dshift) .GT. 1e-8)
				perror_diff = best_error
				do jsp=1,n_sp
					!$acc parallel loop gang vector collapse(2) present(V0,V180,sino_0,sino_180)
					do jx=1,nimage+add_Nx
						do jy=1,nimage
							xr = min(max(cos(best_theta+dtheta*(jsp-n_sp/2))*(jx-midpointx-shift_scale*best_shift) - &
											sin(best_theta+dtheta*(jsp-n_sp/2))*(jy-midpointy)+midpointx,1.0),real(nimage))
				            yr = min(max(sin(best_theta+dtheta*(jsp-n_sp/2))*(jx-midpointx-shift_scale*best_shift) + &
											cos(best_theta+dtheta*(jsp-n_sp/2))*(jy-midpointy)+midpointy,1.0),real(nimage))
				            y1 = min(max(floor(yr),1),nimage)
				            y2 = min(max(y1+1,1),nimage)
				            x1 = min(max(floor(xr),1),nimage)
				            x2 = min(max(x1+1,1),nimage)
				           	V0(jx,jy) = (sino_0(x1,y1)*(x2-xr)*(y2-yr)+sino_0(x2,y1)*(xr-x1)*(y2-yr)+ &
											sino_0(x1,y2)*(x2-xr)*(yr-y1)+sino_0(x2,y2)*(xr-x1)*(yr-y1))
				           	V180(jx,jy) = (sino_180(x1,y1)*(x2-xr)*(y2-yr)+sino_180(x2,y1)*(xr-x1)*(y2-yr)+ &
											sino_180(x1,y2)*(x2-xr)*(yr-y1)+sino_180(x2,y2)*(xr-x1)*(yr-y1))
						enddo
					enddo
					error_diff = 0
					!$acc parallel loop gang vector collapse(2) present(V0,V180) reduction(+:error_diff)
					do jx=1,nimage+add_Nx
						do jy=1,nimage
							error_diff = error_diff +abs(V0(jx,jy)-V180(nimage+add_Nx-jx+1,jy))
						enddo
					enddo
					if (error_diff .LT. best_error) then
				        best_theta = best_theta+dtheta*(jsp-n_sp/2)
				        best_error = error_diff
					endif
				enddo
				
				if (best_error .LT. perror_diff) then
			        dtheta = n_mag*dtheta
				else
					dtheta = n_min*dtheta
				endif
				
				perror_diff = best_error
				do jsp=1,n_sp
					!$acc parallel loop gang vector collapse(2) present(V0,V180,sino_0,sino_180)
					do jx=1,nimage+add_Nx
						do jy=1,nimage
							xr = min(max(cos(best_theta)*(jx-midpointx-shift_scale*(best_shift+dshift*(jsp-n_sp/2))) - &
											sin(best_theta)*(jy-midpointy)+midpointx,1.0),real(nimage))
				            yr = min(max(sin(best_theta)*(jx-midpointx-shift_scale*(best_shift+dshift*(jsp-n_sp/2))) + &
											cos(best_theta)*(jy-midpointy)+midpointy,1.0),real(nimage))
				            y1 = min(max(floor(yr),1),nimage)
				            y2 = min(max(y1+1,1),nimage)
				            x1 = min(max(floor(xr),1),nimage)
				            x2 = min(max(x1+1,1),nimage)
				           	V0(jx,jy) = (sino_0(x1,y1)*(x2-xr)*(y2-yr)+sino_0(x2,y1)*(xr-x1)*(y2-yr)+ &
											sino_0(x1,y2)*(x2-xr)*(yr-y1)+sino_0(x2,y2)*(xr-x1)*(yr-y1))
				           	V180(jx,jy) = (sino_180(x1,y1)*(x2-xr)*(y2-yr)+sino_180(x2,y1)*(xr-x1)*(y2-yr)+ &
											sino_180(x1,y2)*(x2-xr)*(yr-y1)+sino_180(x2,y2)*(xr-x1)*(yr-y1))
						enddo
					enddo
					error_diff = 0
					!$acc parallel loop gang vector collapse(2) present(V0,V180) reduction(+:error_diff)
					do jx=1,nimage+add_Nx
						do jy=1,nimage
							error_diff = error_diff +abs(V0(jx,jy)-V180(nimage+add_Nx-jx+1,jy))
						enddo
					enddo
					if (error_diff .LT. best_error) then
				        best_shift = best_shift+dshift*(jsp-n_sp/2)
				        best_error = error_diff
					endif
				enddo
				
				if (best_error .LT. perror_diff) then
			        dshift = n_mag*dshift
				else
					dshift = n_min*dshift
				endif
				!write(*,*),best_error,best_shift,best_theta	
			enddo
			if (rank .EQ. head_node_id+1) then
				write(*,*),"Best theta: ",best_theta, " Best shift: ",best_shift
			endif
		endif
	
	
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		! Begin Master Slave Load Balance Sinogram Calculations !
		rtmp2 = real(nimage)/real(n_files_sino)
		rtmp = real(nfiles)/real(n_files_sino)
		if (rank .EQ. head_node_id) then
			
			do jlist = 1,nfiles
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
			do while (job_id .GT. 0)
				tag = 1
				call MPI_Send ( rank, 1, MPI_INTEGER, head_node_id, tag, MPI_COMM_WORLD, ierr )
				tag = 2
				call MPI_Recv ( job_id, 1, MPI_INTEGER, head_node_id, tag,MPI_COMM_WORLD, status, ierr )
				if ( job_id .GT. 0) then
					j = ceiling(job_id/rtmp)
					offset = (job_id-floor((j-1)*rtmp)-1)*nimage*nimage
					
					call MPI_FILE_READ_AT(fh_pro(j), offset,datac,nimage*nimage, MPI_REAL4, status, ierr)
					if(job_id .GT. 1) then
						j = ceiling((job_id-1)/rtmp)
						offset = (job_id-floor((j-1)*rtmp)-2)*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,datal,nimage*nimage, MPI_REAL4, status, ierr)
					else
						datal(:,:)=datac(:,:)
					endif
					if(job_id .LT. nfiles) then				
						j = ceiling((job_id+1)/rtmp)
						offset = (job_id-floor((j-1)*rtmp))*nimage*nimage
						
						call MPI_FILE_READ_AT(fh_pro(j), offset,datar,nimage*nimage, MPI_REAL4, status, ierr)
					else
						datar(:,:) = datac(:,:)
					endif
					! Do Modified 3D Median Filter with unit box length !
					!call cpu_time(start_time)
					
				    !$omp barrier
				    !$omp master
					!$acc update device(datal,datac,datar) 
				    !$omp end master
				    !$omp barrier
					
					!$acc parallel loop gang vector collapse(2) present(datal,datac,datar,image_buffer,register) 
					do jx=1,nimage
						do jy=1,nimage
							do jxs =-1,1
								do jys =-1,1
									cjx = min(max(jx+jxs,1),nimage)
									cjy = min(max(jy+jys,1),nimage)
									register(jx,jy,3*(jxs+1-1)+jys+2) = datal(cjx,cjy);
									register(jx,jy,3*(jxs+1-1)+jys+2+9) = datac(cjx,cjy);
									register(jx,jy,3*(jxs+1-1)+jys+2+18) = datar(cjx,cjy);
								end do
							end do
							! sort register !
						    do jxs=1, 26
								ra = register(jx,jy,jxs)
								do jys = jxs+1,27
									rb = register(jx,jy,jys)
									if (rb .GT. ra) then
										register(jx,jy,jys) = ra
										register(jx,jy,jxs) = rb
										ra = rb
									end if
								end do
							end do
							image_buffer(jx,jy)  = register(jx,jy,14)
						enddo
					enddo
					
					! Correct Projection Rotation using Bi-linear Interpolation !
					
					!$acc parallel loop gang vector collapse(2) present(sino,image_buffer)
					do jx=1,nimage+add_Nx
						do jy=1,nimage
							xr = min(max(cos(best_theta)*(jx-midpointx-shift_scale*best_shift) - &
											sin(best_theta)*(jy-midpointy)+midpointx,1.0),real(nimage))
				            yr = min(max(sin(best_theta)*(jx-midpointx-shift_scale*best_shift) + &
											cos(best_theta)*(jy-midpointy)+midpointy,1.0),real(nimage))
				            y1 = min(max(floor(yr),1),nimage)
				            y2 = min(max(y1+1,1),nimage)
				            x1 = min(max(floor(xr),1),nimage)
				            x2 = min(max(x1+1,1),nimage)
				           	sino(jx,jy) = (image_buffer(x1,y1)*(x2-xr)*(y2-yr)+image_buffer(x2,y1)*(xr-x1)*(y2-yr)+ &
											image_buffer(x1,y2)*(x2-xr)*(yr-y1)+image_buffer(x2,y2)*(xr-x1)*(yr-y1))
						enddo
					enddo
					
				    !$omp barrier
				    !$omp master
					!$acc update host(sino) 
				    !$omp end master
				    !$omp barrier

					
					! Write HPC Sinogram File !

					do jy=1,nimage
						j = ceiling(jy/rtmp2)
						
						offset = (jy-floor((j-1)*rtmp2)-1)*(nimage+add_Nx)*nfiles+(job_id-1)*(nimage+add_Nx)
						
						call MPI_FILE_WRITE_AT(fh_sino(j), offset, sino(1:nimage+add_Nx,jy),nimage+add_Nx, MPI_REAL4, status, ierr)
					enddo
				endif
			enddo
		endif
		
		! Finalize Program !
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		do js=1,n_files_sino
			call MPI_FILE_CLOSE(fh_pro(js), ierr)
			call MPI_FILE_CLOSE(fh_sino(js), ierr)
		end do
		
		deallocate(fh_pro,fh_sino)
		
		if (rank.EQ.head_node_id) then
			call cpu_time(finish_time)
			write(*,*) "Total Compute time: ",finish_time-start_time
		else
			deallocate(image_buffer,f_theta,sino_0,sino_180,sino_a,sino_b,datal,datac,datar,sino)
		endif
	
	endif

	call MPI_FINALIZE(ierr)
	
end	




!-----------------------------------------------------------------------------!
module BL_sub_variables
	use define_precision_mod
	use type_of_data

	type(real4d),pointer :: mb_FY(:)
	type(int4d),pointer :: mb_indexwall(:)
	integer,pointer :: indexwall(:,:,:,:)

end module BL_sub_variables
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
subroutine allocate_bl_sub_variables
  use global_variables,only : mb_dim,nblocks
  use bl_sub_variables

	allocate (mb_fy(nblocks),mb_indexwall(nblocks))

	do nb=1,nblocks
		allocate( mb_fy(nb)%a4d(mb_dim(nb,1),mb_dim(nb,2),mb_dim(nb,3),3), &
		  mb_indexwall(nb)%a4d(mb_dim(nb,1),mb_dim(nb,2),mb_dim(nb,3),4) )
	enddo

end subroutine allocate_bl_sub_variables

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine deallocate_bl_sub_variables
  use bl_sub_variables
  use global_variables,only : nblocks

	do nb=1,nblocks
		deallocate( mb_fy(nb)%a4d,mb_indexwall(nb)%a4d )
	enddo
	deallocate (mb_fy,mb_indexwall)

end subroutine deallocate_bl_sub_variables

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine cal_vist_bl
	use global_variables
	use define_precision_mod
    use bl_sub_variables
	implicit none
	integer :: nb,i,j,k,m,nmg
	integer :: nbw,iw,jw,kw
	real(prec)    :: fmax,ymax,taow,yw,ro,vsl,vst,vsti,vsto
	real(prec)    :: yplus,omg,Fwake,udif2,Fkelb,um,vm,wm,argexp,yplim
	real(prec),parameter :: Cwk=0.25,Ckelb=0.3,oAplus=1.0/26.0,KCcp=0.02688,kk=0.4, &
	                  logsmall = log(1.0e-24)

!	call allocate_bl_sub_variables  !调试程序结束后把该语句移到预处理中以便节约计算时间

	call cal_wall_info

!	nmg = nmultigrid

	do nb=1,nblocks

		call recast_grid(nb)
		call recast_field(nb)

!		fy => mb_fy(nb)%a4d

		indexwall => mb_indexwall(nb)%a4d

		call allocate_variables_temp_turbulence
  
!		call dfidxyz_center(-1,1,ni,nj,nk,qprv,dudxyz,method,5,2,2)
!		call dfidxyz_center(-1,1,ni,nj,nk,qprv,dvdxyz,method,5,3,3)
!		call dfidxyz_center(-1,1,ni,nj,nk,qprv,dwdxyz,method,5,4,4)

		call dfidxyz_center_1(-1,1,ni,nj,nk,u,dudxyz,method)
		call dfidxyz_center_1(-1,1,ni,nj,nk,v,dvdxyz,method)
		call dfidxyz_center_1(-1,1,ni,nj,nk,w,dwdxyz,method)

		do k=1,nk !-1
			do j=1,nj !-1
				do i=1,ni !-1

					iw    = mb_indexwall(nb)%a4d(i,j,k,1)
					jw    = mb_indexwall(nb)%a4d(i,j,k,2)
					kw    = mb_indexwall(nb)%a4d(i,j,k,3)
					nbw   = mb_indexwall(nb)%a4d(i,j,k,4)

					ro    = r(i,j,k) !qprv(i,j,k,1)
					um    = u(i,j,k) !qprv(i,j,k,2)
					vm    = v(i,j,k) !qprv(i,j,k,3)
					wm    = w(i,j,k) !qprv(i,j,k,4)
					vsl   = visl(i,j,k)
					yw    = mb_distance(nb)%a3d(i,j,k)

					ymax  = mb_FY(nbw)%a4d(iw,jw,kw,1)
					fmax  = mb_FY(nbw)%a4d(iw,jw,kw,2)
					taow  = mb_FY(nbw)%a4d(iw,jw,kw,3)

					udif2 = 1.0 !um*um+vm*vm+wm*wm

					yplus = sqrt(ro*taow*reynolds)*yw/vsl

			        omg   = sqrt((dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				               (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
										   (dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

!					Fwake = ymax*amin1(fmax,Cwk*udif2/fmax)
                    Fwake = ymax*min(fmax,Cwk*udif2/fmax)

!					fkelb = Ckelb*yw/amax1(ymax,1.0e-14)
                    fkelb = Ckelb*yw/max(ymax,1.0e-14)

!					fkelb = amin1(fkelb,1.0e5)
                    fkelb = min(fkelb,1.0e5)

					fkelb = 1.0/(1.0+5.5*(fkelb**6.0))

	                argexp = max(logsmall,-oAplus*yplus)

					yplim  = 1.0 - exp(argexp)

					vsti  = ro*(kk*yw*yplim)**2.0*omg*reynolds

					vsto  = ro*KCcp*Fwake*Fkelb*reynolds

!					vist(i,j,k) = amin1(vsti,vsto)

					vist(i,j,k) = vsto*tanh(vsti/(vsto+1.0e-20))

!					if(vsti<vsto)write(*,*)i,j,k

				enddo
			enddo
		enddo
	call deallocate_variables_temp_turbulence
	enddo

!	call deallocate_bl_sub_variables  !调试程序结束后把该语句移到预处理中以便节约计算时间
	
	return
end subroutine cal_vist_bl

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

subroutine cal_wall_info
!* tgh BL代数湍流模型使用，仅对有限差分适用

  use define_precision_mod
  use global_variables !,only:mb_bc,nblocks,nmultigrid,mb_dim,reynolds
  use bl_sub_variables
	implicit none

	integer :: s_st(3),s_ed(3),s_lr3d(3)
	integer :: s_nd,s_lr,s_fix,bctype,no
	integer :: nb,i,j,k,m,nr,nrmax,iw,jw,kw
	integer :: id,jd,kd
	real(prec)    :: xwc,ywc,zwc,xc,yc,zc,npt,yw,yplus,omg,fy
	real(prec)    :: uw,vw,ww,uvww,sx,sy,sz,volw,sw,vsl,taow,ro
	real(prec)    :: fmax,ymax,uvw2,uvwn,nx,ny,nz,nxyz,utao,argexp,yplim
	real(prec),parameter :: Cwk=0.25,Ckelb=0.3,oAplus=1.0/26.0,KCcp=0.02688,kk=0.4, &
	                  logsmall = log(1.0e-24)

	do nb=1,nblocks

		call recast_grid(nb)
		call recast_field(nb)
!		fy => mb_fy(nb)%a4d

	    call allocate_variables_temp_turbulence

        nrmax = mb_bc(nb)%nregions         
  
!	    call dfidxyz_center(-1,1,ni,nj,nk,qprv,dudxyz,method,5,2,2)
!	    call dfidxyz_center(-1,1,ni,nj,nk,qprv,dvdxyz,method,5,3,3)
!	    call dfidxyz_center(-1,1,ni,nj,nk,qprv,dwdxyz,method,5,4,4)

	    call dfidxyz_center_1(-1,1,ni,nj,nk,u,dudxyz,method)
	    call dfidxyz_center_1(-1,1,ni,nj,nk,v,dvdxyz,method)
	    call dfidxyz_center_1(-1,1,ni,nj,nk,w,dwdxyz,method)

	    do nr = 1,nrmax
            bctype = mb_bc(nb)%bc(nr)%bctype

            if( bctype==2 .or. (bctype>=20 .and. bctype<=29) ) then               

			    do m=1,3
				    s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)   
				    s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)   
				    s_lr3d(m)=mb_bc(nb)%bc(nr)%s_lr3d(m)
			    enddo

			    s_nd  = mb_bc(nb)%bc(nr)%s_nd           
			    s_lr  = mb_bc(nb)%bc(nr)%s_lr           
			    s_fix = mb_bc(nb)%bc(nr)%s_fix 

			    id = (1+s_lr3d(1))/2
			    jd = (1+s_lr3d(2))/2
			    kd = (1+s_lr3d(3))/2
			
			    if(s_nd==1)then
				    do k = s_st(3),s_ed(3)
					do j = s_st(2),s_ed(2)

						iw = s_st(1) - s_lr  !* 离开物面第一层的位置
						i  = s_st(1)         !* 物面的位置

						xwc = x(i,j,k)
						ywc = y(i,j,k)
						zwc = z(i,j,k)

					    nx = kcx(i,j,k)
					    ny = kcy(i,j,k)
					    nz = kcz(i,j,k)

						nxyz = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

					    nx = nx/nxyz
					    ny = ny/nxyz
					    nz = nz/nxyz

						vsl = visl(i,j,k)
						uw  = u(iw,j,k)
						vw  = v(iw,j,k)
						ww  = w(iw,j,k)

						uvw2 = uw*uw+vw*vw+ww*ww

						uvwn = uw*nx+vw*ny+ww*nz

						uvww = sqrt(uvw2-uvwn*uvwn)+1.0e-20

                        xc = x(iw,j,k)
						yc = y(iw,j,k)
						zc = z(iw,j,k)

						yw = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

						taow = vsl*uvww/(yw+1.0e-20)

						if(s_lr==-1)then

							fmax = -1.0e10

							do i = s_st(1),mb_dim(nb,1)

							    no  = 0

						        vsl = visl(i,j,k)

						        ro  = r(i,j,k)

								xc=x(i,j,k)
							    yc=y(i,j,k)
							    zc=z(i,j,k)

!		       	                call get_center_cor(nb,nmg,i,j,k,xc,yc,zc)

								yw    = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

								yplus = sqrt(ro*abs(taow)*reynolds)*yw/vsl

			    	            omg   = sqrt( (dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                              (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
											  (dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

	                            argexp = max(logsmall,-oAplus*yplus)

								yplim  = 1.0 - exp(argexp)

								fy    = yw*omg*yplim !* F_wake

								if(fy>fmax)then
									fmax = fy
									ymax = yw
								endif

								if(fy<0.9*fmax)goto 10

							enddo

10						    mb_FY(nb)%a4d(iw,j,k,1) = ymax ! ?????
							mb_FY(nb)%a4d(iw,j,k,2) = fmax ! ?????
							mb_FY(nb)%a4d(iw,j,k,3) = taow  ! ?????

						else

							fmax = -1.0e10

							do i = s_st(1),1,-1
							    no  = 0

						        vsl = visl(i,j,k) !mb_visl(nb)%a3d(i,j,k)

						        ro  = r(i,j,k) !mb_qprv(nb)%a4d(i,j,k,1)

								xc=x(i,j,k)
							    yc=y(i,j,k)
							    zc=z(i,j,k)

								yw    = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

								yplus = sqrt(ro*abs(taow)*reynolds)*yw/vsl

			    	            omg   = sqrt( (dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                              (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
											  (dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

								fy    = yw*omg*(1.0-exp(-yplus*oAplus))

								if(fy>fmax)then

									fmax = fy
									ymax = yw
									no   = 1
									
								endif

								if(fy<0.9*fmax)goto 20

							enddo

20						    mb_FY(nb)%a4d(iw,j,k,1) = ymax
							mb_FY(nb)%a4d(iw,j,k,2) = fmax
							mb_FY(nb)%a4d(iw,j,k,3) = taow

						endif
					enddo
				    enddo

			    elseif(s_nd==2)then

				    do k = s_st(3),s_ed(3)
					do i = s_st(1),s_ed(1)

						jw = s_st(2)-s_lr
						j  = s_st(2)

						xwc = x(i,j,k)
						ywc = y(i,j,k)
						zwc = z(i,j,k)

					    nx = etx(i,j,k)
					    ny = ety(i,j,k)
					    nz = etz(i,j,k)

						nxyz = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

					    nx = nx/nxyz
					    ny = ny/nxyz
					    nz = nz/nxyz

						vsl = visl(i,j,k)
						uw  = u(i,jw,k)
						vw  = v(i,jw,k)
						ww  = w(i,jw,k)


						uvw2 = uw*uw+vw*vw+ww*ww

						uvwn = uw*nx+vw*ny+ww*nz

						uvww = sqrt(uvw2-uvwn*uvwn)

!		       	        call get_center_cor(nb,nmg,i,s_st(2),k,xc,yc,zc)

                        xc = x(i,jw,k)
						yc = y(i,jw,k)
						zc = z(i,jw,k)

						yw = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

						taow = vsl*uvww/(yw+1.0e-20)

						if(s_lr==-1)then

							fmax = -1.0e10

							do j = s_st(2),mb_dim(nb,2)

						        vsl = visl(i,j,k)

						        ro  = r(i,j,k)

!		       	                call get_center_cor(nb,nmg,i,j,k,xc,yc,zc)

								xc=x(i,j,k)
							    yc=y(i,j,k)
							    zc=z(i,j,k)

								yw    = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

								yplus = sqrt(ro*abs(taow)*reynolds)*yw/vsl

			    	            omg   = sqrt( (dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                              (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
											  (dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

								fy    = yw*omg*(1.0-exp(-yplus/26.0))

								if(fy>fmax)then

									fmax = fy
									ymax = yw
									
								endif

								if(fy<0.9*fmax)goto 30

							enddo

30						    mb_FY(nb)%a4d(i,jw,k,1) = ymax
							mb_FY(nb)%a4d(i,jw,k,2) = fmax
							mb_FY(nb)%a4d(i,jw,k,3) = taow

						else

							fmax = -1.0e10

							do j = s_st(2),1,-1
							    no  = 0

						        vsl =visl(i,j,k) ! mb_visl(nb)%a3d(i,j,k)

						        ro  = r(i,j,k) !mb_qprv(nb)%a4d(i,j,k,1)

!		       	                call get_center_cor(nb,nmg,i,j,k,xc,yc,zc)

								xc=x(i,j,k)
							    yc=y(i,j,k)
							    zc=z(i,j,k)

								yw    = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

								yplus = sqrt(ro*abs(taow)*reynolds)*yw/vsl

			    	            omg   = sqrt( (dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                              (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
											  (dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

								fy    = yw*omg*(1.0-exp(-yplus/26.0))

								if(fy>fmax )then
									fmax = fy
									ymax = yw
									no   = 1
								endif

								if(fy<0.9*fmax)goto 40

							enddo

40						    mb_FY(nb)%a4d(i,jw,k,1) = ymax
							mb_FY(nb)%a4d(i,jw,k,2) = fmax
							mb_FY(nb)%a4d(i,jw,k,3) = taow

						endif
					enddo
				    enddo
			    else
				    do j = s_st(2),s_ed(2)
					do i = s_st(1),s_ed(1)

						kw = s_st(3)-s_lr
						k  = s_st(3)

!                       call facecoor(nb,nmg,i,j,kw,s_nd,xwc,ywc,zwc)

						xwc = x(i,j,k)
						ywc = y(i,j,k)
						zwc = z(i,j,k)

					    nx = ctx(i,j,k)
					    ny = cty(i,j,k)
					    nz = ctz(i,j,k)

!					    nx = mb_sxyz(nb,nmg)%a5d(i+id,j+jd,k+kd,1,s_nd)
!					    ny = mb_sxyz(nb,nmg)%a5d(i+id,j+jd,k+kd,2,s_nd)
!					    nz = mb_sxyz(nb,nmg)%a5d(i+id,j+jd,k+kd,3,s_nd)

						nxyz = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

					    nx = nx/nxyz
					    ny = ny/nxyz
					    nz = nz/nxyz

						vsl = visl(i,j,k)
						uw  = u(i,j,kw)
						vw  = v(i,j,kw)
						ww  = w(i,j,kw)

!						vsl = mb_visl(nb,nmg)%a3d(i,j,k)
!						uw = mb_qprv(nb,nmg)%a4d(i,j,k,2)
!						vw = mb_qprv(nb,nmg)%a4d(i,j,k,3)
!						ww = mb_qprv(nb,nmg)%a4d(i,j,k,4)


						uvw2 = uw*uw+vw*vw+ww*ww

						uvwn = uw*nx+vw*ny+ww*nz

						uvww = sqrt(uvw2-uvwn*uvwn)

!		       	        call get_center_cor(nb,nmg,i,j,s_st(3),xc,yc,zc)

                        xc = x(i,j,kw)
						yc = y(i,j,kw)
						zc = z(i,j,kw)

						yw = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

						taow = vsl*uvww/(yw+1.0e-20)

						if(s_lr==-1)then

							fmax = -1.0e10

							do k = s_st(3),mb_dim(nb,3)

							    no  = 0

						        vsl = visl(i,j,k)

						        ro  = r(i,j,k)

!		       	                call get_center_cor(nb,nmg,i,j,k,xc,yc,zc)

								xc=x(i,j,k)
							    yc=y(i,j,k)
							    zc=z(i,j,k)

								yw    = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

								yplus = sqrt(ro*abs(taow)*reynolds)*yw/vsl

			    	            omg   = sqrt( (dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                              (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
											  (dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

								fy    = yw*omg*(1.0-exp(-yplus/26.0))

								if(fy>fmax)then
									fmax = fy
									ymax = yw
									no   = 1
								endif

								if(fy<0.9*fmax)goto 50

							enddo

50						    mb_FY(nb)%a4d(i,j,kw,1) = ymax
							mb_FY(nb)%a4d(i,j,kw,2) = fmax
							mb_FY(nb)%a4d(i,j,kw,3) = taow

						else

							fmax = -1.0e10

							do k = s_st(3),1,-1
							    no  = 0

						        vsl = visl(i,j,k)

						        ro  = r(i,j,k)

!		       	                call get_center_cor(nb,nmg,i,j,k,xc,yc,zc)

								xc=x(i,j,k)
							    yc=y(i,j,k)
							    zc=z(i,j,k)

								yw    = sqrt((xc-xwc)**2.0+(yc-ywc)**2.0+(zc-zwc)**2.0)

								yplus = sqrt(ro*taow*reynolds)*yw/vsl

			    	            omg   = sqrt( (dudxyz(i,j,k,2)-dvdxyz(i,j,k,1))**2 + &
				                              (dvdxyz(i,j,k,3)-dwdxyz(i,j,k,2))**2 + &
											  (dwdxyz(i,j,k,1)-dudxyz(i,j,k,3))**2 )

								fy    = yw*omg*(1.0-exp(-yplus/26.0))

								if(fy>fmax)then
									fmax = fy
									ymax = yw
									no   = 1
								endif

								if(fy<0.9*fmax)goto 60

							enddo

60						    mb_FY(nb)%a4d(i,j,kw,1) = ymax
							mb_FY(nb)%a4d(i,j,kw,2) = fmax
							mb_FY(nb)%a4d(i,j,kw,3) = taow

						endif
					enddo
				    enddo
			    endif       
	        endif

	    enddo
	    call deallocate_variables_temp_turbulence
	enddo

	return
end subroutine cal_wall_info

!-----------------------------------------------------------------------------!

subroutine dfidxyz_center_1(nsub,nplus,ni,nj,nk,fi,dfidxyz,method)
	use global_variables,only:kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol
	implicit none

	integer :: nsub,nplus,ni,nj,nk,method,i,j,k,n
	real    :: ddx,ddy,ddz,vol2
	real    :: dkc,det,dct
	real    :: fi(nsub:ni+nplus,nsub:nj+nplus,nsub:nk+nplus)
	real    :: dfidxyz(ni,nj,nk,3)
	real,pointer,dimension(:,:,:,:) :: work


	do n=1,method  ! finite difference method
	do k=1,nk
	do j=1,nj
	do i=1,ni

			dkc = fi(i+1,j,k)-fi(i-1,j,k)
			det = fi(i,j+1,k)-fi(i,j-1,k)
			dct = fi(i,j,k+1)-fi(i,j,k-1)

			ddx = dkc*kcx(i,j,k)+det*etx(i,j,k)+dct*ctx(i,j,k)
			ddy = dkc*kcy(i,j,k)+det*ety(i,j,k)+dct*cty(i,j,k)
			ddz = dkc*kcz(i,j,k)+det*etz(i,j,k)+dct*ctz(i,j,k)

			vol2 =0.5/vol(i,j,k)

			dfidxyz(i,j,k,1) = ddx*vol2
			dfidxyz(i,j,k,2) = ddy*vol2
			dfidxyz(i,j,k,3) = ddz*vol2

	enddo
	enddo
	enddo
	enddo


	return
end subroutine dfidxyz_center_1


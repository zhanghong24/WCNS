
!_____________________________________________________________________!
subroutine post_processing(nout,nbg,ngzone,ngs,ngstep)
	use global_const, only : nwerror,ndisk,nforce,NLAMTUR,nbgmax,ntemp_W,nplot,nchem
    use global_variables,only : ni,nj,nk,x,y,z,fs,rref,r,nust,ndualtst
#ifdef PARALLEL
    use mod_parallels,only : myid,master,write_flow_file_parallel,write_flow_dual_parallel, &
                             write_turbulence_parallel
#endif		
    implicit none
	integer :: nout,nbg,ngzone,ngs,ngstep

!!    if ( nout/nwerror*nwerror == nout )then
!!       call residual_curve(nbg,ngzone)
!!    endif
	
    if ( nout/ndisk*ndisk == nout .or. ngs== ngstep )then
#ifdef PARALLEL
      if (myid == master) then
#endif
      write(*,*)'...... writing data file now......'
#ifdef PARALLEL
      end if
#endif

#ifdef PARALLEL
      call write_flow_file_parallel
      
      if (ndualtst > 0) call write_flow_dual_parallel
      
      if (nlamtur >= 0) call write_turbulence_parallel
#else
      !!call tecflow(nplot)
      call savefile(nout)
      
      if(NLAMTUR .GT. 0) call write_turbulence    !д����������
#endif
      
      if (nplot == 0) then
         call tecout_bin_surface
         call tecout_bin_wall
      else if (nplot == 1) then
         call tecout_bin_wall
      else if (nplot == 2) then
         call tecout_bin_volume
      else if (nplot == 3) then
         call tecout_bin_surface
         call tecout_bin_wall
      else if (nplot == 4) then
         call tecout_bin_2d
         call tecout_bin_wall
      else if (nplot == 5) then
         call tecout_bin_2d_serial   !�����ά��������
         call tecout_bin_wall
	  else
	   write(*,*)'û�������������ѡ��'
      end if

#ifdef PARALLEL
      if (myid == master) then
#endif

       write(*,*)'...... writing finshed      ......'
       
#ifdef PARALLEL
      end if
#endif


!!			 if(nplot ==2 )then ! added by clz 2006.3.20
!!				write(*,*)'  ��ѡ����ȫ�������������ģʽ����ռ�ö����CPUʱ�䣡'
!!				write(*,*)'  ������ѡ��nplot��0��ʹ��Ĭ���������ģʽ'
!!				write(*,*)'  ������������ģʽ��'
!!				read (*,*)nplot
!!			 endif

    endif

    if ( nout/nforce*nforce == nout )then
    
#ifdef PARALLEL
      if (myid == master) then
#endif

      write(*,*)'...... compute aerodynamics......'
      
#ifdef PARALLEL
      end if
#endif

      call force
    endif
			
    return
end subroutine post_processing
!_____________________________________________________________________!
subroutine tecout_bin_surface
    use global_variables
    use mod_tecios
#ifdef PARALLEL
    use mod_parallels
#endif		
    implicit none
    integer :: nb,nbt,nr,nrmax,bctype
    integer :: ibeg(3),iend(3),inrout,idir,pid
    integer :: i,j,k,it,jt,kt,iprn(3),m,ndata,ierr
    integer :: idim,jdim,kdim
#ifdef PARALLEL
    integer :: packsize
    integer :: status(MPI_STATUS_SIZE)
#endif
    logical :: ldraw
    integer :: ntecvars
    character(len=256) :: tecvarnames
    character(len=120) :: zonename
    character(len=32 ) :: str1,str2,str3
    real,pointer :: qpv_tmp(:,:,:,:)
    real :: c2,u2
    
    ntecvars = 10
    tecvarnames = "X,Y,Z,R,U,V,W,P,T,M"
    
    if (nlamtur >= 0) then
       ntecvars = ntecvars + nlamtur + 1
       select case(nlamtur)
       case(0)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,muT"
       case(1)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,muT"
       case(2)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,vT2,muT"
       end select
    end if
    
#ifdef PARALLEL
    if (myid == master) then
#endif
       open(101,file=tecname,form='unformatted', access='stream',status='unknown')
    
       call tecio_ini(101,"WCNS_SOLVER",ntecvars,tecvarnames)
       
       do nb=1,nblocks
          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             bctype = mb_bc(nb)%bc(nr)%bctype
             ldraw  = ( bctype/10 /= 7 )  .and. ( bctype /= 4 ) .and. ( bctype /= 5 ) .and. &
                      ( bctype    > 0  )
          
             if ( ldraw ) then
                ibeg = mb_bc(nb)%bc(nr)%s_st
                iend = mb_bc(nb)%bc(nr)%s_ed
                idim = iend(1)-ibeg(1)+1
                jdim = iend(2)-ibeg(2)+1
                kdim = iend(3)-ibeg(3)+1
                write(str1,'(i6)') nb
                write(str2,'(i6)') nr
                write(str3,'(i6)') bctype
                
                zonename = 'BLK'//trim(adjustl(str1))// &
                           'BC'//trim(adjustl(str2))// &
                           'T'//trim(adjustl(str3))

                call tecio_zone(101,trim(zonename),0,idim,jdim,kdim)
                
             end if
          end do
       end do
       
       call tecio_eohmark(101)
       
       if (prec == single_prec) then
          ndata = 1
       else
          ndata = 2
       end if
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
       nrmax = mb_bc(nb)%nregions
       do nr = 1,nrmax
          bctype = mb_bc(nb)%bc(nr)%bctype
          ldraw  = ( bctype/10 /= 7 )  .and. ( bctype /= 4 ) .and. ( bctype /= 5 ) .and. &
                   ( bctype    > 0  )
          
          if ( ldraw ) then
             ibeg = mb_bc(nb)%bc(nr)%s_st
             iend = mb_bc(nb)%bc(nr)%s_ed
             idir = mb_bc(nb)%bc(nr)%s_nd
             inrout = mb_bc(nb)%bc(nr)%s_lr
             iprn(:) = 0
             iprn(idir) = inrout
             
#ifdef PARALLEL
             pid = mb_pids(nb)
             if (master == pid-1) then
                if (myid == master) then
#endif                
                   allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                                    ibeg(3):iend(3),ntecvars),stat=ierr)
              
                   call tecio_data(101,ndata)
                   
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1) 
                      qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                      c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                      u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                      qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                      
                      if (nlamtur >= 0) then
                         do m=1,nlamtur
                            qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                         end do
                         qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                      end if
                   end do
                   end do
                   end do
                   
                   select case(bctype)
                   case(2,20,21)
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1)
                         it = i - iprn(1)
                         jt = j - iprn(2)
                         kt = k - iprn(3)
                         qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(it,jt,kt)
                         qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(it,jt,kt)
                         qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(it,jt,kt)
                      end do
                      end do
                      end do
                   end select
              
                   write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                             
                   deallocate(qpv_tmp,stat=ierr)
#ifdef PARALLEL                   
                end if
             else
                if (myid == master) then
                   allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                                    ibeg(3):iend(3),ntecvars),stat=ierr)
                                    
                   packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                         
                   call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                                 pid-1,nb,MPI_COMM_WORLD,status,ierr)
                                 
                   call tecio_data(101,ndata)
                                 
                   write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                   
                   deallocate(qpv_tmp,stat=ierr)
                endif
                
                if (myid == pid-1) then
                   allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                                    ibeg(3):iend(3),ntecvars),stat=ierr)
                   
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1) 
                      qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                      c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                      u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                      qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                      
                      if (nlamtur >= 0) then
                         do m=1,nlamtur
                            qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                         end do
                         qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                      end if
                   end do
                   end do
                   end do
                   
                   select case(bctype)
                   case(2,20,21)
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1)
                         it = i - iprn(1)
                         jt = j - iprn(2)
                         kt = k - iprn(3)
                         qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(it,jt,kt)
                         qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(it,jt,kt)
                         qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(it,jt,kt)
                      end do
                      end do
                      end do
                   end select
                  
                   packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                   
                   call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                                 master,nb,MPI_COMM_WORLD,ierr)
                   
                   deallocate(qpv_tmp,stat=ierr)
                end if                
             end if
#endif             
          end if
       end do
    end do
    
#ifdef PARALLEL                   
    if (myid == master) then
#endif
       close(101)
#ifdef PARALLEL                   
    end if
#endif
        
end subroutine tecout_bin_surface
!_____________________________________________________________________!
subroutine tecout_bin_wall
    use global_variables
    use stress_module
    use geometry_module
    use duvwt_module,only : tkc,tet,tct
    use mod_tecios
#ifdef PARALLEL
    use mod_parallels
#endif		
    implicit none
    integer :: nb,nbt,nr,nrmax,bctype
    integer :: ibeg(3),iend(3),inrout,idir,pid
    integer :: i,j,k,it,jt,kt,iprn(3),m,ndata,ierr
    integer :: idim,jdim,kdim
#ifdef PARALLEL
    integer :: packsize
    integer :: status(MPI_STATUS_SIZE)
#endif
    logical :: ldraw
    integer :: ntecvars
    character(len=256) :: tecvarnames
    character(len=120) :: zonename,tecwallname
    character(len=32 ) :: str1,str2,str3
    real,pointer :: qpv_tmp(:,:,:,:)
    real :: pn,nx,ny,nz,fnx,fny,fnz
    real :: dtdx,dtdy,dtdz,cp,kcp,qw
    
    ntecvars = 4
    tecvarnames = "X,Y,Z,cp"
    if ( nvis == 1 ) then
       ntecvars = 9
       tecvarnames = "X,Y,Z,cp,cfx,cfy,cfz,cf,qw"
    end if
    
#ifdef PARALLEL
    if (myid == master) then
#endif
       m = index(tecname,'.',back=.TRUE.)
       if (m == 0) then
          tecwallname = trim(adjustl(tecname)) // "_wall.plt"
       else
          tecwallname = trim(adjustl(tecname(1:m-1))) // "_wall.plt"
       end if
       open(101,file=tecwallname,form='unformatted', access='stream',status='unknown')
    
       call tecio_ini(101,"WCNS_SOLVER",ntecvars,tecvarnames)
       
       do nb=1,nblocks
          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             bctype = mb_bc(nb)%bc(nr)%bctype
             select case(bctype)
             case(2,20,21)
                ldraw = .true.
             case default
                ldraw = .false.
             end select
          
             if ( ldraw ) then
                ibeg = mb_bc(nb)%bc(nr)%s_st
                iend = mb_bc(nb)%bc(nr)%s_ed
                idim = iend(1)-ibeg(1)+1
                jdim = iend(2)-ibeg(2)+1
                kdim = iend(3)-ibeg(3)+1
                write(str1,'(i6)') nb
                write(str2,'(i6)') nr
                write(str3,'(i6)') bctype
                
                zonename = 'BLK'//trim(adjustl(str1))// &
                           'BC'//trim(adjustl(str2))// &
                           'T'//trim(adjustl(str3))

                call tecio_zone(101,trim(zonename),0,idim,jdim,kdim)
                
             end if
          end do
       end do
       
       call tecio_eohmark(101)
       
       if (prec == single_prec) then
          ndata = 1
       else
          ndata = 2
       end if
#ifdef PARALLEL
    end if
#endif

    cp = 1.0/( (gama-1.0) * moo * moo )

    do nb=1,nblocks
       call recast_grid(nb)
       call recast_field(nb)
    
       nrmax = mb_bc(nb)%nregions
       do nr = 1,nrmax
          bctype = mb_bc(nb)%bc(nr)%bctype
          select case(bctype)
          case(2,20,21)
             ldraw = .true.
          case default
             ldraw = .false.
          end select
          
          if ( ldraw ) then
             ibeg = mb_bc(nb)%bc(nr)%s_st
             iend = mb_bc(nb)%bc(nr)%s_ed
             idir = mb_bc(nb)%bc(nr)%s_nd
             inrout = mb_bc(nb)%bc(nr)%s_lr
             iprn(:) = 0
             iprn(idir) = inrout
             
#ifdef PARALLEL
             pid = mb_pids(nb)
             if (master == pid-1) then
                if (myid == master) then
#endif                
                   allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                                    ibeg(3):iend(3),ntecvars),stat=ierr)
              
                   call tecio_data(101,ndata)
                   
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1) 
                      qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                      
                      pn = 2.0 * ( p(i,j,k) - poo )
                      qpv_tmp(i,j,k,4) = pn

                      if ( nvis == 1 ) then
                         call getnxyz_mml(nx,ny,nz,i,j,k,iprn(1),iprn(2),iprn(3))
                         call getgeo(i,j,k)
                         call getuvwtder(i,j,k)
                         call getduvwdxyz
                         call stress( visl(i,j,k) )
                         pn = 2.0*inrout
                         fnx = - pn* ( nx*txx + ny*tyx + nz*tzx ) * re
                         fny = - pn* ( nx*txy + ny*tyy + nz*tzy ) * re
                         fnz = - pn* ( nx*txz + ny*tyz + nz*tzz ) * re
                         
                         dtdx = (kx*tkc + ex*tet + cx*tct)/vjacob
                         dtdy = (ky*tkc + ey*tet + cy*tct)/vjacob
                         dtdz = (kz*tkc + ez*tet + cz*tct)/vjacob
                         kcp  =  visl(i,j,k)*cp/prl
                         qw   = - inrout*kcp* ( nx * dtdx + ny * dtdy + nz * dtdz ) * re
                         
                         qpv_tmp(i,j,k,5) = fnx
                         qpv_tmp(i,j,k,6) = fny
                         qpv_tmp(i,j,k,7) = fnz
                         qpv_tmp(i,j,k,8) = sqrt(fnx*fnx + fny*fny + fnz*fnz)
                         qpv_tmp(i,j,k,9) = qw
                      end if                      
                   end do
                   end do
                   end do

                   write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                             
                   deallocate(qpv_tmp,stat=ierr)
#ifdef PARALLEL                   
                end if
             else
                if (myid == master) then
                   allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                                    ibeg(3):iend(3),ntecvars),stat=ierr)
                                    
                   packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                         
                   call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                                 pid-1,nb,MPI_COMM_WORLD,status,ierr)
                                 
                   call tecio_data(101,ndata)
                                 
                   write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                   
                   deallocate(qpv_tmp,stat=ierr)
                endif
                
                if (myid == pid-1) then
                   allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                                    ibeg(3):iend(3),ntecvars),stat=ierr)
                   
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1) 
                      qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                      qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                      
                      pn = 2.0 * ( p(i,j,k) - poo )
                      qpv_tmp(i,j,k,4) = pn

                      if ( nvis == 1 ) then
                         call getnxyz_mml(nx,ny,nz,i,j,k,iprn(1),iprn(2),iprn(3))
                         call getgeo(i,j,k)
                         call getuvwtder(i,j,k)
                         call getduvwdxyz
                         call stress( visl(i,j,k) )
                         pn = 2.0*inrout
                         fnx = - pn* ( nx*txx + ny*tyx + nz*tzx ) * re
                         fny = - pn* ( nx*txy + ny*tyy + nz*tzy ) * re
                         fnz = - pn* ( nx*txz + ny*tyz + nz*tzz ) * re
                         
                         dtdx = (kx*tkc + ex*tet + cx*tct)/vjacob
                         dtdy = (ky*tkc + ey*tet + cy*tct)/vjacob
                         dtdz = (kz*tkc + ez*tet + cz*tct)/vjacob
                         kcp  =  visl(i,j,k)*cp/prl
                         qw   = - inrout*kcp* ( nx * dtdx + ny * dtdy + nz * dtdz ) * re
                         
                         qpv_tmp(i,j,k,5) = fnx
                         qpv_tmp(i,j,k,6) = fny
                         qpv_tmp(i,j,k,7) = fnz
                         qpv_tmp(i,j,k,8) = sqrt(fnx*fnx + fny*fny + fnz*fnz)
                         qpv_tmp(i,j,k,9) = qw
                      end if                      
                   end do
                   end do
                   end do
                  
                   packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                   
                   call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                                 master,nb,MPI_COMM_WORLD,ierr)
                   
                   deallocate(qpv_tmp,stat=ierr)
                end if                
             end if
#endif             
          end if
       end do
    end do
    
#ifdef PARALLEL                   
    if (myid == master) then
#endif
       close(101)
#ifdef PARALLEL                   
    end if
#endif
        
end subroutine tecout_bin_wall
!_____________________________________________________________________!
subroutine tecout_bin_volume
    use global_variables
    use mod_tecios
#ifdef PARALLEL
    use mod_parallels
#endif		
    implicit none
    integer :: nb,nbt
    integer :: ibeg(3),iend(3),pid
    integer :: i,j,k,m,ndata,ierr
    integer :: idim,jdim,kdim
#ifdef PARALLEL
    integer :: packsize
    integer :: status(MPI_STATUS_SIZE)
#endif
    logical :: ldraw
    integer :: ntecvars
    character(len=256) :: tecvarnames
    character(len=120) :: zonename,tecvolname
    character(len=32 ) :: str1
    real,pointer :: qpv_tmp(:,:,:,:)
    real :: c2,u2
    
    ntecvars = 10
    tecvarnames = "X,Y,Z,R,U,V,W,P,T,M"
    
    if (nlamtur >= 0) then
       ntecvars = ntecvars + nlamtur + 1
       select case(nlamtur)
       case(0)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,muT"
       case(1)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,muT"
       case(2)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,vT2,muT"
       end select
    end if
    
#ifdef PARALLEL
    if (myid == master) then
#endif
       m = index(tecname,'.',back=.TRUE.)
       if (m == 0) then
          tecvolname = trim(adjustl(tecname)) // "_vol.plt"
       else
          tecvolname = trim(adjustl(tecname(1:m-1))) // "_vol.plt"
       end if

       open(101,file=tecvolname,form='unformatted', access='stream',status='unknown')
    
       call tecio_ini(101,"WCNS_SOLVER",ntecvars,tecvarnames)
      
       do nb=1,nblocks
          ibeg = 1
          iend = mb_dim(nb,:)
          idim = iend(1)-ibeg(1)+1
          jdim = iend(2)-ibeg(2)+1
          kdim = iend(3)-ibeg(3)+1
          write(str1,'(i6)') nb
          
          zonename = 'BLK'//trim(adjustl(str1))

          call tecio_zone(101,trim(zonename),0,idim,jdim,kdim)
                
       end do
       
       call tecio_eohmark(101)
       
       if (prec == single_prec) then
          ndata = 1
       else
          ndata = 2
       end if
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
       ibeg = 1
       iend = mb_dim(nb,:)
             
#ifdef PARALLEL
       pid = mb_pids(nb)
       if (master == pid-1) then
          if (myid == master) then
#endif          
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
        
             call tecio_data(101,ndata)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
             
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                       
             deallocate(qpv_tmp,stat=ierr)
#ifdef PARALLEL                   
          end if
       else
          if (myid == master) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
                              
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                   
             call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                           pid-1,nb,MPI_COMM_WORLD,status,ierr)
                           
             call tecio_data(101,ndata)
                           
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
             
             deallocate(qpv_tmp,stat=ierr)
          endif
          
          if (myid == pid-1) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
            
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
             
             call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                           master,nb,MPI_COMM_WORLD,ierr)
             
             deallocate(qpv_tmp,stat=ierr)
          end if                
       end if
#endif             
    end do
    
#ifdef PARALLEL                   
    if (myid == master) then
#endif
       close(101)
#ifdef PARALLEL                   
    end if
#endif
        
end subroutine tecout_bin_volume
!_____________________________________________________________________!
subroutine tecflow( plot_strategy )
    use global_variables,only : nm,mb_varname,method,tecname
    implicit none
    character(len=120) :: title_tecplot
    integer :: nv,mtecmax,plot_strategy

    title_tecplot ='variables="x" "y" "z" "r" "u" "v" "w" "p" "t" "m" "gamma" "vist"'

    mtecmax = 12

!    if ( (nm+1) <= mb_varnumber ) then
!       mtecmax = mtecmax + mb_varnumber  - nm
!    endif

    open(1,file=tecname,status='unknown')
        if ( plot_strategy == 2 ) then
           if ( method == 1 ) then
              call diff_tecplot(1,title_tecplot,mtecmax)
           else
!              call  vol_tecplot(1,title_tecplot,mtecmax)
           endif
        else
           call plot_default(1,title_tecplot,mtecmax)
           if ( plot_strategy == 1 ) then
              call plot_user_def(1,title_tecplot,mtecmax)
           endif
        endif
    close(1)

    return
end subroutine tecflow
!_____________________________________________________________________!
subroutine diff_tecplot(file_id,title_tecplot,mtecmax)
    use global_variables,only : mb_dim,ni,nj,nk,nblocks,nm,mb_varname
    implicit none
    integer :: file_id
    integer :: mtecmax,nb,m,i,j,k,nv
    character(len=120) :: title_tecplot
    real,pointer,dimension( : ) :: tec_var
    real,pointer,dimension( :,:,:,: ) :: block_tecplot

    allocate( tec_var(mtecmax) )
    do nb=1,nblocks
       ni = mb_dim(nb,1)
       nj = mb_dim(nb,2)
       nk = mb_dim(nb,3)

       write(file_id,*)trim(title_tecplot)
!       do nv=nm+1,mb_varnumber
!          write(file_id,*) mb_varname(nv) 
!       enddo
       write(file_id,*)'zone i=',ni,'j=',nj,'k=',nk,' f=block'
       allocate( block_tecplot(mtecmax,ni,nj,nk) )

       call get_block_tecplot(nb,block_tecplot,tec_var,mtecmax)
       write(file_id,*)((((block_tecplot(m,i,j,k),i=1,ni),j=1,nj),k=1,nk),m=1,mtecmax)
       deallocate( block_tecplot )
    enddo

    deallocate( tec_var )

    return
end subroutine diff_tecplot
!_____________________________________________________________________!
subroutine plot_default(file_id,title_tecplot,mtecmax)
    use global_variables, &
	only : mb_dim,ni,nj,nk,nblocks,mb_bc,method,nm,mb_varname

    implicit none
    character(len=120) :: title_tecplot
    logical :: ldraw
    integer :: file_id
    integer :: i,j,k
    integer :: nv
    integer :: nfile_par,num_var,nv_s,nv_t
    integer :: nb,nr,nrmax,mml,mtecmax
    integer :: bctype,m
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    integer :: s0(3),ii1,ii2,jj1,jj2,kk1,kk2
    real :: xc,yc,zc
    real,pointer,dimension( : ) :: tec_var
    real,pointer,dimension( :,:,:,: ) :: block_tecplot
    real,pointer,dimension( :,:,:,: ) :: block_whole

    mml = 1 - method

    allocate( tec_var(mtecmax) )
    do nb=1,nblocks
       ni = mb_dim(nb,1)
       nj = mb_dim(nb,2)
       nk = mb_dim(nb,3)
       allocate( block_whole( mtecmax,ni,nj,nk ) )

       call get_block_tecplot(nb,block_whole,tec_var,mtecmax)

       nrmax = mb_bc(nb)%nregions
       do nr=1,nrmax
          bctype = mb_bc(nb)%bc(nr)%bctype
          ldraw  = ( bctype/10 /= 7 )  .and. ( bctype /= 4 ) .and. ( bctype /= 5 ) .and. &
                   ( bctype    > 0  )
          if ( ldraw ) then
             do m=1,3
                s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
                s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
                if ( s_st(m) == s_ed(m) ) then
                   if ( s_st(m) > 1 ) then
                      s_st(m) = mb_bc(nb)%bc(nr)%s_st(m) + mml
                      s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m) + mml
                   endif
                else
                   s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m) + mml
                endif
             enddo

             write(file_id,*)trim(title_tecplot)
 !            do nv=nm+1,mb_varnumber
 !               write(file_id,*) mb_varname(nv) 
 !            enddo

             write(file_id,*)'zone i=',s_ed(1)-s_st(1)+1,'j=',s_ed(2)-s_st(2)+1,'k=',s_ed(3)-s_st(3)+1,' f=block'
             allocate( block_tecplot(mtecmax,s_st(1):s_ed(1),s_st(2):s_ed(2),s_st(3):s_ed(3)) )

             do k=s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)
                   do i=s_st(1),s_ed(1)
                      do m=1,mtecmax
                         block_tecplot(m,i,j,k) = block_whole(m,i,j,k)  !tec_var(m)
                      enddo
                   enddo
                enddo
             enddo
             write(file_id,*)((((block_tecplot(m,i,j,k),i=s_st(1),s_ed(1)),j=s_st(2),s_ed(2)),k=s_st(3),s_ed(3)),m=1,mtecmax)
             deallocate( block_tecplot )
          endif
       enddo
       deallocate( block_whole )
    enddo

    deallocate( tec_var )

    return
end subroutine plot_default
!_____________________________________________________________________!
subroutine plot_user_def(file_id,title_tecplot,mtecmax)
    use global_variables,only : method,mb_dim,ni,nj,nk,nm,mb_varname
    implicit none
    integer :: file_id,nd
    integer :: num_user_plot,n
    integer :: nb,nr,nrmax,mtecmax,nv
    integer :: i,j,k,m,id,jd,kd,s_lr3d(3)
    integer :: imin,imax,jmin,jmax,kmin,kmax
    character(len=120) :: title_tecplot
    real,pointer,dimension( : ) :: tec_var
    real,pointer,dimension( :,:,:,: ) :: block_tecplot
    real,pointer,dimension( :,:,:,: ) :: block_whole

    allocate( tec_var(mtecmax) )

    open(2,file='output/userdef.dat',status='old')
    read(2,*)num_user_plot
    do n=1,num_user_plot
       read(2,*)imin,imax,jmin,jmax,kmin,kmax,nb
       ni = mb_dim(nb,1)
       nj = mb_dim(nb,2)
       nk = mb_dim(nb,3)
       allocate( block_whole( mtecmax,ni,nj,nk ) )

       call get_block_tecplot(nb,block_whole,tec_var,mtecmax)
       !����������
       !call odd_axi_tec(block_whole,nb,mtecmax)
       allocate( block_tecplot(mtecmax,imin:imax,jmin:jmax,kmin:kmax) )
       write(file_id,*)trim(title_tecplot)
!       do nv=nm+1,mb_varnumber
!          write(file_id,*) mb_varname(nv) 
!       enddo
       write(file_id,*)'zone i=',imax-imin+1,'j=',jmax-jmin+1,'k=',kmax-kmin+1,' f=block'
       do k=kmin,kmax
          do j=jmin,jmax
             do i=imin,imax
                do m=1,mtecmax
                   block_tecplot(m,i,j,k) = block_whole(m,i,j,k)
                enddo
             enddo
          enddo
       enddo
       write(file_id,*)((((block_tecplot(m,i,j,k),i=imin,imax),j=jmin,jmax),k=kmin,kmax),m=1,mtecmax)
       deallocate( block_tecplot )
       deallocate( block_whole   )
    enddo
    close(2)

    deallocate( tec_var )

    return
end subroutine plot_user_def
!_____________________________________________________________________!
subroutine get_block_tecplot(nb,block_tecplot,tec_var,mtecmax)
    use global_variables
    implicit none
    integer :: i,j,k,dii,djj,dkk
    integer :: nv
    integer :: nfile_par,num_var,nv_s,nv_t
    integer :: nb,nr,nrmax,mml,mtecmax
    integer :: bctype,m,is,js,ks,s_nd
    integer :: s_st(3),s_ed(3),nd,s_lr3d(3)
    integer :: naxir, bcaxi,nd3,is1,i1
    real :: xc,yc,zc,cn,cm
    real :: tec_var(mtecmax)
    real :: block_tecplot(mtecmax,ni,nj,nk)

    ni = mb_dim(nb,1)
    nj = mb_dim(nb,2)
    nk = mb_dim(nb,3)
    mml = 1 - method
    !�����ڵ�
    do m=1,3
       s_lr3d(m) = 0
    enddo
    do k=2,nk-1
       do j=2,nj-1
          do i=2,ni-1
             call get_tec_var(nb,i,j,k,method,s_lr3d,tec_var,mtecmax)
             do m=1,mtecmax
                block_tecplot(m,i,j,k) = tec_var(m)
             enddo
          enddo
       enddo
    enddo

    !�����߽�
    do m=1,3
       s_lr3d(m) = 0
    enddo

    do k=1,nk,nk-1
       if ( k == 1  ) s_lr3d(3) = -1
       if ( k == nk ) s_lr3d(3) =  1
       do j=1,nj
          do i=1,ni
             call get_tec_var(nb,i,j,k,method,s_lr3d,tec_var,mtecmax)
             do m=1,mtecmax
                block_tecplot(m,i,j,k) = tec_var(m)
             enddo
          enddo
       enddo
    enddo
    do m=1,3
       s_lr3d(m) = 0
    enddo
    do k=2,nk-1
       do j=1,nj,nj-1
          if ( j == 1  ) s_lr3d(2) = -1
          if ( j == nj ) then
		     s_lr3d(2) =  1
		  endif
          do i=1,ni
             call get_tec_var(nb,i,j,k,method,s_lr3d,tec_var,mtecmax)
             do m=1,mtecmax
                block_tecplot(m,i,j,k) = tec_var(m)
             enddo
          enddo
       enddo
    enddo

    do m=1,3
       s_lr3d(m) = 0
    enddo
    do k=2,nk-1
       do j=2,nj-1
          do i=1,ni,ni-1
             if ( i == 1  ) s_lr3d(1) = -1
             if ( i == ni ) s_lr3d(1) =  1
             call get_tec_var(nb,i,j,k,method,s_lr3d,tec_var,mtecmax)
             do m=1,mtecmax
                block_tecplot(m,i,j,k) = tec_var(m)
             enddo
          enddo
       enddo
    enddo


	!�����ٶ�����
	nrmax = mb_bc(nb)%nregions
    do nr=1,nrmax
       dii = 0
       djj = 0
       dkk = 0
       bctype = mb_bc(nb)%bc(nr)%bctype
       if ( ( bctype/10 == 2 ) .or. ( bctype == 2 ) ) then
           do m=1,3
              s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
              s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
              if ( s_st(m) == s_ed(m) ) then
                 if ( m == 1 .and. method == 1 ) dii =  1
                 if ( m == 1 .and. s_st(m) > 1 ) dii = -1
                 if ( m == 2 .and. method == 1 ) djj =  1
                 if ( m == 2 .and. s_st(m) > 1 ) djj = -1
                 if ( m == 3 .and. method == 1 ) dkk =  1
                 if ( m == 3 .and. s_st(m) > 1 ) dkk = -1
              endif
           enddo
		   if( abs(dii+djj+dkk) == 1 ) then
             do k = s_st(3)-dkk*mml , s_ed(3)+mml
                do j = s_st(2)-djj*mml , s_ed(2)+mml
                   do i = s_st(1)-dii*mml , s_ed(1)+mml
                      do m = 5 , 7
                         block_tecplot(m,i,j,k) = block_tecplot(m,i+dii,j+djj,k+dkk)
                      enddo
                   enddo
                enddo
             enddo
		   endif
        endif
     enddo

    return
end subroutine get_block_tecplot
!_____________________________________________________________________!
subroutine get_tec_var(nb,i,j,k,method,s_lr3d,tec_var,mtecmax)
    use global_variables, &
        only : nm,nchem,mb_varnumber,mb_par_content,mb_control,gama,small, &
           mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,mb_fs,qb_seq,mb_vist, &
		   mb_distance,mb_qke,mb_rhs1,mb_vol,mb_etx,mb_ety,mb_etz,mb_ctx,mb_cty,mb_ctz,mb_kcx,mb_kcy,mb_kcz
    implicit none
    integer :: is,js,ks,method,s_lr3d(3)
    integer :: file_id,nd,mtecmax
    integer :: nb,i,j,k,m,ii,jj,kk,id,jd,kd
    integer :: nfile_par,num_var,nv_s,nv_t
    real :: xc,yc,zc,coe,gama1
    real :: rm,um,vm,wm,pm,mm,cm,tm
    real :: tec_var(mtecmax),mag1,mag2,mag3
    integer :: mmm

    if ( method == 1 ) then
       id = 0
       jd = 0
       kd = 0
       is = i - s_lr3d(1)
       js = j - s_lr3d(1)
       ks = k - s_lr3d(1)
    else
       id = 1 - abs( s_lr3d(1) )
       jd = 1 - abs( s_lr3d(2) )
       kd = 1 - abs( s_lr3d(3) )
       if ( s_lr3d(1)==1 .or. s_lr3d(2)==1 .or. s_lr3d(3)==1 ) then
          is = i - s_lr3d(1)
          js = j - s_lr3d(2)
          ks = k - s_lr3d(3)
       else
          is = i
          js = j
          ks = k
       endif
    endif
    coe = 1.0/( ( id + 1) * (1 + jd )* ( 1 + kd ) )
    rm  = 0.0
    um  = 0.0
    vm  = 0.0
    wm  = 0.0
    pm  = 0.0
    tm  = 0.0
	gama1 = 0.0
!    do m=nm+1,mb_varnumber
!       tec_var(12 + m - nm ) = 0.0
!    enddo
!    nfile_par = mb_par_content(nb)
!    num_var   = mb_control( nfile_par )%num_variable
!    nchem     = mb_control( nfile_par )%nchem
    do ii=-id,0
       do jj=-jd,0
          do kk=-kd,0
             rm  = rm + coe * mb_r(nb)%a3d(is+ii,js+jj,ks+kk)
             um  = um + coe * mb_u(nb)%a3d(is+ii,js+jj,ks+kk)
             vm  = vm + coe * mb_v(nb)%a3d(is+ii,js+jj,ks+kk)
             wm  = wm + coe * mb_w(nb)%a3d(is+ii,js+jj,ks+kk)
             pm  = pm + coe * mb_p(nb)%a3d(is+ii,js+jj,ks+kk)
             tm  = tm + coe * mb_t(nb)%a3d(is+ii,js+jj,ks+kk)

!             gama = mb_control( nfile_par )%gama
!             do m=nm+1,mb_varnumber
!                tec_var(12 + m - nm ) = coe*qb_seq(nfile_par,m) &
!                                        + tec_var(12 + m - nm )
!             enddo
          enddo
       enddo
    enddo

!_____

             rm  =  mb_r(nb)%a3d(i,j,k)
             um  =  mb_u(nb)%a3d(i,j,k)
             vm  =  mb_v(nb)%a3d(i,j,k)
             wm  =  mb_w(nb)%a3d(i,j,k)
             pm  =  mb_p(nb)%a3d(i,j,k)
             tm  =  mb_t(nb)%a3d(i,j,k)
!_____

    cm = sqrt(gama*pm/rm)
    tec_var( 11 ) = gama
    mm = sqrt(um*um + vm*vm + wm*wm)/( cm + small )
    call xyzcoor(nb,i,j,k,xc,yc,zc)
	mag1 = sqrt(mb_kcx(nb)%a3d(i,j,k)**2.0 + mb_kcy(nb)%a3d(i,j,k)**2.0 + mb_kcz(nb)%a3d(i,j,k)**2.0)
	mag2 = sqrt(mb_etx(nb)%a3d(i,j,k)**2.0 + mb_ety(nb)%a3d(i,j,k)**2.0 + mb_etz(nb)%a3d(i,j,k)**2.0)
	mag3 = sqrt(mb_ctx(nb)%a3d(i,j,k)**2.0 + mb_cty(nb)%a3d(i,j,k)**2.0 + mb_ctz(nb)%a3d(i,j,k)**2.0)
    tec_var( 1  ) = xc
    tec_var( 2  ) = yc
    tec_var( 3  ) = zc
    tec_var( 4  ) = rm
    tec_var( 5  ) = um
    tec_var( 6  ) = vm
    tec_var( 7  ) = wm
    tec_var( 8  ) = pm
    tec_var( 9  ) = tm
    tec_var( 10 ) = mm
    tec_var( 11 ) = gama
    tec_var( 12 ) = mb_vist(nb)%a3d(i,j,k)


    return
end subroutine get_tec_var
!_____________________________________________________________________!
subroutine residual_curve(nbg,ngzone)
    use global_variables
    use resdual_module
#ifdef PARALLEL
    use mod_parallels,only : master,myid,numprocs,pnblocks,pnbindexs, &
                             mpi_inprec,mpi_reprec,MPI_COMM_WORLD,MPI_SUM
#endif
    implicit none
    integer :: nb,m,np,np_block,m1,i,j,k
    integer :: nbg,ngz,ngzone
    real :: dtime1,dq0,rdqn,rdqm,res_dts(nl),res_dts0(nl)
    real :: dres(nblocks),dresm,resall_dts,tol0
#ifdef PARALLEL
    integer :: pnb,np_tmp,id_max,ierr
    real :: dres_add_tmp,resmax_tmp1(numprocs),resmax_tmp2(numprocs)
#endif
    
    m1 = 1 + method
    do m=1,nblocks
       res(m) = 0.0
       resmax(m) = 0.0
       nres(m,1) = 2
       nres(m,2) = 2
       nres(m,3) = 2
       nres(m,4) = 1
    enddo
    np = 0
    dres_add = 0.0
    resmax_total = 0.0
    
    do m=1,nl
      res_dts0(m) = 0.0
      res_dts(m)  = 0.0
    end do
    
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
!    do ngz=1,ngzone
       do nb =1,nblocks   !mb_group( nbg )%s_nb(ngz), mb_group( nbg )%e_nb(ngz)
#endif
!          call recast_method(nb)
          call recast_grid(nb)
          call recast_field(nb)
          do k=m1,nk-1   !2,nk-1
             do j=m1,nj-1   !2,nj-1
                do i=m1,ni-1   !2,ni-1
                   dtime1 = 1.0 /( dtdt(i,j,k)*vol(i,j,k) )
                   do m=1,nl    !nm
                      dq0 = abs(dq(m,i,j,k))
                      if (ndualtst > 0) then
                         rdqn = dq0
                         res_dts(m) = res_dts(m) + rdqn*rdqn

                         rdqm = abs(q(m,i,j,k) - mb_qnc(nb)%a4d(m,i,j,k))
                         res_dts0(m) = res_dts0(m) + rdqm*rdqm
                         dresm  = rdqm
                      else
                         dresm  = dq0*dtime1
                      end if
                      
                      if ( dresm > resmax(nb) ) then
                         resmax(nb) = dresm
                         nres(nb,1) = i
                         nres(nb,2) = j
                         nres(nb,3) = k
                         nres(nb,4) = m
                      endif
                      res(nb) = res(nb) + dresm*dresm
                   enddo
                enddo
             enddo
          enddo
          dres_add = dres_add + res(nb)
!         np_block = (ni-2)*(nj-2)*(nk-2)*nl
          np_block = (ni-m1)*(nj-m1)*(nk-m1)*nl ! modified by clz
          res(nb)  = res(nb)/np_block
          np       = np + np_block
          if(resmax_total < resmax(nb) ) then
             resmax_total = resmax(nb)
             do m=1,4
                n_error(m) = nres(nb,m)
             enddo
             n_error(5) = nb
          endif
       enddo
!    enddo

#ifdef PARALLEL
    np_tmp = np
    call MPI_REDUCE(np_tmp,np,1,mpi_inprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(np,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
    
    dres_add_tmp = dres_add
    call MPI_REDUCE(dres_add_tmp,dres_add,1,mpi_reprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dres_add,1,mpi_reprec,master,MPI_COMM_WORLD,ierr)
#endif
    dres_add = sqrt(dres_add/np)
    
    if (ndualtst > 0) then
#ifdef PARALLEL
       call MPI_REDUCE(res_dts(1),resall_dts,1,mpi_reprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(resall_dts,1,mpi_reprec,master,MPI_COMM_WORLD,ierr)
       call MPI_REDUCE(res_dts0(1),resall_dts0,1,mpi_reprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(resall_dts0,1,mpi_reprec,master,MPI_COMM_WORLD,ierr)
#else
       resall_dts = res_dts(1)
       resall_dts0 = res_dts0(1)
#endif
       
       tol0 = resall_dts/(resall_dts0 + small)
       tol0 = sqrt(tol0)

#ifdef PARALLEL
       if (myid == master) then
#endif
          write(*,'(5x,a,i5,3(1x,e12.5))') "DTS:",nsubstep,tol0,tolsub
#ifdef PARALLEL
       end if
#endif
       
       if ((tol0 < tolsub .and. nsubstep > 1) .or. nsubstep == nsubstmx) then
          nstopsub = 1
       else
          nstopsub = 0
       end if
    end if        
    
    
!    call recast_method( n_error(5) )
    call recast_field( n_error(5) )
    call recast_grid( n_error(5) )
    i= n_error(1)
    j= n_error(2)
    k= n_error(3)

!!    dtime1 = dtdt(i,j,k)*vol(i,j,k)
    
#ifdef PARALLEL
    resmax_tmp1(:) = 0.0
    resmax_tmp1(myid + 1) = resmax_total
    call MPI_REDUCE(resmax_tmp1,resmax_tmp2,numprocs,mpi_reprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(resmax_tmp2,numprocs,mpi_reprec,master,MPI_COMM_WORLD,ierr)
    id_max = maxloc(resmax_tmp2,dim=1)
    resmax_total = resmax_tmp2(id_max)
    call MPI_BCAST(n_error,5,mpi_inprec,id_max-1,MPI_COMM_WORLD,ierr)
#endif

    !!call saveres2file(nbg,ngzone)
    
    return
end subroutine residual_curve
!_____________________________________________________________________!
subroutine saveres2file(nbg,ngzone)
    use global_variables
    use resdual_module
#ifdef PARALLEL
    use mod_parallels,only : master,myid
#endif
    implicit none
    integer :: nbg,ngzone,m
    
#ifdef PARALLEL
    if (myid == master) then
#endif
    m = 10*nwerror
    if ( (nout-nwerror)/m*m == (nout-nwerror) ) then
       write(*,*)
       write(*,*)'   iter     CPU        cfl        ave         max       nb   ni   nj   nk   nv'
    endif

    open(1,file=errhis,access='append',status='unknown')
!    write(*,10)nout,CPU_TIME_USED,cfl,dres_add,resmax_total,n_error(5),(n_error(m),m=1,4)
    write(*,10)nout,CPU_TIME_USED/3600.,cfl,dres_add,resmax_total,n_error(5),(n_error(m),m=1,4)
    write(1,10)nout,CPU_TIME_USED,cfl,dres_add,resmax_total,n_error(5),(n_error(m),m=1,4)
10  format(i7,2x,e10.3,f9.3,2x,e11.4,1x,e11.4,1x,5(1x,i4))
    close(1)
    
#ifdef PARALLEL
    end if
#endif

end subroutine saveres2file
!_____________________________________________________________________!
subroutine savefile( it )
    use global_variables,&
    only : wholetime,nblocks,nchem,ns,ni,nj,nk,nout,flowname,x,y,z,r,u,v,w,p,t,fs,&
           roo,uoo,voo,woo,poo,moo,nlamtur,rref,ndualtst,ttdts,mb_qmc,nl
    implicit none
    integer :: nb,i,j,k,m,it

    open(1,file=flowname,status='unknown',form='unformatted')
    write(1)nout,wholetime,roo,uoo,voo,woo,poo,moo
    do nb=1,nblocks
!       call recast_method(nb)
       call recast_grid(nb)
       call recast_field(nb)

       write(1)(((r(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),t(i,j,k), &
                                          i=1,ni),j=1,nj),k=1,nk)
    enddo
	close(1)

    if (ndualtst > 0) then
    
    open(1,file=trim(flowname)//".dts",status='unknown',form='unformatted')
    write(1)ttdts
    do nb=1,nblocks
       call recast_grid(nb)
       call recast_field(nb)
       write(1)((((mb_qmc(nb)%a4d(m,i,j,k),m=1,nl),i=1,ni),j=1,nj),k=1,nk)
    enddo
    close(1)
    
    end if

    return
end subroutine savefile
!_____________________________________________________________________!
subroutine force
   !���ӳ������������������������ϵ��
    use global_variables,&
    only : nwholefield,method,nout,attack,lfref,sref,CPU_TIME_USED &
		     &,  forcename,cfx,cfy,cfz,cmx,cmy,cmz,cl,cd,xcp,nvis &
		     &,  phytime,phydtime,small,sideslip
#ifdef PARALLEL
    use mod_parallels
#endif
    implicit none 
    real :: scf,vcm
    real :: sina,cosa,sinb,cosb
    real :: forces(6),pforces(6)
    integer :: ierr
    
    !��������ϵ������������ϵ��Ϊ0
    cfx = 0.0
    cfy = 0.0
    cfz = 0.0

    cmx = 0.0
    cmy = 0.0
    cmz = 0.0
    if ( method == 1 ) call force_dif

#ifdef PARALLEL
    forces(1) = cfx
    forces(2) = cfy
    forces(3) = cfz
    forces(4) = cmx
    forces(5) = cmy
    forces(6) = cmz
    call MPI_REDUCE(forces,pforces,6,mpi_reprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(pforces,6,mpi_reprec,master,MPI_COMM_WORLD,ierr)
    cfx = pforces(1)
    cfy = pforces(2)
    cfz = pforces(3)
    cmx = pforces(4)
    cmy = pforces(5)
    cmz = pforces(6)
#endif

    scf = 1.0 / sref
    vcm = 1.0 / ( sref * lfref )

    cfx = ( 2 - nwholefield ) * cfx * scf
    cfy = ( 2 - nwholefield ) * cfy * scf
    cfz = ( 2 - nwholefield ) * cfz * scf
    cmx = ( 2 - nwholefield ) * cmx * vcm
    cmy = ( 2 - nwholefield ) * cmy * vcm
    cmz = ( 2 - nwholefield ) * cmz * vcm

    sina = sin(attack)
    cosa = cos(attack)
	sinb = sin(sideslip)
	cosb = cos(sideslip)

    !cl�� cd�ֱ�Ϊ����ϵ��������ϵ��
!    cl = cfy * cosa - cfx * sina
!    cd = cfy * sina + cfx * cosa
    cl = cfy * cosa - cfx * sina
    cd = cosb*(cfy * sina + cfx * cosa)+sinb*cfz
    xcp = sign(1.0,cfy) * cmz/( abs(cfy) + small )

!!    if ( nvis >= 1 ) then
!!       if ( method == 1 ) call heatflux_dif
!!    endif
    
#ifdef PARALLEL
    if (myid == master) then
#endif

    open(1,file=forcename,access='append',status='unknown')
    write(1,10)nout,CPU_TIME_USED,cfx,cfy,cfz,cmx,cmy,cmz,cd,cl,xcp
    close(1)
10  format(1x,i6,2x,e10.4,9(f12.5,1x))

    write(*,20)cfx,cfy,cfz,cmx,cmy,cmz,xcp
20  format(1x,7(f10.6,1x),f10.7)

     write(*,21)"time=",phytime,"   dtime=",phydtime
21  format(1x,(a5,e12.5,a9,e12.5))

#ifdef PARALLEL
    end if
#endif

   return
end subroutine force
!_____________________________________________________________________!
subroutine force_dif
    !���ӳ������������������������ϵ��
    use global_variables
    use stress_module
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none 
    integer :: nb,nr,nrmax,pnb
    integer :: bctype,m,is,js,ks,id,jd,kd,s_nd,s_fix,s_lr
    integer :: ii1,jj1,kk1,ii2,jj2,kk2,ii3,jj3,kk3,ii4,jj4,kk4
    integer :: s_st(3),s_ed(3),s_lr3d(3)
    real :: nx,ny,nz
    real :: pn
    real :: dfx,dfy,dfz,dmx,dmy,dmz
    real :: px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,px4,py4,pz4
    real :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    real,pointer,dimension( :,:,: ) :: pnx,pny,pnz
    
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
!       call recast_method(nb)
       call recast_grid(nb)
       call recast_field(nb)
       nrmax = mb_bc(nb)%nregions                     !���鹲��nrmax���߽���Ҫ����
       do nr=1,nrmax
          bctype = mb_bc(nb)%bc(nr)%bctype
          if ( bctype == 2 .or. bctype/10 == 2 ) then  !�̱ڱ߽磬��Ӧ������
             do m=1,3
                s_st(m)   = mb_bc(nb)%bc(nr)%s_st(m)  !��ʼ������(�ɶ�����)
                s_ed(m)   = mb_bc(nb)%bc(nr)%s_ed(m)  !��ֹ������(�ɶ�����)
                s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)
             enddo

             id = 1 - abs(s_lr3d(1))
             jd = 1 - abs(s_lr3d(2))
             kd = 1 - abs(s_lr3d(3))

             s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
             s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
             s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
             allocate( pnx(s_st(1):s_ed(1),s_st(2):s_ed(2),s_st(3):s_ed(3)) )
             allocate( pny(s_st(1):s_ed(1),s_st(2):s_ed(2),s_st(3):s_ed(3)) )
             allocate( pnz(s_st(1):s_ed(1),s_st(2):s_ed(2),s_st(3):s_ed(3)) )

             !��ȡ����Ӧ��
             do is=s_st(1),s_ed(1)
                do js=s_st(2),s_ed(2)
                   do ks=s_st(3),s_ed(3)
                      call getnxyz_mml(nx,ny,nz,is,js,ks,s_lr3d(1),s_lr3d(2),s_lr3d(3))
!                     pn = p(s1(1),s1(2),s1(3))
                      pn = 2.0 * ( p(is,js,ks) - poo )

                      pnx(is,js,ks) = s_lr * nx * pn
                      pny(is,js,ks) = s_lr * ny * pn
                      pnz(is,js,ks) = s_lr * nz * pn

                      if ( nvis == 1 ) then
                         call getgeo(is,js,ks)
                         call getuvwtder(is,js,ks)
                         call getduvwdxyz
                         call stress( visl(is,js,ks) )
                         pn = 2.0*s_lr
                         pnx(is,js,ks) = pnx(is,js,ks) - pn* ( nx*txx + ny*tyx + nz*tzx ) * re
                         pny(is,js,ks) = pny(is,js,ks) - pn* ( nx*txy + ny*tyy + nz*tzy ) * re
                         pnz(is,js,ks) = pnz(is,js,ks) - pn* ( nx*txz + ny*tyz + nz*tzz ) * re
                      endif
                   enddo
                enddo
             enddo

             do is=s_st(1),s_ed(1)-id
                do js=s_st(2),s_ed(2)-jd
                   do ks=s_st(3),s_ed(3)-kd
                      ii1 = is
                      jj1 = js
                      kk1 = ks
                      if ( id == 0 ) then
                         ii2 = is
                         jj2 = js+1
                         kk2 = ks

                         ii3 = is
                         jj3 = js+1
                         kk3 = ks+1

                         ii4 = is
                         jj4 = js
                         kk4 = ks+1
                      endif
                      if ( jd == 0 ) then
                         ii2 = is+1
                         jj2 = js
                         kk2 = ks

                         ii3 = is+1
                         jj3 = js
                         kk3 = ks+1

                         ii4 = is
                         jj4 = js
                         kk4 = ks+1
                      endif
                      if ( kd == 0 ) then
                         ii2 = is+1
                         jj2 = js
                         kk2 = ks

                         ii3 = is+1
                         jj3 = js+1
                         kk3 = ks

                         ii4 = is
                         jj4 = js+1
                         kk4 = ks
                      endif

                      px1 = pnx(ii1,jj1,kk1)
                      py1 = pny(ii1,jj1,kk1)
                      pz1 = pnz(ii1,jj1,kk1)

                      px2 = pnx(ii2,jj2,kk2)
                      py2 = pny(ii2,jj2,kk2)
                      pz2 = pnz(ii2,jj2,kk2)

                      px3 = pnx(ii3,jj3,kk3)
                      py3 = pny(ii3,jj3,kk3)
                      pz3 = pnz(ii3,jj3,kk3)

                      px4 = pnx(ii4,jj4,kk4)
                      py4 = pny(ii4,jj4,kk4)
                      pz4 = pnz(ii4,jj4,kk4)

                      x1 = x(ii1,jj1,kk1)
                      y1 = y(ii1,jj1,kk1)
                      z1 = z(ii1,jj1,kk1)

                      x2 = x(ii2,jj2,kk2)
                      y2 = y(ii2,jj2,kk2)
                      z2 = z(ii2,jj2,kk2)

                      x3 = x(ii3,jj3,kk3)
                      y3 = y(ii3,jj3,kk3)
                      z3 = z(ii3,jj3,kk3)

                      x4 = x(ii4,jj4,kk4)
                      y4 = y(ii4,jj4,kk4)
                      z4 = z(ii4,jj4,kk4)

                      call getforce(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,px4,py4,pz4, &
                                 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xref,yref,zref, &
                                 dfx,dfy,dfz,dmx,dmy,dmz)
                      cfx = cfx  + dfx
                      cfy = cfy  + dfy
                      cfz = cfz  + dfz

                      cmx = cmx  + dmx
                      cmy = cmy  + dmy
                      cmz = cmz  + dmz
                   enddo
                enddo
             enddo

             deallocate( pnx , pny , pnz )
             
          endif
       enddo
    enddo
    return
end subroutine force_dif
!_____________________________________________________________________!
subroutine getforce(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,px4,py4,pz4, &
                    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xref,yref,zref, &
                    dfx,dfy,dfz,dmx,dmy,dmz)
    implicit none 
    real :: px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,px4,py4,pz4,px5,py5,pz5
    real :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5
    real :: xref,yref,zref
    real :: dfx,dfy,dfz,dmx,dmy,dmz
    real :: r12,r23,r34,r41,r15,r25,r35,r45
    real :: ddfx,ddfy,ddfz,ddmx,ddmy,ddmz
    real :: pxc,pyc,pzc,xc,yc,zc
    real :: s1,s2,s3,s4

    px5 = 0.25 * ( px1 + px2 + px3 + px4 )
    py5 = 0.25 * ( py1 + py2 + py3 + py4 )
    pz5 = 0.25 * ( pz1 + pz2 + pz3 + pz4 )

    x5 = 0.25 * ( x1 + x2 + x3 + x4 )
    y5 = 0.25 * ( y1 + y2 + y3 + y4 )
    z5 = 0.25 * ( z1 + z2 + z3 + z4 )

    call lenth(x1,y1,z1,x2,y2,z2,r12)
    call lenth(x2,y2,z2,x3,y3,z3,r23)
    call lenth(x3,y3,z3,x4,y4,z4,r34)
    call lenth(x4,y4,z4,x1,y1,z1,r41)
    call lenth(x1,y1,z1,x5,y5,z5,r15)
    call lenth(x2,y2,z2,x5,y5,z5,r25)
    call lenth(x3,y3,z3,x5,y5,z5,r35)
    call lenth(x4,y4,z4,x5,y5,z5,r45)

    call halen(r12,r15,r25,s1)
    call halen(r23,r25,r35,s2)
    call halen(r34,r35,r45,s3)
    call halen(r41,r45,r15,s4)
    !��������,���س�ֵ
    dfx = 0.0
    dfy = 0.0
    dfz = 0.0

    dmx = 0.0
    dmy = 0.0
    dmz = 0.0

    !�Ե�һ�������ν�������������

    call aver3(px1,px2,px5,pxc)
    call aver3(py1,py2,py5,pyc)
    call aver3(pz1,pz2,pz5,pzc)
    call aver3(x1 ,x2 ,x5 ,xc )
    call aver3(y1 ,y2 ,y5 ,yc )
    call aver3(z1 ,z2 ,z5 ,zc )

    ddfx = pxc * s1
    ddfy = pyc * s1
    ddfz = pzc * s1

    ddmx = ( yc - yref ) * ddfz - ( zc - zref ) * ddfy
    ddmy = ( zc - zref ) * ddfx - ( xc - xref ) * ddfz
    ddmz = ( xc - xref ) * ddfy - ( yc - yref ) * ddfx

    dfx = dfx + ddfx
    dfy = dfy + ddfy
    dfz = dfz + ddfz

    dmx = dmx + ddmx
    dmy = dmy + ddmy
    dmz = dmz + ddmz

    !�Եڶ��������ν�������������

    call aver3(px2,px3,px5,pxc)
    call aver3(py2,py3,py5,pyc)
    call aver3(pz2,pz3,pz5,pzc)
    call aver3(x2 ,x3 ,x5 ,xc )
    call aver3(y2 ,y3 ,y5 ,yc )
    call aver3(z2 ,z3 ,z5 ,zc )

    ddfx = pxc * s2
    ddfy = pyc * s2
    ddfz = pzc * s2

    ddmx = ( yc - yref ) * ddfz - ( zc - zref ) * ddfy
    ddmy = ( zc - zref ) * ddfx - ( xc - xref ) * ddfz
    ddmz = ( xc - xref ) * ddfy - ( yc - yref ) * ddfx

    dfx = dfx + ddfx
    dfy = dfy + ddfy
    dfz = dfz + ddfz

    dmx = dmx + ddmx
    dmy = dmy + ddmy
    dmz = dmz + ddmz

    !�Ե����������ν�������������

    call aver3(px3,px4,px5,pxc)
    call aver3(py3,py4,py5,pyc)
    call aver3(pz3,pz4,pz5,pzc)
    call aver3(x3 ,x4 ,x5 ,xc )
    call aver3(y3 ,y4 ,y5 ,yc )
    call aver3(z3 ,z4 ,z5 ,zc )

    ddfx = pxc * s3
    ddfy = pyc * s3
    ddfz = pzc * s3

    ddmx = ( yc - yref ) * ddfz - ( zc - zref ) * ddfy
    ddmy = ( zc - zref ) * ddfx - ( xc - xref ) * ddfz
    ddmz = ( xc - xref ) * ddfy - ( yc - yref ) * ddfx

    dfx = dfx + ddfx
    dfy = dfy + ddfy
    dfz = dfz + ddfz

    dmx = dmx + ddmx
    dmy = dmy + ddmy
    dmz = dmz + ddmz

    !�Ե��ĸ������ν�������������

    call aver3(px4,px1,px5,pxc)
    call aver3(py4,py1,py5,pyc)
    call aver3(pz4,pz1,pz5,pzc)
    call aver3(x4 ,x1 ,x5 ,xc )
    call aver3(y4 ,y1 ,y5 ,yc )
    call aver3(z4 ,z1 ,z5 ,zc )

    ddfx = pxc * s4
    ddfy = pyc * s4
    ddfz = pzc * s4

    ddmx = ( yc - yref ) * ddfz - ( zc - zref ) * ddfy
    ddmy = ( zc - zref ) * ddfx - ( xc - xref ) * ddfz
    ddmz = ( xc - xref ) * ddfy - ( yc - yref ) * ddfx

    dfx = dfx + ddfx
    dfy = dfy + ddfy
    dfz = dfz + ddfz

    dmx = dmx + ddmx
    dmy = dmy + ddmy
    dmz = dmz + ddmz

    return
end subroutine getforce
!_____________________________________________________________________!
subroutine heatflux_dif
    !���ӳ��������������
    use global_variables, &
    only : mb_kcx,mb_kcy,mb_kcz,mb_etx,mb_ety,mb_etz,mb_ctx,mb_cty,mb_ctz, &
           mb_x,mb_y,mb_z,mb_bc,nblocks,visl,vist,t,re,ni,nj,nk,nl,nm,ns, &
           gama,prl,prt,moo,nchem_source,nchem,fs,tref,twall,tecname
    use geometry_module
    use duvwt_module,only : tkc,tet,tct
    implicit none
    integer :: m,nb,nr,nrmax,bctype,nsurf,lr
    integer :: s_lr3d(3)
    integer :: i,j,k,ic,jc,kc,id,jd,kd,im,jm,km,ntec
    integer :: i1,j1,k1,i2,j2,k2
    integer :: ip,jp,kp
    integer :: ist,ied,jst,jed,kst,ked
    real    :: nx,ny,nz,on,xc,yc,zc
    real    :: vis_l,vis_t,q0_mml
    real    :: cp,kcp,dtdx,dtdy,dtdz,qx,qy,qz,qodd
    real    :: cs_init(ns),coef,tstoo
    real,pointer,dimension( :,:,:,: ) :: qflux
    character(len=120) :: title_qflux,heatfluxfile
    !integer :: im1,jm1,km1

    ip = len(trim(tecname))
    heatfluxfile=tecname(1:ip-4)//'_heat.dat'

    tstoo = 1.0 + 0.5 * ( gama - 1.0 ) * moo * moo
    coef = ( gama - 1.0 )*moo*moo/(tstoo-twall/tref)

    ntec = 4
    title_qflux ='variables="x" "y" "z" "heatflux" '
    open(1,file=heatfluxfile,status='unknown')
    do nb=1,nblocks
       call recast_grid (nb)
       call recast_field(nb)
       nrmax = mb_bc(nb)%nregions
       do nr=1,nrmax
          bctype = mb_bc(nb)%bc(nr)%bctype
          if ( bctype/10 == 2 ) then 
             ist = mb_bc(nb)%bc(nr)%s_st(1)
             ied = mb_bc(nb)%bc(nr)%s_ed(1)
             jst = mb_bc(nb)%bc(nr)%s_st(2)
             jed = mb_bc(nb)%bc(nr)%s_ed(2)
             kst = mb_bc(nb)%bc(nr)%s_st(3)
             ked = mb_bc(nb)%bc(nr)%s_ed(3)

						 allocate( qflux(ist:ied,jst:jed,kst:ked,1:ntec) )

             s_lr3d = mb_bc(nb)%bc(nr)%s_lr3d

             nsurf  = mb_bc(nb)%bc(nr)%s_nd         
             lr     = mb_bc(nb)%bc(nr)%s_lr          

             do k=kst,ked
                do j=jst,jed
                   do i=ist,ied

                      if ( nsurf == 1 ) then
                         nx = mb_kcx(nb)%a3d(i,j,k)
                         ny = mb_kcy(nb)%a3d(i,j,k)
                         nz = mb_kcz(nb)%a3d(i,j,k)
                      elseif ( nsurf == 2 ) then
                         nx = mb_etx(nb)%a3d(i,j,k)
                         ny = mb_ety(nb)%a3d(i,j,k)
                         nz = mb_etz(nb)%a3d(i,j,k)
                      else
                         nx = mb_ctx(nb)%a3d(i,j,k)
                         ny = mb_cty(nb)%a3d(i,j,k)
                         nz = mb_ctz(nb)%a3d(i,j,k)
                      endif

                      on = 1.0/( sqrt ( nx * nx + ny * ny + nz * nz ) + 1.0e-30 )
                      nx = - lr * nx * on
                      ny = - lr * ny * on
                      nz = - lr * nz * on

                      call getgeo(i,j,k)
                      call getuvwtder(i,j,k)
                      dtdx = (kx*tkc + ex*tet + cx*tct)/vjacob
                      dtdy = (ky*tkc + ey*tet + cy*tct)/vjacob
                      dtdz = (kz*tkc + ez*tet + cz*tct)/vjacob

                      im = i + s_lr3d(1)
                      jm = j + s_lr3d(2)
                      km = k + s_lr3d(3)

                      vis_l = visl(i,j,k)
                      vis_t = 0.0                !!!!!!!!!!!!!!!!

                      qx = 0.0
                      qy = 0.0
                      qz = 0.0

                      cp = 1.0/( (gama-1.0) * moo * moo )

                      kcp = ( vis_l/prl + vis_t/prt ) * cp

                      qx = qx + kcp * dtdx
                      qy = qy + kcp * dtdy
                      qz = qz + kcp * dtdz

                      qflux(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                      qflux(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                      qflux(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
!                      qflux(i,j,k,4) = re*coef*( nx * qx + ny * qy + nz * qz )
                      qflux(i,j,k,4) = re*( nx * qx + ny * qy + nz * qz )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                enddo
             enddo

             !�޸����
             if ( nsurf == 1 ) then
                i  = ist
                k  = kst
                k1 = kst + 1
                k2 = kst + 2
                qodd = 0.0
                do j=jst,jed
                   if( qflux(i,j,k,1) /= qflux(i,jst,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,jst,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,jst,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do j=jst,jed
                      qodd = qodd + ( 4.0 * qflux(i,j,k1,4) - qflux(i,j,k2,4) ) / 3.0
                   enddo
                   do j=jst,jed
                      qflux(i,j,k,4) = qodd*coef/ float( jed - jst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                k  = ked
                k1 = kst-1
                k2 = kst-2
                qodd = 0.0
                do j=jst,jed
                   if( qflux(i,j,k,1) /= qflux(i,jst,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,jst,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,jst,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do j=jst,jed
                      qodd = qodd + ( 4.0 * qflux(i,j,k1,4) - qflux(i,j,k2,4) ) / 3.0
                   enddo
                   do j=jst,jed
                      qflux(i,j,k,4) = qodd*coef/ float( jed - jst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                j  = jst
                j1 = jst + 1
                j2 = jst + 2
                qodd = 0.0
                do k=kst,ked
                   if( qflux(i,j,k,1) /= qflux(i,j,kst,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,j,kst,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,j,kst,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do k=kst,ked
                      qodd = qodd + ( 4.0 * qflux(i,j1,k,4) - qflux(i,j2,k,4) ) / 3.0
                   enddo
                   do k=kst,ked
                      qflux(i,j,k,4) = qodd*coef/ float( ked - kst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                j  = jed
                j1 = jed - 1
                j2 = jed - 2
                qodd = 0.0
                do k=kst,ked
                   if( qflux(i,j,k,1) /= qflux(i,j,kst,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,j,kst,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,j,kst,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do k=kst,ked
                      qodd = qodd + ( 4.0 * qflux(i,j1,k,4) - qflux(i,j2,k,4) ) / 3.0
                   enddo
                   do k=kst,ked
                      qflux(i,j,k,4) = qodd*coef/ float( ked - kst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
             elseif ( nsurf == 2 ) then
                j  = jst
                k  = kst
                k1 = kst + 1
                k2 = kst + 2
                qodd = 0.0
                do i=ist,ied
                   if( qflux(i,j,k,1) /= qflux(ist,j,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(ist,j,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(ist,j,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do i=ist,ied
                      qodd = qodd + ( 4.0 * qflux(i,j,k1,4) - qflux(i,j,k2,4) ) / 3.0
                   enddo
                   do i=ist,ied
                      qflux(i,j,k,4) = qodd*coef/ float( ied - ist + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                k  = ked
                k1 = ked - 1
                k2 = ked - 2
                qodd = 0.0
                do i=ist,ied
                   if( qflux(i,j,k,1) /= qflux(ist,j,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(ist,j,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(ist,j,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do i=ist,ied
                      qodd = qodd + ( 4.0 * qflux(i,j,k1,4) - qflux(i,j,k2,4) ) / 3.0
                   enddo
                   do i=ist,ied
                      qflux(i,j,k,4) = qodd*coef/ float( ied - ist + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                i  = ist
                i1 = ist + 1
                i2 = ist + 2
                qodd = 0.0
                do k=kst,ked
                   if( qflux(i,j,k,1) /= qflux(i,j,kst,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,j,kst,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,j,kst,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do k=kst,ked
                      qodd = qodd + ( 4.0 * qflux(i1,j,k,4) - qflux(i2,j,k,4) ) / 3.0
                   enddo
                   do k=kst,ked
                      qflux(i,j,k,4) = qodd*coef/ float( ked - kst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                i  = ied
                i1 = ied - 1
                i2 = ied - 2
                qodd = 0.0
                do k=kst,ked
                   if( qflux(i,j,k,1) /= qflux(i,j,kst,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,j,kst,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,j,kst,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do k=kst,ked
                      qodd = qodd + ( 4.0 * qflux(i1,j,k,4) - qflux(i2,j,k,4) ) / 3.0
                   enddo
                   do k=kst,ked
                      qflux(i,j,k,4) = qodd*coef/ float( ked - kst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
             else
                k  = kst
                j  = jst
                j1 = jst + 1
                j2 = jst + 2
                qodd = 0.0
                do i=ist,ied
                   if( qflux(i,j,k,1) /= qflux(ist,j,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(ist,j,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(ist,j,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do i=ist,ied
                      qodd = qodd + ( 4.0 * qflux(i,j1,k,4) - qflux(i,j2,k,4) ) / 3.0
                   enddo
                   do i=ist,ied
                      qflux(i,j,k,4) = qodd*coef/ float( ied - ist + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                j  = jed
                j1 = jed - 1
                j2 = jed - 2
                qodd = 0.0
                do i=ist,ied
                   if( qflux(i,j,k,1) /= qflux(ist,j,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(ist,j,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(ist,j,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do i=ist,ied
                      qodd = qodd + ( 4.0 * qflux(i,j1,k,4) - qflux(i,j2,k,4) ) / 3.0
                   enddo
                   do i=ist,ied
                      qflux(i,j,k,4) = qodd*coef/ float( ied - ist + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                i  = ist
                i1 = ist + 1
                i2 = ist + 2
                qodd = 0.0
                do j=jst,jed
                   if( qflux(i,j,k,1) /= qflux(i,jst,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,jst,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,jst,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do j=jst,jed
                      qodd = qodd + ( 4.0 * qflux(i1,j,k,4) - qflux(i2,j,k,4) ) / 3.0
                   enddo
                   do j=jst,jed
                      qflux(i,j,k,4) = qodd*coef/ float( jed - jst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
                i  = ied
                i1 = ied - 1
                i2 = ied - 2
                qodd = 0.0
                do j=jst,jed
                   if( qflux(i,j,k,1) /= qflux(i,jst,k,1) .or. &
                       qflux(i,j,k,2) /= qflux(i,jst,k,2) .or. &
                       qflux(i,j,k,3) /= qflux(i,jst,k,3) ) qodd = 1.0
                enddo
                if ( qodd == 0.0 ) then
                   do j=jst,jed
                      qodd = qodd + ( 4.0 * qflux(i1,j,k,4) - qflux(i2,j,k,4) ) / 3.0
                   enddo
                   do j=jst,jed
                      qflux(i,j,k,4) = qodd*coef/ float( jed - jst + 1 )
!                     qflux(i,j,k,5) = qflux(i,j,k,4)*coef
                   enddo
                endif
             endif

             write(1,*)trim(title_qflux)
             write(1,*)'zone i=',ied-ist+1,'j=',jed-jst+1,'k=',ked-kst+1,' f=block'
             write(1,*)((((qflux(i,j,k,m),i=ist,ied),j=jst,jed),k=kst,ked),m=1,ntec)
						deallocate( qflux )
         endif
       enddo
    enddo
    close(1)

    return
end subroutine heatflux_dif
!------------------------------------------------------------------------!
!------------------------------------------------------------------------!
subroutine tecout_bin_2d
    use global_variables,only : nlamtur,nblocks,tecname,gama,nout,mb_dim &
            ,mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,r,u,v,w,p,t,mb_vist &
            ,mb_x,mb_y,mb_z,x,y,z,mb_qke
    use mod_tecios
#ifdef PARALLEL
    use mod_parallels
#endif		
    implicit none
    integer :: nb,nbt
    integer :: ibeg(3),iend(3),pid
    integer :: i,j,k,m,ndata,ierr
    integer :: idim,jdim,kdim
#ifdef PARALLEL
    integer :: packsize
    integer :: status(MPI_STATUS_SIZE)
#endif
    logical :: ldraw
    integer :: ntecvars
    character(len=256) :: tecvarnames
    character(len=120) :: zonename,tecvolname
    character(len=32 ) :: str1
    real,pointer :: qpv_tmp(:,:,:,:)
    real :: c2,u2
    integer single_prec,prec
    
    prec=kind(c2)
    single_prec=4
    
    
    ntecvars = 10
    tecvarnames = "X,Y,Z,R,U,V,W,P,T,M"
    
    if (nlamtur >= 0) then
       ntecvars = ntecvars + nlamtur + 1
       select case(nlamtur)
       case(0)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,muT"
       case(1)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,muT"
       case(2)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,vT2,muT"
       end select
    end if
    
#ifdef PARALLEL
    if (myid == master) then
#endif
       m = index(tecname,'.',back=.TRUE.)
       if (m == 0) then
          tecvolname = trim(adjustl(tecname)) // "_2d.plt"
       else
          tecvolname = trim(adjustl(tecname(1:m-1))) // "_2d.plt"
       end if

       open(101,file=tecvolname,form='unformatted', access='stream',status='unknown')
    
       call tecio_ini(101,"WCNS_SOLVER",ntecvars,tecvarnames)
      
       do nb=1,nblocks
          ibeg = 1
          iend = mb_dim(nb,:)
          iend(3) = 1
          idim = iend(1)-ibeg(1)+1
          jdim = iend(2)-ibeg(2)+1
          kdim = iend(3)-ibeg(3)+1
          write(str1,'(i6)') nb
          
          zonename = 'BLK'//trim(adjustl(str1))

          call tecio_zone(101,trim(zonename),0,idim,jdim,kdim)
                
       end do
       
       call tecio_eohmark(101)
       
       if (prec == single_prec) then
          ndata = 1
       else
          ndata = 2
       end if
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
       ibeg = 1
       iend = mb_dim(nb,:)
       iend(3) = 1
             
#ifdef PARALLEL
       pid = mb_pids(nb)
       if (master == pid-1) then
          if (myid == master) then
#endif          
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
        
             call tecio_data(101,ndata)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
             
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                       
             deallocate(qpv_tmp,stat=ierr)
#ifdef PARALLEL                   
          end if
       else
          if (myid == master) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
                              
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                   
             call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                           pid-1,nb,MPI_COMM_WORLD,status,ierr)
                           
             call tecio_data(101,ndata)
                           
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
             
             deallocate(qpv_tmp,stat=ierr)
          endif
          
          if (myid == pid-1) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
            
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
             
             call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                           master,nb,MPI_COMM_WORLD,ierr)
             
             deallocate(qpv_tmp,stat=ierr)
          end if                
       end if
#endif             
    end do
    
#ifdef PARALLEL                   
    if (myid == master) then
#endif
       close(101)
#ifdef PARALLEL                   
    end if
#endif
        
end subroutine tecout_bin_2d
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine tecout_bin_2d_serial
    use global_variables,only : nlamtur,nblocks,tecname,gama,nout,mb_dim &
            ,mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,r,u,v,w,p,t,mb_vist &
            ,mb_x,mb_y,mb_z,x,y,z,mb_qke
    use mod_tecios
#ifdef PARALLEL
    use mod_parallels
#endif		
    implicit none
    integer :: nb,nbt
    integer :: ibeg(3),iend(3),pid
    integer :: i,j,k,m,ndata,ierr
    integer :: idim,jdim,kdim
#ifdef PARALLEL
    integer :: packsize
    integer :: status(MPI_STATUS_SIZE)
#endif
    logical :: ldraw
    integer :: ntecvars
    character(len=256) :: tecvarnames
    character(len=120) :: zonename,tecvolname
    character(len=32 ) :: str1
    character(len=7 )  :: str_step
	integer :: teclength,step,step0,step1,step2,step3,step4,step5,step6
    real,pointer :: qpv_tmp(:,:,:,:)
    real :: c2,u2
    integer single_prec,prec
    
	str_step = '0000000'
	step =nout
	step6 = step/1000000
	step5 = (step - step6*1000000)/100000
	step4 = (step - step6*1000000 - step5*100000)/10000
	step3 = (step - step6*1000000 - step5*100000 - step4*10000)/1000
	step2 = (step - step6*1000000 - step5*100000 - step4*10000 - step3*1000)/100
	step1 = (step - step6*1000000 - step5*100000 - step4*10000 - step3*1000 - step2*100)/10
	step0 =  step - step6*1000000 - step5*100000 - step4*10000 - step3*1000 - step2*100 - step1*10
	str_step = char(48 + step6)//char(48 + step5)//char(48 + step4)//char(48 + step3)// &
      char(48 + step2)//char(48 + step1)//char(48 + step0)

    prec=kind(c2)
    single_prec=4
    
    
    ntecvars = 10
    tecvarnames = "X,Y,Z,R,U,V,W,P,T,M"
    
    if (nlamtur >= 0) then
       ntecvars = ntecvars + nlamtur + 1
       select case(nlamtur)
       case(0)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,muT"
       case(1)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,muT"
       case(2)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,vT2,muT"
       end select
    end if
    
#ifdef PARALLEL
    if (myid == master) then
#endif
       m = index(tecname,'.',back=.TRUE.)
       if (m == 0) then
          tecvolname = trim(adjustl(tecname))//str_step// ".plt"
       else
          tecvolname = trim(adjustl(tecname(1:m-1)))//str_step// ".plt"
       end if
  
       open(101,file=tecvolname,form='unformatted', access='stream',status='unknown')
   
       call tecio_ini(101,"WCNS_SOLVER",ntecvars,tecvarnames)
    
       do nb=1,nblocks
          ibeg = 1
          iend = mb_dim(nb,:)
          iend(3) = 1
          idim = iend(1)-ibeg(1)+1
          jdim = iend(2)-ibeg(2)+1
          kdim = iend(3)-ibeg(3)+1
          write(str1,'(i6)') nb
          
          zonename = 'BLK'//trim(adjustl(str1))

          call tecio_zone(101,trim(zonename),0,idim,jdim,kdim)
                
       end do
       
       call tecio_eohmark(101)
       
       if (prec == single_prec) then
          ndata = 1
       else
          ndata = 2
       end if
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
       ibeg = 1
       iend = mb_dim(nb,:)
       iend(3) = 1
             
#ifdef PARALLEL
       pid = mb_pids(nb)
       if (master == pid-1) then
          if (myid == master) then
#endif          
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
        
             call tecio_data(101,ndata)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
             
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                       
             deallocate(qpv_tmp,stat=ierr)
#ifdef PARALLEL                   
          end if
       else
          if (myid == master) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
                              
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                   
             call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                           pid-1,nb,MPI_COMM_WORLD,status,ierr)
                           
             call tecio_data(101,ndata)
                           
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
             
             deallocate(qpv_tmp,stat=ierr)
          endif
          
          if (myid == pid-1) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
            
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
             
             call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                           master,nb,MPI_COMM_WORLD,ierr)
             
             deallocate(qpv_tmp,stat=ierr)
          end if                
       end if
#endif             
    end do
    
#ifdef PARALLEL                   
    if (myid == master) then
#endif
       close(101)
#ifdef PARALLEL                   
    end if
#endif
        
end subroutine tecout_bin_2d_serial
!_____________________________________________________________________!

!_____________________________________________________________________!
subroutine tecout_bin_3d_serial
    use global_variables,only : nlamtur,nblocks,tecname,gama,nout,mb_dim &
            ,mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,r,u,v,w,p,t,mb_vist &
            ,mb_x,mb_y,mb_z,x,y,z,mb_qke
    use mod_tecios
#ifdef PARALLEL
    use mod_parallels
#endif		
    implicit none
    integer :: nb,nbt
    integer :: ibeg(3),iend(3),pid
    integer :: i,j,k,m,ndata,ierr
    integer :: idim,jdim,kdim
#ifdef PARALLEL
    integer :: packsize
    integer :: status(MPI_STATUS_SIZE)
#endif
    logical :: ldraw
    integer :: ntecvars
    character(len=256) :: tecvarnames
    character(len=120) :: zonename,tecvolname
    character(len=32 ) :: str1
    character(len=7 )  :: str_step
	integer :: teclength,step,step0,step1,step2,step3,step4,step5,step6
    real,pointer :: qpv_tmp(:,:,:,:)
    real :: c2,u2
    integer single_prec,prec
    
	str_step = '0000000'
	step =nout
	step6 = step/1000000
	step5 = (step - step6*1000000)/100000
	step4 = (step - step6*1000000 - step5*100000)/10000
	step3 = (step - step6*1000000 - step5*100000 - step4*10000)/1000
	step2 = (step - step6*1000000 - step5*100000 - step4*10000 - step3*1000)/100
	step1 = (step - step6*1000000 - step5*100000 - step4*10000 - step3*1000 - step2*100)/10
	step0 =  step - step6*1000000 - step5*100000 - step4*10000 - step3*1000 - step2*100 - step1*10
	str_step = char(48 + step6)//char(48 + step5)//char(48 + step4)//char(48 + step3)// &
      char(48 + step2)//char(48 + step1)//char(48 + step0)

    prec=kind(c2)
    single_prec=4
    
    
    ntecvars = 10
    tecvarnames = "X,Y,Z,R,U,V,W,P,T,M"
    
    if (nlamtur >= 0) then
       ntecvars = ntecvars + nlamtur + 1
       select case(nlamtur)
       case(0)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,muT"
       case(1)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,muT"
       case(2)
          tecvarnames = "X,Y,Z,R,U,V,W,P,T,M,vT1,vT2,muT"
       end select
    end if
    
#ifdef PARALLEL
    if (myid == master) then
#endif
       m = index(tecname,'.',back=.TRUE.)
       if (m == 0) then
          tecvolname = trim(adjustl(tecname))//str_step// ".plt"
       else
          tecvolname = trim(adjustl(tecname(1:m-1)))//str_step// ".plt"
       end if
  
       open(101,file=tecvolname,form='unformatted', access='stream',status='unknown')
   
       call tecio_ini(101,"WCNS_SOLVER",ntecvars,tecvarnames)
    
       do nb=1,nblocks
          ibeg = 1
          iend = mb_dim(nb,:)
!!          iend(3) = 1
          idim = iend(1)-ibeg(1)+1
          jdim = iend(2)-ibeg(2)+1
          kdim = iend(3)-ibeg(3)+1
          write(str1,'(i6)') nb
          
          zonename = 'BLK'//trim(adjustl(str1))

          call tecio_zone(101,trim(zonename),0,idim,jdim,kdim)
                
       end do
       
       call tecio_eohmark(101)
       
       if (prec == single_prec) then
          ndata = 1
       else
          ndata = 2
       end if
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
       ibeg = 1
       iend = mb_dim(nb,:)
!!       iend(3) = 1
             
#ifdef PARALLEL
       pid = mb_pids(nb)
       if (master == pid-1) then
          if (myid == master) then
#endif          
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
        
             call tecio_data(101,ndata)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
             
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
                       
             deallocate(qpv_tmp,stat=ierr)
#ifdef PARALLEL                   
          end if
       else
          if (myid == master) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
                              
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
                   
             call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                           pid-1,nb,MPI_COMM_WORLD,status,ierr)
                           
             call tecio_data(101,ndata)
                           
             write(101)((((qpv_tmp(i,j,k,m),i=ibeg(1),iend(1)),j=ibeg(2),iend(2)),k=ibeg(3),iend(3)),m=1,ntecvars)
             
             deallocate(qpv_tmp,stat=ierr)
          endif
          
          if (myid == pid-1) then
             allocate(qpv_tmp(ibeg(1):iend(1),ibeg(2):iend(2), &
                              ibeg(3):iend(3),ntecvars),stat=ierr)
             
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                qpv_tmp(i,j,k,1) = mb_x(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,2) = mb_y(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,3) = mb_z(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,4) = mb_r(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,5) = mb_u(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,6) = mb_v(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,7) = mb_w(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,8) = mb_p(nb)%a3d(i,j,k)
                qpv_tmp(i,j,k,9) = mb_t(nb)%a3d(i,j,k)
                c2 = gama*qpv_tmp(i,j,k,8)/qpv_tmp(i,j,k,4)
                u2 = qpv_tmp(i,j,k,5)**2 + qpv_tmp(i,j,k,6)**2 + qpv_tmp(i,j,k,7)**2
                qpv_tmp(i,j,k,10) = sqrt(u2/c2)
                
                if (nlamtur >= 0) then
                   do m=1,nlamtur
                      qpv_tmp(i,j,k,10+m) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(i,j,k,ntecvars) = mb_vist(nb)%a3d(i,j,k)
                end if
             end do
             end do
             end do
            
             packsize = product(iend(:)-ibeg(:)+1,1)*ntecvars
             
             call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                           master,nb,MPI_COMM_WORLD,ierr)
             
             deallocate(qpv_tmp,stat=ierr)
          end if                
       end if
#endif             
    end do
    
#ifdef PARALLEL                   
    if (myid == master) then
#endif
       close(101)
#ifdef PARALLEL                   
    end if
#endif
        
end subroutine tecout_bin_3d_serial
!_____________________________________________________________________!
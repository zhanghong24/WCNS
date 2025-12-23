!_____________________________________________________________________!
subroutine pre_processing
    use global_variables,only: nvis,nlhs,nptt,nscheme
   use mod_parallels,only : master,myid,distrib_parameter,read_bc_parallel, &
                             distrib_grid,distrib_bc,build_nppos_list, &
                             build_nppos_list_ex,exchange_bc_geometric
    implicit none
    
    call init0

   if (myid == master) then
    call read_parameter
    end if

    call distrib_parameter

    call write_param_dump

   if (myid == master) then
       call read_plot3d_grid 
       call read_bc_parallel
    end if

    call distrib_grid 

    call distrib_bc

    call set_bc_index !here

    call allocate_other_variable
    
    call set_grid_derivative
    
    call check_grid_derivative
    
    call init_inflow
    
    call read_bc_par
    
    call initialization

#define CONNECT_PRE_EX
#ifdef PARALLEL
#ifndef CONNECT_PRE_EX

    call connect_pre_new
#endif
#else

    call connect_pre_new
#endif

#ifdef PARALLEL
#ifdef CONNECT_PRE_EX
    call build_nppos_list_ex
#else
    call build_nppos_list
#endif
    call exchange_bc_geometric
#endif

	if(nvis==1)CALL INITIAL_T  !*tgh. �����ʼ�¶�

    return
end subroutine pre_processing

!_____________________________________________________________________!
subroutine init0
    use global_const,only:cyc,iii
    implicit none
    cyc(1,1) = 2 ; cyc(1,2) = 3
    cyc(2,1) = 3 ; cyc(2,2) = 1
    cyc(3,1) = 1 ; cyc(3,2) = 2
    iii(1,1) = 1 ; iii(1,2) = 0 ; iii(1,3) = 0
    iii(2,1) = 0 ; iii(2,2) = 1 ; iii(2,3) = 0
    iii(3,1) = 0 ; iii(3,2) = 0 ; iii(3,3) = 1
    return
end subroutine init0
!_____________________________________________________________________!
subroutine read_parameter
!--------
    use global_variables
    implicit none
    integer :: nfile_par
    namelist /inflow/moo,reynolds,attack,sideslip,tref,twall,pref,rref,vref,ndim,height
    namelist /control/nstart,nmethod,nbgmax,ndisk,newton,nomax,nforce,nwerror,method,nplot
    namelist /step/ntmst,cfl,timedt,timedt_rate,timedt_turb,ndualtst,dtdts,nsubstmx,tolsub
    namelist /flowtype/nvis,nchem,ntmodel
    namelist /turbulence/nlamtur,nameturb,nreset,nwallfun,ntrans,xtrans,nwtmax,ncmpcor,kmaxlim,ndes
    namelist /technic/nlhs,nscheme,nlimiter,efix,csrv,nsmooth,nflux
    namelist /interplate/xk,xb,c_k2,c_k4
    namelist /filename/flowname,gridname,bcname,forcename,errhis,tecname
    namelist /force/nwholefield,sref,lfref,lref,xref,yref,zref
    namelist /chemical/gasmodel,nchem_source,nchem_rad

    namelist /gcl_cic/gcl,cic1,cic2,cic3,cic4,cic5
    namelist /connect_0/connect_point,conpointname,connect_order,dis_tao

    open(1,file='param.dat',form='formatted',status='old')
    
    read(1,nml=inflow)
#ifdef OUT_INPUT_NAMELIST    
    write(*,nml=inflow)
#endif		

    read(1,nml=force)
    
    read(1,nml=filename)
#ifdef OUT_INPUT_NAMELIST    
    write(*,nml=filename)
#endif		
		
    read(1,nml=control)
#ifdef OUT_INPUT_NAMELIST    
    write(*,nml=control)
#endif		

       read(1,nml=step)
       
       read(1,nml=flowtype)
#ifdef OUT_INPUT_NAMELIST    
       write(*,nml=flowtype)
#endif		

       read(1,nml=turbulence)
#ifdef OUT_INPUT_NAMELIST    
       write(*,nml=turbulence)
#endif		
       
       read(1,nml=technic)
#ifdef OUT_INPUT_NAMELIST    
       write(*,nml=technic)
#endif		

       read(1,nml=interplate)
       read(1,nml=chemical)
       if( nchem_source == 0 ) nchem_rad = 0

    read(1,nml=gcl_cic)
#ifdef OUT_INPUT_NAMELIST    
    write(*,nml=gcl_cic)
#endif		
	  
    read(1,nml=connect_0)
#ifdef OUT_INPUT_NAMELIST
    write(*,nml=connect_0)
#endif

    close(1)
    
    if (ndualtst > 0) then
       if (nlhs == 6) call stop_by_error(1000,"method of time integration is incompatible with implicit method")
    end if
    
       call read_perfect_gasmodel !( nfile_par )

    return
!--------
end subroutine read_parameter
!_____________________________________________________________________!
subroutine write_param_dump()
    use global_variables
#ifdef PARALLEL
    use mod_parallels, only: myid, master
#endif
    implicit none
    integer :: unitno, ierr
    logical :: do_write

#ifdef PARALLEL
    do_write = (myid == master)
#else
    do_write = .true.
#endif

    if (.not. do_write) return

    call execute_command_line('mkdir -p output', wait=.true., exitstat=ierr)
    if (ierr /= 0) then
       call stop_by_error(9001, 'failed to create output directory for param_dump')
    end if

    open(newunit=unitno, file='output/param_dump.txt', status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
       call stop_by_error(9002, 'failed to open output/param_dump.txt')
    end if

    call write_real_entry(unitno, 'moo', moo)
    call write_real_entry(unitno, 'reynolds', reynolds)
    call write_real_entry(unitno, 'attack', attack)
    call write_real_entry(unitno, 'sideslip', sideslip)
    call write_real_entry(unitno, 'tref', tref)
    call write_real_entry(unitno, 'twall', twall)
    call write_real_entry(unitno, 'pref', pref)
    call write_real_entry(unitno, 'rref', rref)
    call write_real_entry(unitno, 'vref', vref)
    call write_integer_entry(unitno, 'ndim', ndim)
    call write_real_entry(unitno, 'height', height)

    call write_integer_entry(unitno, 'nwholefield', nwholefield)
    call write_real_entry(unitno, 'sref', sref)
    call write_real_entry(unitno, 'lfref', lfref)
    call write_real_entry(unitno, 'lref', lref)
    call write_real_entry(unitno, 'xref', xref)
    call write_real_entry(unitno, 'yref', yref)
    call write_real_entry(unitno, 'zref', zref)

    call write_string_entry(unitno, 'flowname', flowname)
    call write_string_entry(unitno, 'gridname', gridname)
    call write_string_entry(unitno, 'bcname', bcname)
    call write_string_entry(unitno, 'forcename', forcename)
    call write_string_entry(unitno, 'errhis', errhis)
    call write_string_entry(unitno, 'tecname', tecname)

    call write_integer_entry(unitno, 'nstart', nstart)
    call write_integer_entry(unitno, 'nmethod', nmethod)
    call write_integer_entry(unitno, 'nbgmax', nbgmax)
    call write_integer_entry(unitno, 'ndisk', ndisk)
    call write_integer_entry(unitno, 'newton', newton)
    call write_integer_entry(unitno, 'nomax', nomax)
    call write_integer_entry(unitno, 'nforce', nforce)
    call write_integer_entry(unitno, 'nwerror', nwerror)
    call write_integer_entry(unitno, 'method', method)
    call write_integer_entry(unitno, 'nplot', nplot)

    call write_integer_entry(unitno, 'ntmst', ntmst)
    call write_real_entry(unitno, 'cfl', cfl)
    call write_real_entry(unitno, 'timedt', timedt)
    call write_real_entry(unitno, 'timedt_rate', timedt_rate)
    call write_real_entry(unitno, 'timedt_turb', timedt_turb)
    call write_integer_entry(unitno, 'ndualtst', ndualtst)
    call write_real_entry(unitno, 'dtdts', dtdts)
    call write_integer_entry(unitno, 'nsubstmx', nsubstmx)
    call write_real_entry(unitno, 'tolsub', tolsub)

    call write_integer_entry(unitno, 'nvis', nvis)
    call write_integer_entry(unitno, 'nchem', nchem)
    call write_integer_entry(unitno, 'ntmodel', ntmodel)

    call write_integer_entry(unitno, 'nlamtur', nlamtur)
    call write_string_entry(unitno, 'nameturb', nameturb)
    call write_integer_entry(unitno, 'nreset', nreset)
    call write_integer_entry(unitno, 'nwallfun', nwallfun)
    call write_integer_entry(unitno, 'ntrans', ntrans)
    call write_real_entry(unitno, 'xtrans', xtrans)
    call write_integer_entry(unitno, 'nwtmax', nwtmax)
    call write_integer_entry(unitno, 'ncmpcor', ncmpcor)
    call write_real_entry(unitno, 'kmaxlim', kmaxlim)
    call write_integer_entry(unitno, 'ndes', ndes)

    call write_integer_entry(unitno, 'nlhs', nlhs)
    call write_integer_entry(unitno, 'nscheme', nscheme)
    call write_integer_entry(unitno, 'nlimiter', nlimiter)
    call write_real_entry(unitno, 'efix', efix)
    call write_real_entry(unitno, 'csrv', csrv)
    call write_integer_entry(unitno, 'nsmooth', nsmooth)
    call write_integer_entry(unitno, 'nflux', nflux)

    call write_real_entry(unitno, 'xk', xk)
    call write_real_entry(unitno, 'xb', xb)
    call write_real_entry(unitno, 'c_k2', c_k2)
    call write_real_entry(unitno, 'c_k4', c_k4)

    call write_string_entry(unitno, 'gasmodel', gasmodel)
    call write_integer_entry(unitno, 'nchem_source', nchem_source)
    call write_integer_entry(unitno, 'nchem_rad', nchem_rad)

    call write_integer_entry(unitno, 'gcl', gcl)
    call write_integer_entry(unitno, 'cic1', cic1)
    call write_integer_entry(unitno, 'cic2', cic2)
    call write_integer_entry(unitno, 'cic3', cic3)
    call write_integer_entry(unitno, 'cic4', cic4)
    call write_integer_entry(unitno, 'cic5', cic5)

    call write_integer_entry(unitno, 'connect_point', connect_point)
    call write_string_entry(unitno, 'conpointname', conpointname)
    call write_integer_entry(unitno, 'connect_order', connect_order)
    call write_real_entry(unitno, 'dis_tao', dis_tao)

    close(unitno)

contains
    function real_to_string(value) result(out)
        real, intent(in) :: value
        character(len=64) :: out
        integer :: idx

        write(out, '(ES24.17E3)') value
        out = adjustl(out)
        do idx = 1, len(out)
            if (out(idx:idx) == 'E' .or. out(idx:idx) == 'D') out(idx:idx) = 'e'
        end do
    end function real_to_string

    subroutine write_real_entry(unit_id, name, value)
        integer, intent(in) :: unit_id
        character(len=*), intent(in) :: name
        real, intent(in) :: value
        character(len=64) :: formatted

        formatted = real_to_string(value)
        write(unit_id, '(A, '' real '', A)') trim(name), trim(formatted)
    end subroutine write_real_entry

    subroutine write_integer_entry(unit_id, name, value)
        integer, intent(in) :: unit_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: value

        write(unit_id, '(A, '' integer '', I0)') trim(name), value
    end subroutine write_integer_entry

    subroutine write_string_entry(unit_id, name, value)
        integer, intent(in) :: unit_id
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: value

        write(unit_id, '(A, '' character '', A)') trim(name), trim(value)
    end subroutine write_string_entry
end subroutine write_param_dump
!___________________________________________________________________!
subroutine grid_generation
    use global_variables,only : ndim,nblocks,nmax,gridname,mb_dim,mb_x,mb_y,mb_z
    implicit none
    integer :: i,j,k
    integer :: nb,idim,jdim,kdim
    integer :: ngridpoint
    nmax = 0
    ngridpoint = 0
    write(*,*)'��������',gridname
    open(1,file=gridname,status='unknown')
    read(1,*)nblocks !��������ֿ�����
    write(*,*)'����ֿ���',nblocks
    if ( ndim == 3 ) then
       allocate( mb_x(nblocks), mb_y(nblocks), mb_z(nblocks) )
       allocate( mb_dim(nblocks,3) )

       do nb=1,nblocks
          read(1,*)idim,jdim,kdim
          mb_dim(nb,1) = idim
          mb_dim(nb,2) = jdim
          mb_dim(nb,3) = kdim
          nmax = max(idim,jdim,kdim,nmax)
          ngridpoint = ngridpoint + idim*jdim*kdim
          write(*,*)nb,mb_dim(nb,1),mb_dim(nb,2),mb_dim(nb,3)
          allocate( mb_x(nb)%a3d(1:idim,1:jdim,1:kdim) )
          allocate( mb_y(nb)%a3d(1:idim,1:jdim,1:kdim) )
          allocate( mb_z(nb)%a3d(1:idim,1:jdim,1:kdim) )

          read(1,*) (((mb_x(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim), &
                    (((mb_y(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim), &
                    (((mb_z(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)


       enddo
    endif

    close(1)
    write(*,*)'�����ܵ���',ngridpoint


    return
end subroutine grid_generation
!___________________________________________________________________!
subroutine read_plot3d_grid
    use global_variables,only : ndim,nblocks,nmax,gridname,mb_dim,mb_x,mb_y,mb_z
    implicit none
    integer :: i,j,k
    integer :: nb,idim,jdim,kdim
    integer :: ngridpoint
    nmax = 0
    ngridpoint = 0
    write(*,*)'read mesh filename: ',gridname
    open(101,file=gridname,form='unformatted',status='old')
    read(101) nblocks 
    write(*,*)'num of blocks: ',nblocks
    if ( ndim == 3 ) then
       allocate( mb_x(nblocks), mb_y(nblocks), mb_z(nblocks) )
       allocate( mb_dim(nblocks,3) )

       do nb=1,nblocks
          read(101)idim,jdim,kdim
          mb_dim(nb,1) = idim
          mb_dim(nb,2) = jdim
          mb_dim(nb,3) = kdim
       end do

       do nb=1,nblocks
          idim = mb_dim(nb,1)
          jdim = mb_dim(nb,2)
          kdim = mb_dim(nb,3)
          nmax = max(idim,jdim,kdim,nmax)
          ngridpoint = ngridpoint + idim*jdim*kdim
          write(*,*)nb,mb_dim(nb,1),mb_dim(nb,2),mb_dim(nb,3)
          allocate( mb_x(nb)%a3d(1:idim,1:jdim,1:kdim) )
          allocate( mb_y(nb)%a3d(1:idim,1:jdim,1:kdim) )
          allocate( mb_z(nb)%a3d(1:idim,1:jdim,1:kdim) )

          read(101) (((mb_x(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim), &
                    (((mb_y(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim), &
                    (((mb_z(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)


       enddo
    endif

    close(101)
    write(*,*)'Total num of grids: ',ngridpoint


    return
end subroutine read_plot3d_grid
!_____________________________________________________________________!
subroutine read_bc
    use global_variables
    implicit none
    integer :: flow_solver_id,number_of_blocks,nb,nbt,bctype
    integer :: ndif,nr,nrmax,ntemp,js1,js2,ks1,ks2,ls1,ls2
    integer :: m,n,i,j,k,imax,jmax,kmax,tmp
    integer :: co,op(2),t_coor(2)
    integer :: s_nd,s_lr,s_fix,idelt,jdelt,kdelt
    integer :: t_nd,t_lr,t_fix,s_t_dirction(3,3)
    integer :: s_index,t_index,s_sign(3),t_sign(3),st_sign(3)
    integer,dimension(1:3) :: s_st,s_ed,t_st,t_ed
    character(len=120) :: blockname
    open(1,file=bcname,status='unknown')

    read(1,*)flow_solver_id
    read(1,*)number_of_blocks
    write(*,*)'���߽��ļ���Ϣ'
    if ( number_of_blocks /= nblocks ) then
       write(*,*)'�߽������ļ��������ļ��зֿ�����һ��!!'
       stop
    endif
    allocate( mb_bc(nblocks) )
    do nb=1,nblocks
       read(1,*)imax,jmax,kmax         !��ÿ�������ά��
       ndif = abs(imax - mb_dim(nb,1)) + abs(jmax - mb_dim(nb,2)) + abs(kmax - mb_dim(nb,3))

       if ( ndif /= 0 ) then
          write(*,*)'��',nb,'��߽������ļ�������ά���������ļ�������ά����һ��!'
          stop
       endif
       read(1,*)blockname              !��ÿ��Ŀ���
       read(1,*)nrmax                  !��ÿ�������(�߽�)��
       !write(*,*)'nrmax=',nrmax
       mb_bc(nb)%nregions = nrmax
       allocate( mb_bc(nb)%bc(nrmax) )
       do nr=1,nrmax
          read(1,*)s_st(1),s_ed(1),s_st(2),s_ed(2),s_st(3),s_ed(3),bctype
          mb_bc(nb)%bc(nr)%bctype = bctype
          mb_bc(nb)%bc(nr)%nbs = nb     !������е���࣬����������д�ϰ�
          do m=1,3
             mb_bc(nb)%bc(nr)%s_st(m) = min( abs(s_st(m)),abs(s_ed(m)) )
             mb_bc(nb)%bc(nr)%s_ed(m) = max( abs(s_st(m)),abs(s_ed(m)) )
             mb_bc(nb)%bc(nr)%s_lr3d(m) = 0
             if( s_st(m) == s_ed(m) ) then
                 s_nd = m
                 if( s_st(m) == 1 )then
                     s_lr = -1
                     mb_bc(nb)%bc(nr)%s_lr3d(m) = -1
                 else
                     s_lr =  1
                     mb_bc(nb)%bc(nr)%s_lr3d(m) =  1
                 endif
                 s_fix = s_st(m)
             else
                mb_bc(nb)%bc(nr)%s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m) + method - 1
             endif
          enddo

          mb_bc(nb)%bc(nr)%s_lr  = s_lr
          mb_bc(nb)%bc(nr)%s_nd  = s_nd
          mb_bc(nb)%bc(nr)%s_fix = s_fix

          if ( bctype < 0 ) then !�����ƴ����
			    call read_bc_connect(1,nb,nr,s_nd,s_st,s_ed)
		  endif

       enddo
    enddo

    close(1)
    write(*,*)'finished reading bc info'
    return
end subroutine read_bc
!_____________________________________________________________________!
subroutine read_bc_mbs
    use global_variables
    implicit none
    integer :: flow_solver_id,number_of_blocks,nb,nbt,bctype
    integer :: ndif,nr,nrmax,ntemp,js1,js2,ks1,ks2,ls1,ls2
    integer :: m,n,i,j,k,imax,jmax,kmax,tmp
    integer :: co,op(2),t_coor(2)
    integer :: s_nd,s_lr,s_fix,idelt,jdelt,kdelt
    integer :: t_nd,t_lr,t_fix,s_t_dirction(3,3)
    integer :: s_index,t_index,s_sign(3),t_sign(3),st_sign(3)
    integer,dimension(1:3) :: s_st,s_ed,t_st,t_ed
    character(len=120) :: blockname
    open(1,file=bcname,status='unknown')

    read(1,*)flow_solver_id
    if(flow_solver_id == 1976)read(1,*)tmp
    read(1,*)number_of_blocks
    write(*,*)'���߽��ļ���Ϣ'
    if ( number_of_blocks /= nblocks ) then
       write(*,*)'�߽������ļ��������ļ��зֿ�����һ��!!'
       stop
    endif
    allocate( mb_bc(nblocks) )
    do nb=1,nblocks
       if(flow_solver_id == 1976)read(1,*)tmp
       read(1,*)imax,jmax,kmax         !��ÿ�������ά��
       ndif = abs(imax - mb_dim(nb,1)) + abs(jmax - mb_dim(nb,2)) + abs(kmax - mb_dim(nb,3))

       if ( ndif /= 0 ) then
          write(*,*)'��',nb,'��߽������ļ�������ά���������ļ�������ά����һ��!'
          stop
       endif
       read(1,*)blockname              !��ÿ��Ŀ���
       read(1,*)nrmax                  !��ÿ�������(�߽�)��
       !write(*,*)'nrmax=',nrmax
       mb_bc(nb)%nregions = nrmax
       allocate( mb_bc(nb)%bc(nrmax) )
       do nr=1,nrmax
          read(1,*)s_st(1),s_ed(1),s_st(2),s_ed(2),s_st(3),s_ed(3),bctype
          mb_bc(nb)%bc(nr)%bctype = bctype
          mb_bc(nb)%bc(nr)%nbs = nb     !������е���࣬����������д�ϰ�
          do m=1,3
             mb_bc(nb)%bc(nr)%s_st(m) = min( abs(s_st(m)),abs(s_ed(m)) )
             mb_bc(nb)%bc(nr)%s_ed(m) = max( abs(s_st(m)),abs(s_ed(m)) )
             mb_bc(nb)%bc(nr)%s_lr3d(m) = 0
             if( s_st(m) == s_ed(m) ) then
                 s_nd = m
                 if( s_st(m) == 1 )then
                     s_lr = -1
                     mb_bc(nb)%bc(nr)%s_lr3d(m) = -1
                 else
                     s_lr =  1
                     mb_bc(nb)%bc(nr)%s_lr3d(m) =  1
                 endif
                 s_fix = s_st(m)
             else
                mb_bc(nb)%bc(nr)%s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m) + method - 1
             endif
          enddo

          mb_bc(nb)%bc(nr)%s_lr  = s_lr
          mb_bc(nb)%bc(nr)%s_nd  = s_nd
          mb_bc(nb)%bc(nr)%s_fix = s_fix

          if ( bctype < 0 ) then !�����ƴ����
			    call read_bc_connect_mbs(1,nb,nr,s_nd,s_st,s_ed,flow_solver_id)
		  endif

       enddo
    enddo

    close(1)
    write(*,*)'finished reading bc info'
    return
end subroutine read_bc_mbs
!_____________________________________________________________________!
subroutine read_bc_connect_mbs(fileid,nb,nr,s_nd,s_st,s_ed,flow_solver_id)
    use global_variables
    implicit none
    integer :: fileid,flow_solver_id
    integer :: ndif,nr,nrmax,ntemp,js1,js2,ks1,ks2,ls1,ls2
    integer :: m,n,i,j,k,imax,jmax,kmax,nbt,nb,tmp
    integer :: co,op(2),t_coor(2)
    integer :: s_nd,s_lr,s_fix,idelt,jdelt,kdelt
    integer :: t_nd,t_lr,t_fix,s_t_dirction(3,3)
    integer :: s_index,t_index,s_sign(3),t_sign(3),st_sign(3)
    integer,dimension(1:3) :: s_st,s_ed,t_st,t_ed

    if (flow_solver_id == 1976) then
       read(1,*)t_st(1),t_ed(1),t_st(2),t_ed(2),t_st(3),t_ed(3),nbt,tmp
    else
       read(1,*)t_st(1),t_ed(1),t_st(2),t_ed(2),t_st(3),t_ed(3),nbt
    end if
    do m=1,3
       mb_bc(nb)%bc(nr)%t_lr3d(m) = 0
       if( t_st(m) == t_ed(m) ) then
           t_nd = m
           if( t_st(m) == 1 )then
               t_lr = -1
               mb_bc(nb)%bc(nr)%t_lr3d(m) = -1
           else
               t_lr =  1
               mb_bc(nb)%bc(nr)%t_lr3d(m) = 1
           endif
           t_fix = t_st(m)
       endif
    enddo

    mb_bc(nb)%bc(nr)%nbt   = nbt
    mb_bc(nb)%bc(nr)%t_lr  = t_lr
    mb_bc(nb)%bc(nr)%t_nd  = t_nd
    mb_bc(nb)%bc(nr)%t_fix = t_fix

    !ȷ���Խӷ���
    do m=1,3
       do n=1,3
          s_t_dirction(m,n) = 0
       enddo
    enddo

    do m=1,3
       do n=1,3
          if ( m /= s_nd .and. n /= t_nd ) then
             js1 = s_st(m)
             if ( abs(js1) < abs(s_ed(m)) ) js1 = s_ed(m)
             js2 = t_st(n)
             if ( abs(js2) < abs(t_ed(n)) ) js2 = t_ed(n)
             if ( js1*js2 > 0 ) then
                s_t_dirction(n        , m        ) = 1
                s_t_dirction(t_nd     , s_nd     ) = 1
                s_t_dirction(6-n-t_nd , 6-m-s_nd ) = 1
                goto 10
             endif
          endif
       enddo
    enddo

10  continue

!    ȷ���ű�������1 ����-1 �������ű��Ϊ��
    do m=1,3
       s_sign(m) = 1
       s_st(m) =  abs(s_st(m) )
       s_ed(m) =  abs(s_ed(m) )
       if ( m /= s_nd ) then
          if ( s_st(m) > s_ed(m) ) s_sign(m) = -1
       endif

       t_sign(m) = 1
       t_st(m) =  abs(t_st(m) )
       t_ed(m) =  abs(t_ed(m) )
       if ( m /= t_nd ) then
          if ( t_st(m) > t_ed(m) ) t_sign(m) = -1
       endif
    enddo

    do m=1,3  !s_t_dirction(m)  �ű�ת������
       co = 0
       do n=1,3
          if(s_t_dirction(n,m)==1) co = n
       enddo
       st_sign(m) = t_sign(co)*s_sign(m)
    enddo

    js1 = mb_bc(nb)%bc(nr)%s_st(1)
    js2 = mb_bc(nb)%bc(nr)%s_ed(1)

    ks1 = mb_bc(nb)%bc(nr)%s_st(2)
    ks2 = mb_bc(nb)%bc(nr)%s_ed(2)

    ls1 = mb_bc(nb)%bc(nr)%s_st(3)
    ls2 = mb_bc(nb)%bc(nr)%s_ed(3)

    allocate( mb_bc(nb)%bc(nr)%image(js1:js2,ks1:ks2,ls1:ls2) )
    allocate( mb_bc(nb)%bc(nr)%jmage(js1:js2,ks1:ks2,ls1:ls2) )
    allocate( mb_bc(nb)%bc(nr)%kmage(js1:js2,ks1:ks2,ls1:ls2) )
    do i = s_st(1),s_ed(1),s_sign(1)
       idelt = (i-s_st(1))*st_sign(1)
       do j = s_st(2),s_ed(2),s_sign(2)
          jdelt = (j-s_st(2))*st_sign(2)
          do k = s_st(3),s_ed(3),s_sign(3)
             kdelt = (k - s_st(3))*st_sign(3)
             co    = s_t_dirction(1,1)*idelt + s_t_dirction(1,2)*jdelt + s_t_dirction(1,3)*kdelt
             mb_bc(nb)%bc(nr)%image(i,j,k) = t_st(1) + co
             co    = s_t_dirction(2,1)*idelt + s_t_dirction(2,2)*jdelt + s_t_dirction(2,3)*kdelt
             mb_bc(nb)%bc(nr)%jmage(i,j,k) = t_st(2) + co   !*st_sign(2)
             co    = s_t_dirction(3,1)*idelt + s_t_dirction(3,2)*jdelt + s_t_dirction(3,3)*kdelt
             mb_bc(nb)%bc(nr)%kmage(i,j,k) = t_st(3) + co   !*st_sign(3)
          enddo
       enddo
    enddo
    return
end subroutine read_bc_connect_mbs
!_____________________________________________________________________!
subroutine allocate_other_variable
    use global_variables
    use bl_sub_variables
    use mod_parallels,only : pnblocks,pnbindexs
    implicit none
    integer :: nb,idim,jdim,kdim,npdi,npdj,idimp1,jdimp1,kdimp1,idimp3,jdimp3,kdimp3
    integer :: nfile_par,i,j,k,left1,right1,left3,right3,pnb

    allocate( mb_vol(nblocks),mb_volt(nblocks) )
    allocate( mb_kcx(nblocks), mb_kcy(nblocks), mb_kcz(nblocks), mb_kct(nblocks) )
    allocate( mb_etx(nblocks), mb_ety(nblocks), mb_etz(nblocks), mb_ett(nblocks) )
    allocate( mb_ctx(nblocks), mb_cty(nblocks), mb_ctz(nblocks), mb_ctt(nblocks) )
    
    allocate( mb_flg(nblocks,6) )               !  lhy
    
    allocate( mb_r(nblocks), mb_u(nblocks), mb_v(nblocks), mb_w(nblocks), mb_p(nblocks) )
    allocate( mb_t(nblocks), mb_c(nblocks) )
    allocate( mb_sra(nblocks), mb_srb(nblocks), mb_src(nblocks))
    allocate( mb_srva(nblocks), mb_srvb(nblocks), mb_srvc(nblocks))
    allocate( mb_dtdt(nblocks) )
    allocate( mb_visl(nblocks) )
    allocate( mb_vist(nblocks) )  !��������ճ��ϵ��,��ṹ��ϵ,����Ҳ����
    allocate( mb_q(nblocks)   , mb_dq(nblocks)   )

    ! *** ˫ʱ�䲽���� ***
    allocate( mb_qnc(nblocks)   , mb_qmc(nblocks)   )

    allocate( mb_rhs1(nblocks) )

    !���ǻ�ѧ��Ӧ
    allocate( mb_fs(nblocks)  )  !����������ַ��̼���������
    allocate( mb_srs(nblocks) )  !���仯ѧ��ӦԴ���װ뾶

    !��������
    allocate( mb_qke(nblocks)) !, mb_dqke(nblocks) ) !�����غ������������
    allocate( mb_distance(nblocks))                !  Normal Distence from Solid Surface

    
	left1  = -1
	right1 = 1
	
	left3  = -2
	right3 = 3

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
        idim = mb_dim(nb,1)
        jdim = mb_dim(nb,2)
        kdim = mb_dim(nb,3)

        idimp1 = idim + right1
        jdimp1 = jdim + right1
        kdimp1 = kdim + right1

        idimp3 = idim + right3
        jdimp3 = jdim + right3
        kdimp3 = kdim + right3

        allocate( mb_vol(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )

        allocate( mb_kcx(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_kcy(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_kcz(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_kct(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )

        allocate( mb_etx(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_ety(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_etz(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_ett(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )

        allocate( mb_ctx(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_cty(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_ctz(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_ctt(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        
        allocate( mb_flg(nb,1)%a3d(   1:1   ,1:jdim,1:kdim) )
        allocate( mb_flg(nb,2)%a3d(idim:idim,1:jdim,1:kdim) )
        allocate( mb_flg(nb,3)%a3d(1:idim,   1:1   ,1:kdim) )
        allocate( mb_flg(nb,4)%a3d(1:idim,jdim:jdim,1:kdim) )
        allocate( mb_flg(nb,5)%a3d(1:idim,1:jdim,   1:1   ) )
        allocate( mb_flg(nb,6)%a3d(1:idim,1:jdim,kdim:kdim) )

        allocate( mb_r(nb)%a3d(left3:idimp3,left3:jdimp3,left3:kdimp3) )
        allocate( mb_u(nb)%a3d(left3:idimp3,left3:jdimp3,left3:kdimp3) )
        allocate( mb_v(nb)%a3d(left3:idimp3,left3:jdimp3,left3:kdimp3) )
        allocate( mb_w(nb)%a3d(left3:idimp3,left3:jdimp3,left3:kdimp3) )
        allocate( mb_p(nb)%a3d(left3:idimp3,left3:jdimp3,left3:kdimp3) )
        allocate( mb_t(nb)%a3d(left3:idimp3,left3:jdimp3,left3:kdimp3) )

        allocate( mb_c(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_sra(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_srb(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_src(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_srva(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_srvb(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_srvc(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )
        allocate( mb_visl(nb)%a3d(left1:idimp1,left1:jdimp1,left1:kdimp1) )

        allocate( mb_dtdt(nb)%a3d(1:idim,1:jdim,1:kdim) )
        allocate( mb_q(nb)%a4d (1:nl,left3:idimp3,left3:jdimp3,left3:kdimp3) )
        allocate( mb_dq(nb)%a4d(1:nl,left3:idimp3,left3:jdimp3,left3:kdimp3) )

       ! *** ˫ʱ�䲽���� ***
       if (ndualtst > 0) then
          allocate( mb_qnc(nb)%a4d(1:nl,left3:idimp3,left3:jdimp3,left3:kdimp3) )
          allocate( mb_qmc(nb)%a4d(1:nl,left3:idimp3,left3:jdimp3,left3:kdimp3) )
       end if
        
        allocate( mb_vist(nb)%a3d( left3:idimp3,left3:jdimp3,left3:kdimp3 ) ) !��������ճ��ϵ��
        do k=left3,kdimp3
        do j=left3,jdimp3
        do i=left3,idimp3
           mb_vist(nb)%a3d( i,j,k ) = 0.0
        enddo 
        enddo
        enddo
        if ( nlamtur >= 0 ) then       
           allocate( mb_qke(nb)%a4d (left3:idimp3,left3:jdimp3,left3:kdimp3,nlamtur) )
        end if
      
		if ( (nvis>0) .and. (nlamtur>=0) )then
           allocate( mb_distance(nb)%a3D (1:idim,1:jdim,1:kdim))
		endif
    enddo

      if( (nvis>0) .and. (nlamtur>=0) .and. nameturb == "BL" )then
         call allocate_bl_sub_variables
      endif
	
	call fill_bc_flag  !!lhy

    return
end subroutine allocate_other_variable
!_____________________________________________________________________!
subroutine fill_bc_flag
    use global_variables,only : nblocks,mb_bc,mb_flg
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
    integer :: nb,pnb,ib,nr,ibctype,m,i,j,k
    integer :: nregions,ibeg(3),iend(3),idir,inrout
    
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
       nregions = mb_bc(nb)%nregions

       do ib=1,nregions
          nr = mb_bc(nb)%bcindexs(ib)
          ibctype = mb_bc(nb)%bc(nr)%bctype
          ibeg = mb_bc(nb)%bc(nr)%s_st
          iend = mb_bc(nb)%bc(nr)%s_ed
          idir = mb_bc(nb)%bc(nr)%s_nd
          inrout = mb_bc(nb)%bc(nr)%s_lr
          
          m = 2*idir + (inrout-1)/2
          
          do k=ibeg(3),iend(3)
          do j=ibeg(2),iend(2)
          do i=ibeg(1),iend(1) 
             mb_flg(nb,m)%a3d(i,j,k) = nr
          end do
          end do
          end do
       end do
    end do
    
end subroutine fill_bc_flag
!_____________________________________________________________________!
subroutine read_bc_par
    use global_variables,&
    only : nblocks,mb_bc,twall,bc_par_len,bc_par,roo,poo,rref,vref,pref
    implicit none
    integer :: nb,nr,bctype,m,nrec,nrecord,bc_disf
    integer :: s_st(3),s_ed(3),i,j,k
	real    :: abc
    nrec = 0
    bc_disf = 0
    do nb=1,nblocks
       do nr=1,mb_bc(nb)%nregions
          bctype = mb_bc(nb)%bc(nr)%bctype
          mb_bc(nb)%bc(nr)%bc_par_len = 0
          if ( bctype == 20 .or. bctype == 22 ) then
             mb_bc(nb)%bc(nr)%bc_par_len = 1
             allocate( mb_bc(nb)%bc(nr)%bc_par( mb_bc(nb)%bc(nr)%bc_par_len ) )
             mb_bc(nb)%bc(nr)%bc_par(1)= twall
          endif
       enddo
    enddo
    if ( nrec > 0 ) then
       open(1,file= 'grid/bc_par.dat' ,form='formatted',status='old')
       read(1,*)nrecord
       do nrec=1,nrecord
          read(1,*)nb, nr, bctype, bc_par_len
		  if(mb_bc(nb)%bc(nr)%bctype /= bctype ) then
		     write(6,*)'bc_par.dat ��,��',nb,'��ĵ�',nr,'���ڴ��ڴ��� ������'
		  endif
          allocate( bc_par( bc_par_len ) )
          read(1,*)(bc_par(m),m=1,bc_par_len)
          write(*,*)'bc_par=',bc_par
          if ( bctype == 51 .or. bctype == 53 ) then     !��ȫ�������ȷֲ�
             bc_par(1) = bc_par(1)*roo/rref
             bc_par(2) = bc_par(2)/vref
             bc_par(3) = bc_par(3)/vref
             bc_par(4) = bc_par(4)/vref
             bc_par(5) = bc_par(5)*poo/pref
          endif

          write(*,*)'bc_par=',bc_par
          mb_bc(nb)%bc(nr)%bc_par_len =  bc_par_len
          mb_bc(nb)%bc(nr)%bc_par     => bc_par
          write(*,*) mb_bc(nb)%bc(nr)%bc_par,' nb=',nb,' nr=',nr
       enddo
       close(1)
    endif

    if ( bc_disf > 0 ) then
       open(1,file= 'bc_par1.dat' ,form='formatted',status='old')
       read(1,*)nrecord
       do nrec=1,nrecord
          read(1,*)nb, nr, bctype, bc_par_len
          read(1,*)(s_st(i),s_ed(i),i=1,3)
          !�жϵ�����߽��ļ��������Ƿ�һ��
          do m=1,3
             i = s_st(m) - mb_bc(nb)%bc(nr)%s_st(m)
             j = s_ed(m) - mb_bc(nb)%bc(nr)%s_ed(m)
             if ( i<0 .or. j<0 ) then
                write(*,*)'���棺 �߽�����ļ�bc_par1.dat �ĵ�',nb
                write(*,*)'�����ĵ�',nr,'���߽����������'
             endif
          enddo
          mb_bc(nb)%bc(nr)%bc_par_len =  bc_par_len

          allocate( mb_bc(nb)%bc(nr)%bc_par1( bc_par_len, &
                    s_st(1):s_ed(1), s_st(2):s_ed(2), s_st(3):s_ed(3)))
          do i=s_st(1),s_ed(1)
             do j=s_st(2),s_ed(2)
                do k=s_st(3),s_ed(3)
                   read(1,*)( mb_bc(nb)%bc(nr)%bc_par1(m,i,j,k), m = 1,bc_par_len )
                enddo
             enddo
          enddo
          if ( bctype == 52 .or. bctype == 54 ) then
             do i=s_st(1),s_ed(1)
             do j=s_st(2),s_ed(2)
                do k=s_st(3),s_ed(3)
				   mb_bc(nb)%bc(nr)%bc_par1(1,i,j,k) = mb_bc(nb)%bc(nr)%bc_par1(1,i,j,k)*roo/rref
				   mb_bc(nb)%bc(nr)%bc_par1(2,i,j,k) = mb_bc(nb)%bc(nr)%bc_par1(2,i,j,k)/vref
				   mb_bc(nb)%bc(nr)%bc_par1(3,i,j,k) = mb_bc(nb)%bc(nr)%bc_par1(3,i,j,k)/vref
				   mb_bc(nb)%bc(nr)%bc_par1(4,i,j,k) = mb_bc(nb)%bc(nr)%bc_par1(4,i,j,k)/vref
				   mb_bc(nb)%bc(nr)%bc_par1(5,i,j,k) = mb_bc(nb)%bc(nr)%bc_par1(5,i,j,k)*poo/pref
                enddo
             enddo
             enddo
          endif

       enddo
       close(1)

    endif
    return
end subroutine read_bc_par
!_____________________________________________________________________!
subroutine init_n
    use global_variables,&
        only : nblocks,nchem,mb_par_content,mb_control,ntmodel,ns,ne,nc,nl,nm
    implicit none
    integer :: nb,nfile_par

    do nb = 1, nblocks
       nfile_par = mb_par_content(nb)
       nchem = mb_control( nfile_par )%nchem
!          mb_control( nfile_par )%nl = nm
    enddo

    return
end subroutine init_n
!_____________________________________________________________________!
subroutine init_inflow
    use global_variables
    use resdual_module
#ifdef PARALLEL
    use mod_parallels,only : master,myid
#endif
    implicit none
    integer :: is
    real :: mavd1,cp,cv,ccoo,ae,aaaa,ccon
    real :: cps(ns),hs(ns)

!    call rev_comp_par_nochem( 1 )
!____________
!   �󹥽ǡ��໬�Ǽ������ٶ�
    attack = attack * pai / 180.0
    sideslip = sideslip * pai / 180.0

    uoo = cos(attack)*cos(sideslip)
    voo = sin(attack)*cos(sideslip)
    woo = sin(sideslip)

!   ����������
    roo = 1.0
    too = 1.0
!____________
!   ����reynolds=������������ŵ��/��������x����� added by clz 2006.5.24
    reynolds = (reynolds/lfref)*lref

    do is=1,ns
       ws(is) = ws(is) * 1.0e-3        !��Ԫ������(������)
    enddo

    do is=1,ns
       ws1(is) = 1.0/ws(is)            !��Ԫ�������ĵ���(������)
    enddo

    mref1 = 0.0                        !ƽ���������ĵ���(������)
    do is=1,ns
       mref1 = mref1 + cn_init(is+ne) * ws1(is)
    enddo

    mref = 1.0/mref1                   !ƽ��������(������)

    do is=1,ns
      ms(is) = ws(is) * mref1          !����ַ����������ٻ�
    enddo

    do is=1,ns
      ms1(is) = 1.0/ms(is)             !����ַ������ĵ���(������)
    enddo
!___________
    !��ο�ѹ���Ͳο��ܶ�
    if( height < 0.0000 ) then
       if(pref > 0.0 .and. tref > 0.0 .and. rref > 0.0 )then
          write(*,*)'ѹ�����¶Ⱥ��ܶȲ���ͬʱ������'
          stop
       endif
       if ( pref == 0.0 .and. rref == 0.0 ) then
          write(*,*)'ѹ�����ܶȲ���ͬʱ������'
       endif
       if    ( pref > 0.0 .and. tref > 0.0 .and. rref <= 0.0 ) then
          rref = pref/(tref*rjmk*mref1)
       elseif( rref > 0.0 .and. tref > 0.0 .and. pref <= 0.0 ) then
          pref = rref*tref*rjmk*mref1
       elseif( pref > 0.0 .and. rref > 0.0 .and. tref <= 0.0 ) then
          tref = pref/(rref*rjmk*mref1)
       else
          write(*,*)'����ѹ�����¶Ⱥ��ܶȲ�����'
          stop
       endif
       ccon   = (1.+110.4/273.15)/(1.+110.4/tref)
       visloo = sqrt(tref/273.15)*ccon*1.715e-5
       ccoo   = sqrt(gama*pref/rref)
    else
       call air(height,tref,pref,rref,ccoo)
       rref   = pref/(tref*rjmk*mref1)
       ccon   = (1.+110.4/273.15)/(1.+110.4/tref)
       visloo = sqrt(tref/273.15)*ccon*1.715e-5
    endif

#ifdef PARALLEL
    if (myid == master) then
#endif

    write(*,*)'      高度         温度         压力         密度         音速'
    write(*,10)height,tref,pref,rref,ccoo
10  format(1x,5e13.4)

#ifdef PARALLEL
    end if
#endif

!___________

    coo = sqrt(gama*rjmk*mref1*tref)  !(������)
       aaaa = rjmk
!	   moo = vref/ccoo
       vref = ccoo * moo                  !����ο��ٶ�(������)
       rvl   = rref * vref / lref
       rvl1  = 1.0/rvl
       vl    = vref / lref
       vl1   = 1.0/vl
       beta1 = rjmk*tref/(vref*vref*mref)
       if ( height >= 0.0 ) reynolds = rref*vref*lref/visloo

       poo = 1.0/(gama*moo*moo)
       coo = sqrt(gama*poo/roo)
       eoo = poo/(gama-1.0) + 0.5*roo*( uoo*uoo + voo*voo + woo*woo )
       hoo = (eoo + poo)/roo
!	   reynolds = rref*vref*lref/visloo  !*tgh. �ر�ע��
       !*tgh. �ر�ע����ŵ�����ر�ע��re

#ifdef PARALLEL
    if (myid == master) then
#endif

    write(*,*)'     雷诺数    参考粘性系数    马赫数      参考速度'
    write(*,10)reynolds,visloo,moo,vref
    write(*,*)'     特征长度  网格长度比'
    write(*,10)lfref,lref

#ifdef PARALLEL
    end if
#endif

    !����������
    q1oo = roo
    q2oo = roo * uoo
    q3oo = roo * voo
    q4oo = roo * woo
    q5oo = eoo
    !��פ��ѹ�����ܶ�
    call pr_stag

    cq_lam = 1.0/((gama-1.0)*moo*moo*prl)
    cq_tur = 1.0/((gama-1.0)*moo*moo*prt)
    visc   = 110.4/tref
    !ע���ʼ��re
    re     = 1.0/reynolds

    allocate( nres(nblocks,4), res(nblocks), resmax(nblocks) )

	call inflow_turbulence  ! added by clz 2006.5.28
	
    return
end subroutine init_inflow
!_____________________________________________________________________!
subroutine pr_stag
    use global_const,only : moo,gama,pst,rst,poo,roo,too
    implicit none
    real :: m2,gam1,ogm1,gam2,arg
    
    m2 = moo * moo
    
    gam1 = gama - 1.0
    ogm1 = 1.0/gam1
    gam2 = gama + 1.0
    if ( moo >= 1.0 ) then
       arg = 1.0 + 0.5*gam1*m2
       pst = poo * (0.5*gam2*m2)**(gama*ogm1) / (2.0*gama*m2/gam2 - gam1/gam2)**ogm1
       rst = gama*m2*pst/(too*arg)
    else
       arg = 1.0 + 0.5*gam1*m2
       pst = poo * arg**(gama*ogm1)
       rst = roo * arg**ogm1
    endif

    return
end subroutine pr_stag
!_____________________________________________________________________!
subroutine initialization
    use global_variables
    implicit none
    real :: ratio_mref
    integer :: nfile_par,is,nb,i,j,k

 !   mb_control( 1 )%ratio_mref = 1.0


	  if( nvis > 0 .and. nlamtur >= 0  ) then
		  call init_turbulence    ! Initialization of Turblence Models
	  endif                     ! Modified by clz

    if ( nstart == 0 ) then
       call init_nstart0
    else
       call init_nstart1
    endif

    return
end subroutine initialization
!_____________________________________________________________________!
subroutine set_bc2_to_p
    use global_variables
    implicit none
    integer :: i,j,k,m

    do k=0,nk+1
       do j=0,nj+1
          r(  -1,j,k) = 10.0
          r(ni+1,j,k) = 10.0
       enddo
    enddo

    do k=0,nk+1
       do i=0,ni+1
          r(i,  -1,k) = 10.0
          r(i,nj+1,k) = 10.0
       enddo
    enddo

    do j=0,nj+1
       do i=0,ni+1
          r(i,j,  -1) = 10.0
          r(i,j,nk+1) = 10.0
       enddo
    enddo

    return
end subroutine set_bc2_to_p
!_____________________________________________________________________!
subroutine set_dq_to_0
    use global_variables
    implicit none
    integer :: i,j,k,m

    do k=-1,nk+1
       do j=-1,nj+1
          do i=-1,ni+1
             do m=1,nl
                dq(m,i,j,k) = 0.0
             enddo
          enddo
       enddo
    enddo

    return
end subroutine set_dq_to_0
!_____________________________________________________________________!
subroutine set_rad_to_0
    use global_variables
    implicit none
    integer :: nb
    call set_spectinv_to_0
    call set_spectvis_to_0

    return
end subroutine set_rad_to_0
!_____________________________________________________________________!
subroutine set_spectinv_to_0
    use global_variables
    implicit none
    integer :: i,j,k
    do k=-1,nk+1
       do j=-1,nj+1
          do i=-1,ni+1
             sra(i,j,k) = 0.0
             srb(i,j,k) = 0.0
             src(i,j,k) = 0.0
          enddo
       enddo
    enddo

    return
end subroutine set_spectinv_to_0
!_____________________________________________________________________!
subroutine set_spectvis_to_0
    use global_variables
    integer :: i,j,k
    do k=-1,nk+1
       do j=-1,nj+1
          do i=-1,ni+1
             srva(i,j,k) = 0.0
             srvb(i,j,k) = 0.0
             srvc(i,j,k) = 0.0
          enddo
       enddo
    enddo

    return
end subroutine set_spectvis_to_0
!_____________________________________________________________________!
subroutine init_nstart0
    use global_variables
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs,myid
#endif
    implicit none
    integer :: nb,i,j,k,m,pnb
    real :: prim(nl),q_q(nl)

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
!       call recast_method(nb)
       call recast_grid(nb)
       call recast_field(nb)
       !�Ƚ���������ı߽��Ϊ��������ֵ
       do k=1,nk
       do j=1,nj
       do i=1,ni
          r(i,j,k) = roo
          u(i,j,k) = uoo
          v(i,j,k) = voo
          w(i,j,k) = woo
          p(i,j,k) = poo
          t(i,j,k) = 1.0
          
          if (ndualtst > 0) then
            prim(1) = r(i,j,k)
            prim(2) = u(i,j,k)
            prim(3) = v(i,j,k)
            prim(4) = w(i,j,k)
            prim(5) = p(i,j,k)
            call prim_to_q(prim,q_q,gama)
            q(1,i,j,k) = q_q(1)
            q(2,i,j,k) = q_q(2)
            q(3,i,j,k) = q_q(3)
            q(4,i,j,k) = q_q(4)
            q(5,i,j,k) = q_q(5)
            do m=1,nl
               mb_qmc(nb)%a4d(m,i,j,k) = q(m,i,j,k)
               mb_qnc(nb)%a4d(m,i,j,k) = q(m,i,j,k)
            end do
          end if

       enddo
       enddo
       enddo
    enddo

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
!       call recast_method(nb)
       call recast_grid(nb)
       call recast_field(nb)
       !call init_primitive(nb)
    enddo

    return
end subroutine init_nstart0
!_____________________________________________________________________!
subroutine init_nstart1
    use global_variables
#ifdef PARALLEL
    use mod_parallels,only : read_flow_file_parallel,read_flow_dual_parallel
#endif
    implicit none
    integer :: nb,i,j,k,is

#ifdef PARALLEL
    call read_flow_file_parallel
    
    if (ndualtst > 0) call read_flow_dual_parallel
#else
    call read_flow_file
#endif

    return
end subroutine init_nstart1
!_____________________________________________________________________!
subroutine init_primitive(nb)
    use global_variables
    implicit none
    integer :: nb,i,j,k,is,nfile_par,nmeth,npdi,npdj

    npdi = -2
    npdj = 3
    nfile_par = mb_par_content(nb)
    nmeth = 1 - method
    !�̱ڰ�ţ����������������
    call newtonflow(nb)

    !��ȫ�����߽����������
    !call knownflow_bc(nb)
    !call boundary(nb)
    !call inner_bc(nb)
     
    !��ȫ�������ݱ߽��ֵ
    call infinite_interplat1(r,ni,nj,nk,npdi,npdj,method)
    call infinite_interplat1(u,ni,nj,nk,npdi,npdj,method)
    call infinite_interplat1(v,ni,nj,nk,npdi,npdj,method)
    call infinite_interplat1(w,ni,nj,nk,npdi,npdj,method)
    call infinite_interplat1(p,ni,nj,nk,npdi,npdj,method)
    if ( nchem == 1 ) call infinite_interplat2(fs,ns,ni,nj,nk,npdi,npdj,method)
    
    return
end subroutine init_primitive
!_____________________________________________________________________!
subroutine newtonflow(nb)
    use global_variables
    implicit none
    integer :: i,j,k
    integer :: nb,nr,nrmax,need_process
    integer :: bctype,m,is,js,ks,s_nd,s_lr,s_fix,t_n0
    integer :: s_st(3),s_ed(3),cmethod
    integer :: s0(3),t0(3),ijk_fix1,ijk_fix2,ijk_fix
    real    :: pn,rn,cst,nx,ny,nz,cgm

    call pr_stag
    nrmax = mb_bc(nb)%nregions     !���鹲��nrmax���߽���Ҫ����
    do nr=1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if ( bctype == 2 .or. bctype/10 == 2 ) then  !�̱ڱ߽磬��Ӧ������
          !�Թ̱ڽ���ţ��������
          do m=1,3
             s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
             s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
          enddo
          s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
          s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
          s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
          t_n0 = s_fix                             !t_n0Ϊ���굼�����ڵ���
          !�Թ̱ڵ�ѹ�����ܶȡ��ٶȸ�ֵ
          do is=s_st(1),s_ed(1)
             do js=s_st(2),s_ed(2)
                do ks=s_st(3),s_ed(3)
                   if ( s_nd == 1 ) then
                      nx = mb_kcx(nb)%a3d(is,js,ks)
                      ny = mb_kcy(nb)%a3d(is,js,ks)
                      nz = mb_kcz(nb)%a3d(is,js,ks)
                      ijk_fix2 = (1 + ni )/2
                   elseif ( s_nd == 2 ) then
                      nx = mb_etx(nb)%a3d(is,js,ks)
                      ny = mb_ety(nb)%a3d(is,js,ks)
                      nz = mb_etz(nb)%a3d(is,js,ks)
                      ijk_fix2 = (1 + nj )/2
                   else
                      nx = mb_ctx(nb)%a3d(is,js,ks)
                      ny = mb_cty(nb)%a3d(is,js,ks)
                      nz = mb_ctz(nb)%a3d(is,js,ks)
                      ijk_fix2 = (1 + nk )/2
                   endif
                   cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
                   nx  = (nx*uoo + ny*voo + nz*woo )/cgm
                   cst = nx * nx
                   pn = cst * pst + ( 1.0 - cst ) * poo
                   rn = rst * ( pn/pst )**(1.0/gama)

                   do ijk_fix = 0,ijk_fix2
                      t0(1) = is - ijk_fix * mb_bc(nb)%bc(nr)%s_lr3d(1)
                      t0(2) = js - ijk_fix * mb_bc(nb)%bc(nr)%s_lr3d(2)
                      t0(3) = ks - ijk_fix * mb_bc(nb)%bc(nr)%s_lr3d(3)
                      mb_r(nb)%a3d(t0(1),t0(2),t0(3)) = rn
                      mb_u(nb)%a3d(t0(1),t0(2),t0(3)) = 0.0
                      mb_v(nb)%a3d(t0(1),t0(2),t0(3)) = 0.0
                      mb_w(nb)%a3d(t0(1),t0(2),t0(3)) = 0.0
                      mb_p(nb)%a3d(t0(1),t0(2),t0(3)) = pn

!                      mb_r(nb)%a3d(t0(1),t0(2),t0(3)) =  roo
!                      mb_u(nb)%a3d(t0(1),t0(2),t0(3)) =  uoo
!                      mb_v(nb)%a3d(t0(1),t0(2),t0(3)) =  voo
!                      mb_w(nb)%a3d(t0(1),t0(2),t0(3)) =  woo
!                      mb_p(nb)%a3d(t0(1),t0(2),t0(3)) =  poo
                      if ( nchem == 1 ) then
                         do m=1,ns
                           fs(m,t0(1),t0(2),t0(3)) = cn_init(m)
                         enddo
                      endif
                   enddo
                enddo
             enddo
          enddo
       endif
    enddo

    return
end subroutine newtonflow
!_____________________________________________________________________!
subroutine inner_bc(nb)
    use global_variables
    implicit none
    integer :: i,j,k,ijk_fix,ijk_fix1,ijk_fix2
    integer :: nb,nr,nrmax,need_process,nmeth
    integer :: bctype,m,is,js,ks,s_nd,s_lr,s_fix,t_n0
    integer :: s_st(3),s_ed(3)
    integer :: s0(3),t0(3)
    real    :: pn,rn,cst,nx,ny,nz,cgm
    nrmax = mb_bc(nb)%nregions                  !���鹲��nrmax���߽���Ҫ����
    do nr=1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if ( bctype < 0 ) then                   !�Խӱ߽�
          do m=1,3
             s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
             s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
          enddo
          do is=s_st(1),s_ed(1)
             do js=s_st(2),s_ed(2)
                do ks=s_st(3),s_ed(3)
                   mb_r(nb)%a3d(is,js,ks) =  roo
                   mb_u(nb)%a3d(is,js,ks) =  uoo
                   mb_v(nb)%a3d(is,js,ks) =  voo
                   mb_w(nb)%a3d(is,js,ks) =  woo
                   mb_p(nb)%a3d(is,js,ks) =  poo
                   if ( nchem == 1 ) then
                      do m=1,ns
                         fs(m,is,js,ks) = cn_init(m)
                      enddo
                   endif
                enddo
             enddo
           enddo
        endif
    enddo
    return
end subroutine inner_bc
!_____________________________________________________________________!
subroutine knownflow_bc(nb)
    use global_variables
    implicit none
    integer :: i,j,k,ijk_fix,ijk_fix1,ijk_fix2,nmeth
    integer :: nb,nr,nrmax,need_process
    integer :: bctype,m,is,js,ks,s_nd,s_lr,s_fix,t_n0
    integer :: s_st(3),s_ed(3)
    integer :: s0(3),t0(3)
    real :: pn,rn,cst,nx,ny,nz,cgm
    nrmax = mb_bc(nb)%nregions     !���鹲��nrmax���߽���Ҫ����
    do nr=1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if ( bctype == 51 ) then
              !�����߽�����
          do m=1,3
             s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
             s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
          enddo
          s_nd  = mb_bc(nb)%bc(nr)%s_nd
          s_lr  = mb_bc(nb)%bc(nr)%s_lr
          s_fix = mb_bc(nb)%bc(nr)%s_fix
          t_n0  = s_fix
          if ( s_nd == 1 ) then
             ijk_fix  = (1 + ni )/2 - s_fix
             ijk_fix1 = s_fix + (1 + s_lr )*ijk_fix/2
             ijk_fix2 = s_fix + (1 - s_lr )*ijk_fix/2
          elseif ( s_nd == 2 ) then
             ijk_fix  = (1 + nj )/2 - s_fix
             ijk_fix1 = s_fix + (1 + s_lr )*ijk_fix/2
             ijk_fix2 = s_fix + (1 - s_lr )*ijk_fix/2
          else
             ijk_fix  = (1 + nk )/2 - s_fix
             ijk_fix1 = s_fix + (1 + s_lr )*ijk_fix/2
             ijk_fix2 = s_fix + (1 - s_lr )*ijk_fix/2
          endif

          do is=s_st(1),s_ed(1)
             do js=s_st(2),s_ed(2)
                do ks=s_st(3),s_ed(3)
                   do ijk_fix = ijk_fix1,ijk_fix2
                      t0(1) = is - ijk_fix * mb_bc(nb)%bc(nr)%s_lr3d(1)
                      t0(2) = js - ijk_fix * mb_bc(nb)%bc(nr)%s_lr3d(2)
                      t0(3) = ks - ijk_fix * mb_bc(nb)%bc(nr)%s_lr3d(3)
                      mb_r(nb)%a3d(t0(1),t0(2),t0(3)) = bc_par(1)
                      mb_u(nb)%a3d(t0(1),t0(2),t0(3)) = bc_par(2)
                      mb_v(nb)%a3d(t0(1),t0(2),t0(3)) = bc_par(3)
                      mb_w(nb)%a3d(t0(1),t0(2),t0(3)) = bc_par(4)
                      mb_p(nb)%a3d(t0(1),t0(2),t0(3)) = bc_par(5)
                      do m=1,nl-5
                         fs(m,t0(1),t0(2),t0(3)) = bc_par(m+5)
                      enddo
                   enddo
                enddo
               enddo
          enddo
       endif
    enddo
    return
end subroutine knownflow_bc
!_____________________________________________________________________!
subroutine read_flow_file
    use global_variables,&
    only : nblocks,ni,nj,nk,ns,nstart,nchem,nout,fs,wholetime,flowname, &
	       r,u,v,w,p,t,nlamtur,ndualtst,nl,gama,ttdts,q,mb_qnc,mb_qmc
    implicit none
    integer :: nb,i,j,k,m
    real :: tmp,prim(nl),q_q(nl)

    open(1,file=flowname,status='unknown',form='unformatted')
    read(1)nout,wholetime,tmp,tmp,tmp,tmp,tmp,tmp
    do nb=1,nblocks
!       call recast_method(nb)
       call recast_grid(nb)
       call recast_field(nb)

       read(1)(((r(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k),t(i,j,k), &
                                          i=1,ni),j=1,nj),k=1,nk)

       if (ndualtst > 0) then
          do k=1,nk
          do j=1,nj
          do i=1,ni
             prim(1) = r(i,j,k)
             prim(2) = u(i,j,k)
             prim(3) = v(i,j,k)
             prim(4) = w(i,j,k)
             prim(5) = p(i,j,k)
             call prim_to_q(prim,q_q,gama)
             q(1,i,j,k) = q_q(1)
             q(2,i,j,k) = q_q(2)
             q(3,i,j,k) = q_q(3)
             q(4,i,j,k) = q_q(4)
             q(5,i,j,k) = q_q(5)
             
             do m=1,nl
                mb_qnc(nb)%a4d(m,i,j,k) = q(m,i,j,k)
             end do
          
          enddo
          enddo
          enddo
       end if
       
    enddo
    close(1)
    
    if (ndualtst > 0) then
       if (ndualtst == 2) then
         open(1,file=trim(flowname)//".dts",status='old',form='unformatted')
         read(1)ttdts
         do nb=1,nblocks
           call recast_grid(nb)
           call recast_field(nb)
           read(1)((((mb_qmc(nb)%a4d(m,i,j,k),m=1,nl),i=1,ni),j=1,nj),k=1,nk)
         enddo
         close(1)
       else
         do nb=1,nblocks
           call recast_grid(nb)
           call recast_field(nb)
           do k=1,nk
           do j=1,nj
           do i=1,ni
              do m=1,nl
                 mb_qmc(nb)%a4d(m,i,j,k) = mb_qnc(nb)%a4d(m,i,j,k)
              end do
           enddo
           enddo
           enddo
         enddo
       end if
    end if
        
    return
end subroutine read_flow_file
!_____________________________________________________________________!
subroutine set_grid_derivative
    use global_variables,only : nblocks
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
    integer :: nb,pnb

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
       call recast_grid(nb)
       call grid_derivative_nb(nb)
    enddo

    return
end subroutine set_grid_derivative
!_____________________________________________________________________!
subroutine grid_derivative_nb(nb)
    use global_variables, &
	     only : method,ndim,ni,nj,nk,kct,ett,ctt,nscheme,gcl
    implicit none
    integer :: nb,i,j,k

    if ( ndim == 3 ) then
        call GRID_DERIVATIVE_gcl(nb)
    end if

    do k=1,nk
    do j=1,nj
    do i=1,ni
       kct(i,j,k) = 0.0
       ett(i,j,k) = 0.0
       ctt(i,j,k) = 0.0
    enddo
    enddo
    enddo

    return
end subroutine grid_derivative_nb
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------!
!   subroutine connect_pre_new                                                         !
!	���ܣ�ȷ���Խ����ϵ��λ����Ϣ���������ļ�������鹫�У���ÿ�����ϵı�ţ�         !
!         ����Ԥ�������֣��ڵõ����������ͱ߽������ļ���͸�������ֵ֮ǰ���ô�ģ�顣 !
!	���룺�߽���������Ϣ��nblocks,mb_bc��                                              !
!	������Խ����ϵ��λ����Ϣ��nptt,cdate,single,singles                              !
!	��1��nptt �Խ����ϵ㣻��2��cdateΪ�ṹ������npp��Nbijk(:,:)����Ԫ��                !
!	��3��npp��ʾ���жԽ����ϵ�����������NbijkΪλ����Ϣ                              !
!	��4��Nbijk(n,1)��Nb��Nbijk(n,2)��i��Nbijk(n,3)��i��Nbijk(n,4)��k��nΪ1��npp        !
!	��ƣ�ëö��                                                                       !
!	���ԣ�Ϳ����                ������                                                !
!   ʱ�䣺2009��2��                                                                    !
!--------------------------------------------------------------------------------------!
subroutine connect_pre_new
    use global_variables,only: nblocks,mb_bc,cdate,nptt,single,singles,connect_point,conpointname
#ifdef PARALLEL
    use mod_parallels,only : master,myid
#endif
	implicit none
	integer :: nb,nrmax,nr,bctype,npti,i,j,k,npp,ii,jj,kk,nbt,it,jt,kt,npp0
	integer,dimension(1:3) :: s_st,s_ed
	integer,pointer :: Nbijk (:,:)   !��ʱ���飬��Խ����������ռ�һ�����Ӧ�ļ�������ϵ�µĵ�������Ϣ
    integer,pointer :: Mijk (:,:,:)  !��ʱ���飬��Խ����������ռ����е��Ӧ�ļ�������ϵ�µĵ�������Ϣ
	integer :: m,nbb !,IBC_TEM

!    WRITE(*,*)'�����룺1---��ȡ�Խ����ã�0---�������Խ�����'
!	READ(*,*)IBC_TEM

	IF(connect_point ==1)THEN
	   write(*,*)'��ȡ�Խӱ߽�������������Ա�ƽ��ʱʹ��'
	   OPEN(2,FILE=conpointname,FORM="unformatted")
		   read(2)NPTT
		   allocate( cdate(nptt) )

		   do i=1,nptt
              READ(2)NPP
		      cdate(i)%npp=npp
			  allocate( cdate(i)%nbijk(npp,4) )

		      DO J=1,NPP
		         READ(2)(cdate(i)%Nbijk (j,k),k=1,4)
		      ENDDO
		ENDDO
	   CLOSE(2)

	  RETURN
	ENDIF

#ifdef PARALLEL
    if (myid == master) then
#endif

    write(*,*)'���ڼ���Խӱ߽�������������Ա�ƽ��ʱʹ��'

#ifdef PARALLEL
    end if
#endif


!   ����ȷ��nptt
    nptt=0
    do nb=1,nblocks
	   nrmax = mb_bc(nb)%nregions
	   do nr = 1,nrmax
	      bctype = mb_bc(nb)%bc(nr)%bctype
		  if( bctype < 0 ) then               !�Խӱ߽��ϸ����
			  i=mb_bc(nb)%bc(nr)%s_st(1) - mb_bc(nb)%bc(nr)%s_ed(1)
			  j=mb_bc(nb)%bc(nr)%s_st(2) - mb_bc(nb)%bc(nr)%s_ed(2)
			  k=mb_bc(nb)%bc(nr)%s_st(3) - mb_bc(nb)%bc(nr)%s_ed(3)
			  nptt = nptt + ( ABS(i)+1 ) * ( ABS(j)+1 ) * ( ABS(k)+1 )
		  endif
	   enddo
    enddo

	ALLOCATE (cdate(nptt/2+1))   !�Խ����ϵĵ㣬�����������������,��ˣ��������ռ俴��ֻ��һ�롣
	ALLOCATE (Mijk(Nptt,2,4))
    ALLOCATE (Nbijk(Nptt,4))

    Nptt = 0
 !  Mijk��¼���жԽ����ϵ�Ŀ�ź�����ָ�� i,j,k
	do nb=1,nblocks
	    nrmax = mb_bc(nb)%nregions
	    do nr = 1,nrmax
		    bctype = mb_bc(nb)%bc(nr)%bctype
			if( bctype < 0 ) then               !�Խӱ߽�
			    do m=1,3
                   s_st(m)  = mb_bc(nb)%bc(nr)%s_st(m)   !��ʼ������(�ɶ�����)
                   s_ed(m)  = mb_bc(nb)%bc(nr)%s_ed(m)   !��ֹ������(�ɶ�����)
                enddo

			    do i= s_st(1),s_ed(1)
                    do j= s_st(2),s_ed(2)
                        do k=s_st(3),s_ed(3)
					        Nptt= Nptt+1
                            Mijk (Nptt,1,1) = nb
                            Mijk (Nptt,1,2) = i
                            Mijk (Nptt,1,3) = j
                            Mijk (Nptt,1,4) = k

                            Mijk (Nptt,2,1) = mb_bc(nb)%bc(nr)%nbt
                            Mijk (Nptt,2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                            Mijk (Nptt,2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                            Mijk (Nptt,2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)
					    enddo
				    enddo
			    enddo
			endif
		enddo
	enddo

!  ������3����������ϵĶԽ���ĵ���кϲ����޳���׼ȷ�õ��������ռ俴�ĶԽ����ϵ��ܵ�����λ����Ϣ
	npti=0
	do i=1,nptt
	   if(Mijk(i,1,1)>0) then
		  do j=1,4
		     nbijk(1,j)=mijk(i,1,j)
			 nbijk(2,j)=mijk(i,2,j)
		  enddo
		  npp=2           !nppΪ���иõ���������
	!======================================================================================
	!  Ϊ�˲�©������4��ʱ���������飬�Ҳ����ø�����겹���жϣ���ˣ������һ��ѭ��kt  !
	!======================================================================================
		  do kt=1,nptt
		     npp0=npp
		     do j=i+1,nptt   !�����ж�������ʶ���ظ���û�У�ȫ�ظ����޳��������ظ���
			                 !��û���ظ��ļ�¼���������кϲ�,Ȼ���޳�����Mijk (Npt,1,1)��0��ǡ�
                if(mijk(j,1,1)>0) then
			       do k=1,npp
				      If( Nbijk (k,1) == Mijk (j,1,1) .and. &
					      Nbijk (k,2) == Mijk (j,1,2) .and. &
                          Nbijk (k,3) == Mijk (j,1,3) .and. &
						  Nbijk (k,4) == Mijk (j,1,4) )   Mijk (j,1,1) = 0  !�˵��Ѽ��룬Ӧ���޳�

				      If( Nbijk (k,1) == Mijk (j,2,1) .and.  &
					      Nbijk (k,2) == Mijk (j,2,2) .and.  &
                          Nbijk (k,3) == Mijk (j,2,3) .and.  &
						  Nbijk (k,4) == Mijk (j,2,4)  )  Mijk (j,2,1) = 0  !�˵��Ѽ��룬Ӧ���޳�
				   Enddo

				   if(Mijk (j,1,1)==0 .and. Mijk (j,2,1) > 0 ) then
				      Npp=npp+1   !��j��δ��i�г��ֵģ��ϲ���i��
                      Do k =1,4
                         Nbijk (Npp, k)= Mijk (j,2, k)
                      Enddo
                   end If

				   if(Mijk (j,2,1)==0 .and. Mijk (j,1,1) > 0 ) then
                      Npp=npp+1    !��j��δ��i�г��ֵģ��ϲ���i��
                      Do k =1,4
                         Nbijk (Npp, k)= Mijk (j,1, k)
                      Enddo
                      Mijk (j,1,1) = 0
                   end If
				endif
			 enddo
			 if (npp == npp0 ) goto 10000  !�������û�����ӣ����ʾ�������µ����������˵㣬��ktѭ��������
		  enddo
10000     continue
		  npti = npti + 1
          ALLOCATE ( cdate(npti)%Nbijk (Npp,4))   !���ֵ�npti���λ����Ϣ����š������ţ�
		  cdate(npti)%Npp=Npp
          Do j=1, Npp
             Do k =1,4
                cdate(npti)%Nbijk (j,k)=Nbijk (j, k)
             Enddo
          Enddo
	   endif
	enddo
    nptt = npti
    DeALlocate(Nbijk, Mijk)
!	do i=1,nptt
!	    n=cdate(i)%Npp
!	    do j=1,n
!		    write(*,*) i,j,cdate(i)%nbijk(j,1),cdate(i)%nbijk(j,2),cdate(i)%nbijk(j,3),cdate(i)%nbijk(j,4)
!		enddo
!	enddo
    singles=0
    do i=1,nptt
	    npp=cdate(i)%npp
		if(npp>2) then
		    singles=singles+1
		endif
	enddo
	allocate(single(1:singles))

	singles=0
	do ii=1,nptt
	    npp=cdate(ii)%npp
		if(npp>2) then
		    singles=singles+1
			single(singles)%singularity=ii
		    allocate(single(singles)%singularities(1:npp,1:2))
		    do i=1,npp
			    do j=1,2
			        single(singles)%singularities(i,j)=-10000
				enddo
	        enddo

		    do nb=1,nblocks
                nrmax = mb_bc(nb)%nregions         !���鹲��nrmax���߽���Ҫ����
                do nr = 1,nrmax
	                bctype = mb_bc(nb)%bc(nr)%bctype
                    if( bctype == 2 ) then          !��ȫ����/��ȫ�Ǵ߻� ���ȱ���
	                    do m=1,3
                            s_st(m)   = mb_bc(nb)%bc(nr)%s_st(m)   !��ʼ������(�ɶ�����)
                            s_ed(m)   = mb_bc(nb)%bc(nr)%s_ed(m)   !��ֹ������(�ɶ�����)
                        enddo

                        do i = s_st(1),s_ed(1)
                            do j = s_st(2),s_ed(2)
                                do k = s_st(3),s_ed(3)
										do jj=1,npp
										    nbb=cdate(ii)%Nbijk (jj,1)
			                                it=cdate(ii)%Nbijk (jj,2)
			                                jt=cdate(ii)%Nbijk (jj,3)
			                                kt=cdate(ii)%Nbijk (jj,4)
											if(nbb==nb.and.i==it.and.j==jt.and.k==kt) then
											    single(singles)%singularities(jj,1)=jj
												single(singles)%singularities(jj,2)=2
											endif
										enddo
								enddo
							enddo
						enddo
                    endif

					if( bctype == 4 ) then          !��ȫ����/��ȫ�Ǵ߻� ���ȱ���
	                    do m=1,3
                            s_st(m)   = mb_bc(nb)%bc(nr)%s_st(m)   !��ʼ������(�ɶ�����)
                            s_ed(m)   = mb_bc(nb)%bc(nr)%s_ed(m)   !��ֹ������(�ɶ�����)
                        enddo

                        do i = s_st(1),s_ed(1)
                            do j = s_st(2),s_ed(2)
                                do k = s_st(3),s_ed(3)
								   do jj=1,npp
									  nbb=cdate(ii)%Nbijk (jj,1)
			                          it=cdate(ii)%Nbijk (jj,2)
			                          jt=cdate(ii)%Nbijk (jj,3)
			                          kt=cdate(ii)%Nbijk (jj,4)
									  if(nbb==nb.and.i==it.and.j==jt.and.k==kt) then
									     single(singles)%singularities(jj,1)=jj
										 single(singles)%singularities(jj,2)=4
									  endif
								   enddo
							   enddo
							enddo
						enddo
                    endif

                enddo
            enddo
		endif
	enddo
!	do i=1,singles
!	    ii=single(i)%singularity
!	    do j=1,cdate(ii)%npp
!	        jj=single(i)%singularities(j,1)
!			kk=single(i)%singularities(j,2)
!		    write(*,*) ii,jj,kk,cdate(ii)%Nbijk(j,1),cdate(ii)%Nbijk(j,2),cdate(ii)%Nbijk(j,3),cdate(ii)%Nbijk(j,4)
!		enddo
!	enddo

#ifndef PARALLEL
	OPEN(2,FILE=conpointname,FORM="unformatted")

		WRITE(2)NPTT
		do i=1,nptt
		   NPP=cdate(i)%npp
           WRITE(2)NPP

		   DO J=1,NPP
		      WRITE(2)(cdate(i)%Nbijk (j,k),k=1,4)
		   ENDDO
		ENDDO
	CLOSE(2)

!!	OPEN(3,FILE="grid/BC_CONNECT_ASCII.DAT",FORM="FORMATTED")
!!
!!		WRITE(3,*)NPTT
!!		do i=1,nptt
!!		   NPP=cdate(i)%npp
!!           WRITE(3,*)NPP,', POINTS'
!!
!!		   DO J=1,NPP
!!		      WRITE(3,*)(cdate(i)%Nbijk (j,k),k=1,4)
!!		   ENDDO
!!		ENDDO
!!	CLOSE(3)
#endif

	return
end subroutine connect_pre_new
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------!
subroutine boundary_match_dq
!	��1��nptt �Խ����ϵ㣻��2��cdateΪ�ṹ������npp��Nbijk(:,:)����Ԫ��                !
!	��3��npp��ʾ���жԽ����ϵ�����������NbijkΪλ����Ϣ                              !
!	��4��Nbijk(n,1)��Nb��Nbijk(n,2)��i��Nbijk(n,3)��i��Nbijk(n,4)��k��nΪ1��npp        !
!    (5) �ýṹ�ķ�ʽ���������λ�õĽǶȼ�¼�Խ����ϵ��λ��

    use global_variables,only:mb_vol,mb_dq,cdate,nptt,nl
	use define_precision_mod
    implicit none
    integer :: nb,i,j,k,it,jt,kt,npp,m
	real(prec):: dqq(1:nl),vol_p

    do i=1,nptt
	    npp=cdate(i)%npp
		do m=1,nl
		    dqq(m)=0.0_prec
		enddo
	    do j=1,npp
            nb=cdate(i)%Nbijk (j,1)
			it=cdate(i)%Nbijk (j,2)
			jt=cdate(i)%Nbijk (j,3)
			kt=cdate(i)%Nbijk (j,4)
			vol_p=mb_vol(nb)%a3d(it,jt,kt)
			do m=1,nl
			    dqq(m)=dqq(m)+mb_dq(nb)%a4d(m,it,jt,kt)/vol_p
			enddo
		enddo

		do m=1,nl
		    dqq(m)=dqq(m)/(npp*1.0_prec)
		enddo

		do j=1,npp
		    nb=cdate(i)%Nbijk (j,1)
			it=cdate(i)%Nbijk (j,2)
			jt=cdate(i)%Nbijk (j,3)
			kt=cdate(i)%Nbijk (j,4)
			vol_p=mb_vol(nb)%a3d(it,jt,kt)
		    do m=1,nl
			    mb_dq(nb)%a4d(m,it,jt,kt)=dqq(m)*vol_p
			enddo
		enddo
	enddo
	return

end subroutine boundary_match_dq
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------!
subroutine boundary_match_pv
!	��1��nptt �Խ����ϵ㣻��2��cdateΪ�ṹ������npp��Nbijk(:,:)����Ԫ��                !
!	��3��npp��ʾ���жԽ����ϵ�����������NbijkΪλ����Ϣ                              !
!	��4��Nbijk(n,1)��Nb��Nbijk(n,2)��i��Nbijk(n,3)��i��Nbijk(n,4)��k��nΪ1��npp        !
!    (5) �ýṹ�ķ�ʽ���������λ�õĽǶȼ�¼�Խ����ϵ��λ��

    use global_variables,only : mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,cdate,nptt,nl
	use define_precision_mod
    implicit none
    integer :: nb,i,j,k,it,jt,kt,npp,m
	real(prec):: dqq(1:6)

    do i=1,nptt
	    npp=cdate(i)%npp
		do m=1,6
		    dqq(m)=0.0_prec
		enddo
	    do j=1,npp
            nb=cdate(i)%Nbijk (j,1)
			it=cdate(i)%Nbijk (j,2)
			jt=cdate(i)%Nbijk (j,3)
			kt=cdate(i)%Nbijk (j,4)

			    dqq(1)=dqq(1)+mb_r(nb)%a3d(it,jt,kt)
			    dqq(2)=dqq(2)+mb_u(nb)%a3d(it,jt,kt)
			    dqq(3)=dqq(3)+mb_v(nb)%a3d(it,jt,kt)
			    dqq(4)=dqq(4)+mb_w(nb)%a3d(it,jt,kt)
			    dqq(5)=dqq(5)+mb_p(nb)%a3d(it,jt,kt)
			    dqq(6)=dqq(6)+mb_t(nb)%a3d(it,jt,kt)
		enddo

		do m=1,6
		    dqq(m)=dqq(m)/(npp*1.0_prec)
		enddo

		do j=1,npp
		    nb=cdate(i)%Nbijk (j,1)
			it=cdate(i)%Nbijk (j,2)
			jt=cdate(i)%Nbijk (j,3)
			kt=cdate(i)%Nbijk (j,4)
			    mb_r(nb)%a3d(it,jt,kt)=dqq(1)
			    mb_u(nb)%a3d(it,jt,kt)=dqq(2)
			    mb_v(nb)%a3d(it,jt,kt)=dqq(3)
			    mb_w(nb)%a3d(it,jt,kt)=dqq(4)
			    mb_p(nb)%a3d(it,jt,kt)=dqq(5)
			    mb_t(nb)%a3d(it,jt,kt)=dqq(6)
		enddo
	enddo
	return

end subroutine boundary_match_pv
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------!
!   subroutine connect_pre_new1                                                             !
!	���ܣ�ȷ���Խ����������������λ����Ϣ���������ļ�������鹫�У���ÿ�����ϵı�ţ� !
!         ����Ԥ�������֣��ڵõ����������ͱ߽������ļ���͸�������ֵ֮ǰ���ô�ģ�顣 !
!	���룺����������ͱ߽���������Ϣ��nblocks,mb_bc��                                  !
!	������Խ����ϵ��λ����Ϣ��nptt,cdate,single,singles��                            !
!	     ��1��nptt �Խ����ϵ㣻��2��cdateΪ�ṹ������npp��Nbijk(:,:)����Ԫ��           !
!	     ��3��npp��ʾ���жԽ����ϵ�����������NbijkΪλ����Ϣ                         !
!	     ��4��Nbijk(n,1)��Nb��Nbijk(n,2)��i��Nbijk(n,3)��i��Nbijk(n,4)��k��nΪ1��npp   !
!	��ƣ�ëö��                                                                       !
!	���ԣ� Ϳ����                ������                                                !
!   ʱ�䣺2009��2��                                                                    !
!--------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------!
subroutine connect_pre_new1
    use global_variables,only: nblocks,mb_bc,cdate,nptt,single,singles,connect_point,conpointname
#ifdef PARALLEL
    use mod_parallels,only : master,myid
#endif
	implicit none
	integer :: nb,nrmax,nr,bctype,npti,i,j,k,npp,ii,jj,kk,nbt,it,jt,kt,npp0
	integer,dimension(1:3) :: s_st,s_ed,s_lr3d
	integer,pointer :: Nbijk (:,:)   !��ʱ���飬��Խ����������ռ�һ�����Ӧ�ļ�������ϵ�µĵ�������Ϣ
    integer,pointer :: Mijk (:,:,:)  !��ʱ���飬��Խ����������ռ����е��Ӧ�ļ�������ϵ�µĵ�������Ϣ
    integer,pointer :: NWindow (:,:)  !��ʱ���飬��Խ���λ����Ϣ
	integer :: m,nbb,nb0s,nb0t,nb1s,nb1t

!    WRITE(*,*)'�����룺1---��ȡ�Խ����ã�0---�������Խ�����'
!	READ(*,*)IBC_TEM
	IF(connect_point ==1)THEN
	   write(*,*)'��ȡ�Խӱ߽�������������Ա�ƽ��ʱʹ��'
	   OPEN(2,FILE=conpointname,FORM="unformatted")
		   read(2)NPTT
		   allocate( cdate(nptt) )

		   do i=1,nptt
              READ(2)NPP
		      cdate(i)%npp=npp
			  allocate( cdate(i)%nbijk(npp,4) )

		      DO J=1,NPP
		         READ(2)(cdate(i)%Nbijk (j,k),k=1,4)
		      ENDDO
		ENDDO
	   CLOSE(2)

	  RETURN
	ENDIF

#ifdef PARALLEL
    if (myid == master) then
#endif

    write(*,*)'���ڼ���Խӱ߽�������������Ա�ƽ��ʱʹ��'

#ifdef PARALLEL
    end if
#endif

!   ����ȷ��nptt
    nptt=0
	npp0=0
    do nb=1,nblocks
	   nrmax = mb_bc(nb)%nregions
	   do nr = 1,nrmax
	      bctype = mb_bc(nb)%bc(nr)%bctype
		  if( bctype < 0 ) then               !�Խӱ߽��ϸ����
			  i=mb_bc(nb)%bc(nr)%s_st(1) - mb_bc(nb)%bc(nr)%s_ed(1)
			  j=mb_bc(nb)%bc(nr)%s_st(2) - mb_bc(nb)%bc(nr)%s_ed(2)
			  k=mb_bc(nb)%bc(nr)%s_st(3) - mb_bc(nb)%bc(nr)%s_ed(3)
			  nptt = nptt + ( ABS(i)+1 ) * ( ABS(j)+1 ) * ( ABS(k)+1 )
			  npp0 = npp0 + 1  !*TGH. �����һ���ж��ٸ��Խ���
		  endif
	   enddo
    enddo

	ALLOCATE (cdate(nptt/2+1))   !�Խ����ϵĵ㣬�����������������,��ˣ��������ռ俴��ֻ��һ�롣
	ALLOCATE (Mijk(Nptt,2,4))
    ALLOCATE (Nbijk(Nptt,4))
    ALLOCATE (NWindow(npp0,9))

    npp0 = 0  ! ��¼�Խ��洰����
    do nb=1,nblocks
	   nrmax = mb_bc(nb)%nregions
	   do nr = 1,nrmax
	      bctype = mb_bc(nb)%bc(nr)%bctype
		  if( bctype < 0 ) then               !�ԽӴ���λ����Ϣ����š����ں�
			  npp0 = npp0 + 1
			  NWindow(npp0,1) = nb
			  NWindow(npp0,2) = nr
			  NWindow(npp0,3) = mb_bc(nb)%bc(nr)%nbt
			  i=( mb_bc(nb)%bc(nr)%s_st(1) + mb_bc(nb)%bc(nr)%s_ed(1) )/2
			  j=( mb_bc(nb)%bc(nr)%s_st(2) + mb_bc(nb)%bc(nr)%s_ed(2) )/2
			  k=( mb_bc(nb)%bc(nr)%s_st(3) + mb_bc(nb)%bc(nr)%s_ed(3) )/2
			  NWindow(npp0,4) = i                              ! �Խ��洰�����ĵ�λ��
			  NWindow(npp0,5) = j
			  NWindow(npp0,6) = k

			  NWindow(npp0,7) = mb_bc(nb)%bc(nr)%image(i,j,k)   ! �Խ��洰�����ĵ�λ��
			  NWindow(npp0,8) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
			  NWindow(npp0,9) = mb_bc(nb)%bc(nr)%kmage(i,j,k)
		  endif
	   enddo
    enddo

!   �޳��ظ���¼�ĶԽ��洰�ڣ����һ��һ���޳�Ҫ��
    do nr=1,npp0
	   nb0s= NWindow(nr,1)
	   nb0t= NWindow(nr,3)
	   if( nb0s > 0 ) then  !*tgh. ���ǶԽ�����û�б��޳�
	      i  = NWindow(nr,4)
	      j  = NWindow(nr,5)
	      k  = NWindow(nr,6)
	      ii = NWindow(nr,7)
	      jj = NWindow(nr,8)
	      kk = NWindow(nr,9)
          do npp = nr+1,npp0
	         nb1s = NWindow(npp,1)
	         nb1t = NWindow(npp,3)

		     if( (nb0s == nb1t) .and. (nb0t == nb1s) ) then  !*tgh. �����Ӧ��
			     ! �Խ��洰���ظ��������������š����ĵ����

	            it = NWindow(npp,7)
	            jt = NWindow(npp,8)
	            kt = NWindow(npp,9)
			    if((i ==it) .and. (j ==jt) .and. (k ==kt) ) then
				   NWindow(npp,1) = -100
				   NWindow(npp,3) = -100
				   goto 0100
			    else
					write(*,*)'���Ӧ���ˣ����ǵ�û�ж�Ӧ��,��������֮�����������������϶Խ�'
				endif

		     endif
	      enddo
0100      continue
	   endif
	enddo


    Nptt = 0   !      �Խӱ߽�����ı��ϵĵ���
	npti = 0   !      �Խӱ߽����ϵ��ڵ�ĵ���
 !  cdate��Mijk��¼���жԽ����ڵ�ͱ��ϵ�Ŀ�ź�����ָ�� i,j,k
    do nrmax=1,npp0
!	   nb= NWindow(nr,1) !TGH. OLD. �Ƿ���ȷ��Ӧ�ø�Ϊ��ʽ
	   nb = nwindow(nrmax,1) !TGH. �޸�

	   if(nb > 0 ) then
		  nr = NWindow(nrmax,2)
!           =======================���洢һ���Խ���߽�==================================
			    do m=1,3
                   s_st(m)  = mb_bc(nb)%bc(nr)%s_st(m)   !��ʼ������(�ɶ�����)
                   s_ed(m)  = mb_bc(nb)%bc(nr)%s_ed(m)   !��ֹ������(�ɶ�����)
                   s_lr3d(m)= abs( mb_bc(nb)%bc(nr)%s_lr3d(m) )
                enddo

				ii = 1 - s_lr3d(1)
				jj = 1 - s_lr3d(2)
				kk = 1 - s_lr3d(3)
!       �Խӱ߽����ϵ��ڵ㣬ֻ����������������飬ֱ�Ӵ���cdate
			    do i= s_st(1)+ii,s_ed(1)-ii
                    do j= s_st(2)+jj,s_ed(2)-jj
                        do k=s_st(3)+kk,s_ed(3)-kk
                  		    npti = npti + 1
                            ALLOCATE ( cdate(npti)%Nbijk (2,4) )   !���ֵ�npti���λ����Ϣ����š������ţ�
		                    cdate(npti)%npp        = 2
                            cdate(npti)%Nbijk(1,1) = nb
                            cdate(npti)%Nbijk(1,2) = i
                            cdate(npti)%Nbijk(1,3) = j
                            cdate(npti)%Nbijk(1,4) = k

                            cdate(npti)%Nbijk(2,1) = mb_bc(nb)%bc(nr)%nbt
                            cdate(npti)%Nbijk(2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                            cdate(npti)%Nbijk(2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                            cdate(npti)%Nbijk(2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)

          			    enddo
				    enddo
			    enddo
!*tgh. �ڵ�Խ�ֱ�Ӵ��� cdate(npti)%Nbijk(1,:)��cdate(npti)%Nbijk(2,:)��

!      �Խӱ߽�����ı��ϵĵ㣬�ſ����г�����������鹲�д˵����,��Mijk�洢
                do it = 1,s_lr3d(1)           ! i�ű�Ϊ�����ĶԽ�������
				   i  =  s_st(1)
				   jt = max(1,s_ed(2)-s_st(2))
				   kt = max(1,s_ed(3)-s_st(3))

                   do j= s_st(2),s_ed(2)
                      do k=s_st(3),s_ed(3),kt
				         Nptt= Nptt+1
                         Mijk (Nptt,1,1) = nb
                         Mijk (Nptt,1,2) = i
                         Mijk (Nptt,1,3) = j
                         Mijk (Nptt,1,4) = k

                         Mijk (Nptt,2,1) = mb_bc(nb)%bc(nr)%nbt
                         Mijk (Nptt,2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                         Mijk (Nptt,2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                         Mijk (Nptt,2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)

					  enddo
				   enddo
                   do j= s_st(2),s_ed(2),jt
                      do k=s_st(3)+1,s_ed(3)-1
					     Nptt= Nptt+1
                         Mijk (Nptt,1,1) = nb
                         Mijk (Nptt,1,2) = i
                         Mijk (Nptt,1,3) = j
                         Mijk (Nptt,1,4) = k

                         Mijk (Nptt,2,1) = mb_bc(nb)%bc(nr)%nbt
                         Mijk (Nptt,2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                         Mijk (Nptt,2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                         Mijk (Nptt,2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)

					  enddo
				   enddo
				enddo

                do jt = 1,s_lr3d(2)    ! j�ű�Ϊ�����ĶԽ�������
				   j  =  s_st(2)
				   it = max(1,s_ed(1)-s_st(1))
				   kt = max(1,s_ed(3)-s_st(3))

                   do i= s_st(1),s_ed(1)
                      do k=s_st(3),s_ed(3),kt
				         Nptt= Nptt+1
                         Mijk (Nptt,1,1) = nb
                         Mijk (Nptt,1,2) = i
                         Mijk (Nptt,1,3) = j
                         Mijk (Nptt,1,4) = k

                         Mijk (Nptt,2,1) = mb_bc(nb)%bc(nr)%nbt
                         Mijk (Nptt,2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                         Mijk (Nptt,2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                         Mijk (Nptt,2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)

				      enddo
				   enddo
                   do i= s_st(1),s_ed(1),it
                      do k=s_st(3)+1,s_ed(3)-1
					     Nptt= Nptt+1
                         Mijk (Nptt,1,1) = nb
                         Mijk (Nptt,1,2) = i
                         Mijk (Nptt,1,3) = j
                         Mijk (Nptt,1,4) = k

                         Mijk (Nptt,2,1) = mb_bc(nb)%bc(nr)%nbt
                         Mijk (Nptt,2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                         Mijk (Nptt,2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                         Mijk (Nptt,2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)
					  enddo
				   enddo
				enddo

                do kt = 1,s_lr3d(3)   ! k�ű�Ϊ�����ĶԽ�������
				   k  =  s_st(3)
				   it = max(1,s_ed(1)-s_st(1))
				   jt = max(1,s_ed(2)-s_st(2))

                   do i= s_st(1),s_ed(1)
                      do j=s_st(2),s_ed(2),jt
				         Nptt= Nptt+1
                         Mijk (Nptt,1,1) = nb
                         Mijk (Nptt,1,2) = i
                         Mijk (Nptt,1,3) = j
                         Mijk (Nptt,1,4) = k

                         Mijk (Nptt,2,1) = mb_bc(nb)%bc(nr)%nbt
                         Mijk (Nptt,2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                         Mijk (Nptt,2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                         Mijk (Nptt,2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)
					  enddo
				   enddo
                   do i= s_st(1),s_ed(1),it
                      do j=s_st(2)+1,s_ed(2)-1
					     Nptt= Nptt+1
                         Mijk (Nptt,1,1) = nb
                         Mijk (Nptt,1,2) = i
                         Mijk (Nptt,1,3) = j
                         Mijk (Nptt,1,4) = k

                         Mijk (Nptt,2,1) = mb_bc(nb)%bc(nr)%nbt
                         Mijk (Nptt,2,2) = mb_bc(nb)%bc(nr)%image(i,j,k)
                         Mijk (Nptt,2,3) = mb_bc(nb)%bc(nr)%jmage(i,j,k)
                         Mijk (Nptt,2,4) = mb_bc(nb)%bc(nr)%kmage(i,j,k)
				      enddo
				   enddo
				enddo

!*tgh. �ڵ�Խ�ֱ�Ӵ��� cdate(npti)%Nbijk(1,:)��cdate(npti)%Nbijk(2,:)��
!*tgh. npti���������ڵ�Խӵ��ܺ�
!*tgh. ÿ�����4���߶ԽӴ��� Mijk (Nptt,1,:) �� Mijk (Nptt,2,:) ��
!*tgh. nptt�������б߽�Խӵ��ܺ�
!           =======================�洢��һ���Խ���߽�==================================

	   endif
	enddo

!  ������3����������ϵĶԽ���ĵ���кϲ����޳���׼ȷ�õ��������ռ俴�ĶԽ����ϵ��ܵ�����λ����Ϣ
	do i=1,nptt  !*tgh. ���б߽�Խӵ�
	   if(Mijk(i,1,1)>0) then
		  do j=1,4
		     nbijk(1,j)=mijk(i,1,j)
			 nbijk(2,j)=mijk(i,2,j)
		  enddo
		  npp=2           !nppΪ���иõ���������
	!======================================================================================
	!            Ϊ�˲�©������4��ʱ���������飬��ˣ������һ��ѭ��kt                  !
	!======================================================================================
		  do kt=1,nptt !*tgh.  ���б߽�Խӵ�
		     npp0=npp
		     do j=i+1,nptt   !�����ж�������ʶ���ظ���û�У�ȫ�ظ����޳��������ظ���
			                 !��û���ظ��ļ�¼���������кϲ�,Ȼ���޳�����Mijk (Npt,1,1)��0��ǡ�
                if(mijk(j,1,1)>0) then
			       do k=1,npp
				      If( Nbijk (k,1) == Mijk (j,1,1) .and. &
					      Nbijk (k,2) == Mijk (j,1,2) .and. &
                          Nbijk (k,3) == Mijk (j,1,3) .and. &
						  Nbijk (k,4) == Mijk (j,1,4) )   Mijk (j,1,1) = 0  !�˵��Ѽ��룬Ӧ���޳�

				      If( Nbijk (k,1) == Mijk (j,2,1) .and.  &
					      Nbijk (k,2) == Mijk (j,2,2) .and.  &
                          Nbijk (k,3) == Mijk (j,2,3) .and.  &
						  Nbijk (k,4) == Mijk (j,2,4)  )  Mijk (j,2,1) = 0  !�˵��Ѽ��룬Ӧ���޳�
				   Enddo

				   if(Mijk (j,1,1)==0 .and. Mijk (j,2,1) > 0 ) then
				      Npp=npp+1   !��j��δ��i�г��ֵģ��ϲ���i��
                      Do k =1,4
                         Nbijk (Npp, k)= Mijk (j,2, k)
                      Enddo
                   end If

				   if(Mijk (j,2,1)==0 .and. Mijk (j,1,1) > 0 ) then
                      Npp=npp+1    !��j��δ��i�г��ֵģ��ϲ���i��
                      Do k =1,4
                         Nbijk (Npp, k)= Mijk (j,1, k)
                      Enddo
                      Mijk (j,1,1) = 0
                   end If
				endif
			 enddo
			 if (npp == npp0 ) goto 10000  !�������û�����ӣ����ʾ�������µ����������˵㣬��ktѭ��������
		  enddo
10000     continue
		  npti = npti + 1
          ALLOCATE ( cdate(npti)%Nbijk (Npp,4))   !���ֵ�npti���λ����Ϣ����š������ţ�
		  cdate(npti)%Npp=Npp
          Do j=1, Npp
             Do k =1,4
                cdate(npti)%Nbijk (j,k)=Nbijk (j, k)
             Enddo
          Enddo
	   endif
	enddo
    nptt = npti
    DeALlocate(Nbijk, Mijk, NWindow)

!	do i=1,nptt
!	    n=cdate(i)%Npp
!	    do j=1,n
!		    write(*,*) i,j,cdate(i)%nbijk(j,1),cdate(i)%nbijk(j,2),cdate(i)%nbijk(j,3),cdate(i)%nbijk(j,4)
!		enddo
!	enddo
    singles=0
    do i=1,nptt
	    npp=cdate(i)%npp
		if(npp>2) then
		    singles=singles+1
		endif
	enddo
	allocate(single(1:singles))

	singles=0
	do ii=1,nptt
	    npp=cdate(ii)%npp
		if(npp>2) then
		    singles=singles+1
			single(singles)%singularity=ii
		    allocate(single(singles)%singularities(1:npp,1:2))
		    do i=1,npp
			    do j=1,2
			        single(singles)%singularities(i,j)=-10000
				enddo
	        enddo

		    do nb=1,nblocks
                nrmax = mb_bc(nb)%nregions         !���鹲��nrmax���߽���Ҫ����
                do nr = 1,nrmax
	                bctype = mb_bc(nb)%bc(nr)%bctype
                    if( bctype == 2 ) then          !��ȫ����/��ȫ�Ǵ߻� ���ȱ���
	                    do m=1,3
                            s_st(m)   = mb_bc(nb)%bc(nr)%s_st(m)   !��ʼ������(�ɶ�����)
                            s_ed(m)   = mb_bc(nb)%bc(nr)%s_ed(m)   !��ֹ������(�ɶ�����)
                        enddo

                        do i = s_st(1),s_ed(1)
                            do j = s_st(2),s_ed(2)
                                do k = s_st(3),s_ed(3)
										do jj=1,npp
										    nbb=cdate(ii)%Nbijk (jj,1)
			                                it=cdate(ii)%Nbijk (jj,2)
			                                jt=cdate(ii)%Nbijk (jj,3)
			                                kt=cdate(ii)%Nbijk (jj,4)
											if(nbb==nb.and.i==it.and.j==jt.and.k==kt) then
											    single(singles)%singularities(jj,1)=jj
												single(singles)%singularities(jj,2)=2
											endif
										enddo
								enddo
							enddo
						enddo
                    endif

					if( bctype == 4 ) then          !��ȫ����/��ȫ�Ǵ߻� ���ȱ���
	                    do m=1,3
                            s_st(m)   = mb_bc(nb)%bc(nr)%s_st(m)   !��ʼ������(�ɶ�����)
                            s_ed(m)   = mb_bc(nb)%bc(nr)%s_ed(m)   !��ֹ������(�ɶ�����)
                        enddo

                        do i = s_st(1),s_ed(1)
                            do j = s_st(2),s_ed(2)
                                do k = s_st(3),s_ed(3)
								   do jj=1,npp
									  nbb=cdate(ii)%Nbijk (jj,1)
			                          it=cdate(ii)%Nbijk (jj,2)
			                          jt=cdate(ii)%Nbijk (jj,3)
			                          kt=cdate(ii)%Nbijk (jj,4)
									  if(nbb==nb.and.i==it.and.j==jt.and.k==kt) then
									     single(singles)%singularities(jj,1)=jj
										 single(singles)%singularities(jj,2)=4
									  endif
								   enddo
							   enddo
							enddo
						enddo
                    endif

                enddo
            enddo
		endif
	enddo
!	do i=1,singles
!	    ii=single(i)%singularity
!	    do j=1,cdate(ii)%npp
!	        jj=single(i)%singularities(j,1)
!			kk=single(i)%singularities(j,2)
!		    write(*,*) ii,jj,kk,cdate(ii)%Nbijk(j,1),cdate(ii)%Nbijk(j,2),cdate(ii)%Nbijk(j,3),cdate(ii)%Nbijk(j,4)
!		enddo
!	enddo
	OPEN(2,FILE=conpointname,FORM="unformatted")

		WRITE(2)NPTT
		do i=1,nptt
		   NPP=cdate(i)%npp
           WRITE(2)NPP

		   DO J=1,NPP
		      WRITE(2)(cdate(i)%Nbijk (j,k),k=1,4)
		   ENDDO
		ENDDO
	CLOSE(2)

!!	OPEN(3,FILE="grid/BC_CONNECT_ASCII.DAT",FORM="FORMATTED")
!!
!!		WRITE(3,*)NPTT
!!		do i=1,nptt
!!		   NPP=cdate(i)%npp
!!           WRITE(3,*)NPP,', POINTS'
!!
!!		   DO J=1,NPP
!!		      WRITE(3,*)(cdate(i)%Nbijk (j,k),k=1,4)
!!		   ENDDO
!!		ENDDO
!!	CLOSE(3)

	return
end subroutine connect_pre_new1
!==================================================================================
subroutine set_bc_index
    use global_variables,only : nblocks,mb_bc
    implicit none
    integer,parameter :: nbclistmax=6
    integer,parameter :: bclists(nbclistmax)=(/-1,4,5,6,2,3/)
    integer :: nb,nrmax,nr,i,ibl,idbc,bctype
    
    do nb=1,nblocks
       nrmax = mb_bc(nb)%nregions
       allocate(mb_bc(nb)%bcindexs(nrmax))
       
       i = 0
       do ibl=1,nbclistmax
          idbc = bclists(ibl)
          do nr=1,nrmax
             bctype = mb_bc(nb)%bc(nr)%bctype
             if (bctype == idbc) then
                i = i + 1
                mb_bc(nb)%bcindexs(i) = nr
             end if
          end do
       end do
       
       if (i /= nrmax) then
          print *,'ERROR: there are some incorrect BC types!'
          stop
       end if
    end do

end subroutine set_bc_index

!==================================================================================
subroutine check_grid_derivative
    use global_variables,only : nblocks,mb_bc,mb_dim,vol,large,pole_vol,sml_vol,ni,nj,nk,nijk2nd
#ifdef PARALLEL
    use mod_parallels
#endif
    implicit none
    integer :: nb,pnb,nr,bctype,ibeg(3),iend(3)
    integer :: i,j,k,ierr,negvol,negvol_t,is2d,is2d_t,negvol_nb
    real :: minvol,maxvol
    
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
       call recast_grid(nb)
       
       do nr=1,mb_bc(nb)%nregions
          bctype = mb_bc(nb)%bc(nr)%bctype
          ibeg = mb_bc(nb)%bc(nr)%s_st
          iend = mb_bc(nb)%bc(nr)%s_ed
          if (bctype == 71 .or. bctype == 72 .or. bctype == 73) then
             do k=ibeg(3),iend(3)
             do j=ibeg(2),iend(2)
             do i=ibeg(1),iend(1) 
                vol(i,j,k) = pole_vol
             end do
             end do
             end do
          end if
       end do
    end do
    
    
    minvol =  large
    maxvol = -large
    
    negvol = 0
    
    is2d = 0
    
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
       call recast_grid(nb)
       
       if (min(ni,nj,nk) < nijk2nd) then
          is2d = 1
       end if
       
       negvol_nb = 0
       
       do k=1,nk
       do j=1,nj
       do i=1,ni
          minvol = min(minvol, vol(i,j,k))
          maxvol = max(maxvol, vol(i,j,k))
          
          if (vol(i,j,k) < sml_vol) then
             negvol_nb = negvol_nb + 1
             
             if (negvol_nb <= 3) then
                write(*,'(1x,a,4(1x,i5),1x,e12.5)')'Jacobi<0:',nb,i,j,k,vol(i,j,k)
             end if
             
             vol(i,j,k) = sml_vol
          end if
       end do
       end do
       end do
       
       negvol = negvol + negvol_nb
    end do
    
#ifdef PARALLEL
    call MPI_REDUCE(is2d,is2d_t,1,mpi_inprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(negvol,negvol_t,1,mpi_inprec,MPI_SUM,master,MPI_COMM_WORLD,ierr)
#else
    is2d_t = is2d
    negvol_t = negvol
#endif

#ifdef PARALLEL
    if (myid == master) then
#endif
    if (negvol_t > 0) then
       print *,'ERROR(Jacobi<0):',negvol_t
       stop
    else
       print '(1x,a,2(e12.5,a))','$ Interval of Jacobi: (',minvol,',',maxvol,')'
    end if
    
    if (is2d_t > 0) then
       print *,'$ The space dimension is 2D!'
    end if
#ifdef PARALLEL
    end if
#endif

end subroutine check_grid_derivative
!==================================================================================

module mod_tecios
    implicit none
    private
    integer,parameter :: int32=4,float32=4
    character(len=1),parameter :: nullchar=char(0)
    character(len=8),parameter :: cmagicnum='#!TDV102'

    integer(int32),parameter :: irworder=1

    real(float32),parameter :: fzonemark=299.0
    real(float32),parameter :: feohmark=357.0

    integer(int32),pointer :: ititle(:)
    integer(int32) :: inumvar
    integer(int32),pointer :: ivarname(:)

    integer(int32),pointer :: izonename(:)
    integer(int32) :: izonecolor
    integer(int32) :: izonetype
    integer(int32) :: idatapack
    integer(int32) :: ivarloc
    integer(int32) :: ifaceconnect
    integer(int32) :: imax,jmax,kmax
    integer(int32) :: iauxiliar

    integer(int32),pointer :: idataform(:)
    integer(int32) :: ivarshare
    integer(int32) :: ishareconnect

    public :: tecio_open,tecio_ini,tecio_zone
    public :: tecio_eohmark,tecio_data,tecio_close

    contains

    subroutine tecio_open(fileid,filename)
       implicit none
       integer :: fileid,ierr
       character(len=*) :: filename

       open(fileid,file=filename,form='unformatted',iostat=ierr)

    end subroutine tecio_open

    subroutine tecio_ini(fileid,ctitle,nvar,cvarname)
       implicit none
       character(len=*) :: ctitle
       integer :: fileid,nvar
       character(len=*) :: cvarname
       character(len=len(cvarname)) :: cname
       integer :: i,nchar,nsep,ierr

       write(fileid,iostat=ierr)cmagicnum
       write(fileid,iostat=ierr)irworder

       nchar = len_trim(ctitle)
       allocate(ititle(nchar+1),stat=ierr)
       call tecio_ichar(nchar+1,ctitle(1:nchar)//nullchar,ititle(:))
       write(fileid,iostat=ierr)ititle(:)
       deallocate(ititle,stat=ierr)

       inumvar = nvar
       write(fileid,iostat=ierr)inumvar

       nchar = len_trim(cvarname)
       cname = cvarname(1:nchar)
       nsep  = 0
       do i=1,nchar
          if (cname(i:i) == ',') then
             cname(i:i) = nullchar
             nsep = nsep + 1
          end if
       end do
       if (nsep /= inumvar-1) then
          call stop_by_error(201,"tecio varname is invalid")
       end if
       allocate(ivarname(nchar+1),stat=ierr)
       call tecio_ichar(nchar+1,cname(1:nchar)//nullchar,ivarname(:))
       write(fileid,iostat=ierr)ivarname(:)
       deallocate(ivarname,stat=ierr)

    end subroutine tecio_ini

    subroutine tecio_zone(fileid,czonename,ndatapack,ni,nj,nk)
       implicit none
       character(len=*) :: czonename
       integer :: fileid,ndatapack,ni,nj,nk
       integer :: nchar,ierr

       write(fileid,iostat=ierr)fzonemark

       nchar = len_trim(czonename)
       allocate(izonename(nchar+1),stat=ierr)
       call tecio_ichar(nchar+1,czonename(1:nchar)//nullchar,izonename(:))
       write(fileid,iostat=ierr)izonename(:)
       deallocate(izonename,stat=ierr)

       izonecolor = -1          !set to -1 if you want tecplot to determine
       write(fileid,iostat=ierr)izonecolor

       izonetype = 0            !0=ORDERED,1=FELINESEG,2=FETRIANGLE,
                                !3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
       write(fileid,iostat=ierr)izonetype

       idatapack = ndatapack    !0=Block, 1=Point
       write(fileid,iostat=ierr)idatapack

       ivarloc = 0              !0 = Don't specify, all data is
                                !located at the nodes. 1 = Specify
       write(fileid,iostat=ierr)ivarloc

       ifaceconnect = 0         !Number of user defined face neighbor connections (value >= 0)
       write(fileid,iostat=ierr)ifaceconnect

       imax = ni
       jmax = nj
       kmax = nk
       write(fileid,iostat=ierr)imax,jmax,kmax

       iauxiliar = 0            !1=Auxiliary name/value pair to follow
                                !0=No more Auxiliar name/value pairs.
       write(fileid,iostat=ierr)iauxiliar

    end subroutine tecio_zone

    subroutine tecio_eohmark(fileid)
       implicit none
       integer :: fileid,ierr

       write(fileid,iostat=ierr)feohmark

    end subroutine tecio_eohmark

    subroutine tecio_data(fileid,ndataform)
       implicit none
       integer :: ndataform
       integer :: fileid,ierr

       write(fileid,iostat=ierr)fzonemark

       allocate(idataform(inumvar),stat=ierr)
       idataform(:) = ndataform !1=Float,2=Double,3=LongInt,4=ShortInt,
                                !5=Byte,6=Bit
       write(fileid,iostat=ierr)idataform(:)
       deallocate(idataform,stat=ierr)

       ivarshare = 0            !0=no,1=yes
       write(fileid,iostat=ierr)ivarshare

       ishareconnect = -1       !Zone number to share connectivity list with (-1 = no sharing).
       write(fileid,iostat=ierr)ishareconnect

    end subroutine tecio_data

    subroutine tecio_close(fileid)
       implicit none
       integer :: fileid,ierr

       close(fileid,iostat=ierr)

    end subroutine tecio_close

    subroutine tecio_ichar(n,ch,ich)
       implicit none
       integer :: n
       character(len=n) :: ch
       integer(int32) :: ich(n)
       integer :: i

       do i=1,n
          ich(i) = ichar(ch(i:i))
       end do

    end subroutine tecio_ichar

end module mod_tecios

module mod_comp_connects
    implicit none
    
    type simp_pnt_t
       integer :: nbijk(4)
       type(simp_pnt_t),pointer :: next
    end type simp_pnt_t
    
    type comp_pnt_t
       integer :: nsimp_pnts
       type(simp_pnt_t),pointer :: simp_head
       type(comp_pnt_t),pointer :: next
    end type comp_pnt_t
    
    integer :: ncomp_pnts
    type(comp_pnt_t),pointer :: comp_list
    
    contains
    
    subroutine simp_list_insert(simp_head,nbijk,nbijkt,nsimps,success)
       implicit none
       type(simp_pnt_t), pointer :: simp_head
       integer :: nbijk(4),nbijkt(4),nsimps
       logical :: success
       type(simp_pnt_t), pointer :: current,last
       logical :: exist,exist_t
       
       exist = .false.
       current => simp_head
       do while ( associated(current) )
          if (sum( abs(current%nbijk - nbijk) ) == 0) then
             exist = .true.
             exit
          end if
          current => current%next
       end do
       
       exist_t = .false.
       current => simp_head
       do while ( associated(current) )
          if (sum( abs(current%nbijk - nbijkt) ) == 0) then
             exist_t = .true.
             exit
          end if
          current => current%next
       end do
       
       if (exist .or. exist_t) then
          last => simp_head
          do while ( associated(last%next) )
             last => last%next
          end do
       
          if (.not. exist) then
             allocate( last%next )
             current => last%next
             current%next => null()
             current%nbijk = nbijk
             nsimps = nsimps + 1
             
             last => current
          end if
          
          if (.not. exist_t) then
             allocate( last%next )
             current => last%next
             current%next => null()
             current%nbijk = nbijkt
             
             nsimps = nsimps + 1
          end if
          
          success = .true.
       else
          success = .false.
       end if
       
    end subroutine simp_list_insert
    
    subroutine simp_list_destroy( simp_head )
       implicit none
       type(simp_pnt_t), pointer :: simp_head
       type(simp_pnt_t), pointer :: current,next

       current => simp_head
       do while ( associated(current%next) )
          next => current%next
          deallocate( current )
          current => next
       end do
       
       if ( associated(current) ) then
          deallocate( current )
       end if
       
       simp_head => null()
       
    end subroutine simp_list_destroy
    
    
    subroutine comp_list_insert( comp_head,nbijk,nbijkt,ncomps)
       implicit none
       type(comp_pnt_t), pointer :: comp_head
       integer :: nbijk(4), nbijkt(4),ncomps
       logical :: success
       type(comp_pnt_t), pointer :: current,last
       
       success = .false.
       
       current => comp_head
       do while( associated(current) )
          call simp_list_insert(current%simp_head,nbijk,nbijkt,current%nsimp_pnts,success)
          if (success) exit
          current => current%next
       end do
       
       if (.not. success) then
          allocate( current )
          current%next => null()
          
          allocate( current%simp_head )
          current%simp_head%nbijk = nbijk
          allocate( current%simp_head%next )
          current%simp_head%next%next => null()
          current%simp_head%next%nbijk = nbijkt
          
          current%nsimp_pnts = 2
          
          if (.not. associated(comp_head)) then
             comp_head => current
             ncomps = 1
          else
             last => comp_head
             do while ( associated(last%next) )
                last => last%next
             end do
             last%next => current
             ncomps = ncomps + 1
          end if
       end if
       
    end subroutine comp_list_insert
    
    subroutine comp_list_print( comp_head,ncomps)
       implicit none
       type(comp_pnt_t), pointer :: comp_head
       integer :: ncomps
       type(comp_pnt_t), pointer :: current
       type(simp_pnt_t), pointer :: simpcur
       integer :: i,j,nbijk(4)
       
       i = 0
       current => comp_head
       do while ( associated(current) )
          if (current%nsimp_pnts > 2) then
          i = i + 1
          write(*,*)i,current%nsimp_pnts
          
          j = 0
          simpcur => current%simp_head
          do while ( associated(simpcur) )
             j = j + 1
             write(*,*) '    ',j,simpcur%nbijk
             
             simpcur => simpcur%next
          end do
          end if
          
          current => current%next
       end do
       
    end subroutine comp_list_print
    
    function comp_list_count(comp_head) result(ncount)
       implicit none
       type(comp_pnt_t), pointer :: comp_head
       integer :: ncount(2)
       type(comp_pnt_t), pointer :: current
       
       ncount = 0
       
       current => comp_head
       do while ( associated(current) )
          ncount(1) = ncount(1) + 1
          if (current%nsimp_pnts > 2) then
             ncount(2) = ncount(2) + 1
          end if
          current => current%next
       end do
       
    end function comp_list_count    
        
    subroutine comp_list_destroy( comp_head )
       implicit none
       type(comp_pnt_t), pointer :: comp_head
       type(comp_pnt_t), pointer :: current,next

       current => comp_head
       do while ( associated(current%next) )
          next => current%next
          
          call simp_list_destroy( current%simp_head )
          deallocate( current )
          
          current => next
       end do
       
       if ( associated(current) ) then
          call simp_list_destroy( current%simp_head )
          deallocate( current )
       end if
       
       comp_head => null()
       
    end subroutine comp_list_destroy    
   
end module mod_comp_connects

module mod_parallels
#ifdef PARALLEL
    use mpi
    implicit none
    
    integer,parameter :: rksingle = 4  !kind(1.0e0)
    integer,parameter :: rkdouble = 8  !kind(1.0d0)
    integer,parameter :: iksingle = 2  !selected_int_kind(4)
    integer,parameter :: ikdouble = 4  !selected_int_kind(8)

    integer,parameter :: mbsmarker=1976
    
    integer :: reprec,inprec
    integer :: mpi_reprec,mpi_inprec
    integer :: mpi_chprec

    integer,parameter :: master = 0
    
    integer :: myid,numprocs
    integer :: namelen
    character(len=MPI_MAX_PROCESSOR_NAME) :: procname
    
    integer,pointer :: mb_pids(:)
    
    integer :: pnblocks              ! number of computional blocks on processor myid
    integer,pointer :: pnbindexs(:)  ! index of computional blocks on processor myid

    integer :: nvtot
    
    integer :: nppos
    integer,pointer :: ppos_nbijks(:,:)                !! nb,i,j,k
    integer,pointer :: ppos_icdate(:)                  !! ppos_icdate(4,j) <===> cdate(i)
    integer,pointer :: nppos_local(:),ipp_st_local(:)
    
    type cdate_ex_t
       integer :: npp
       integer,pointer :: nbijk(:,:) 
       integer,pointer :: ppos_indexs(:)
    end type cdate_ex_t
    
    integer :: nptt_ex
    type(cdate_ex_t),pointer :: cdate_ex(:)
    
    real,pointer :: dq_npp(:,:)             !! dimension(5,nppos)
    real,pointer :: dq_npp_local(:,:)       !! dimension(5,nppos_local(myid+1))
    real,pointer :: pv_npp(:,:)             !! dimension(5,nppos)
    real,pointer :: pv_npp_local(:,:)       !! dimension(5,nppos_local(myid+1))
    
    type sblock_t
       integer :: nb
       integer :: idim,jdim,kdim
       
       real,pointer :: coor(:,:,:,:)
       real,pointer :: vol(:,:,:)
       real,pointer,dimension(:,:,:) :: kcx,kcy,kcz,kct
       real,pointer,dimension(:,:,:) :: etx,ety,etz,ett
       real,pointer,dimension(:,:,:) :: ctx,cty,ctz,ctt
       
       integer :: nsb
       integer,pointer :: nbsbs(:)
       integer,pointer :: iminsb(:,:),imaxsb(:,:)
    end type sblock_t

    contains
    
    subroutine parallel_init
       implicit none
       real :: val_real
       integer :: val_int
       integer :: ierr
       
       reprec = kind(val_real)
       if (reprec == rksingle) then
          mpi_reprec = MPI_REAL
       else if (reprec == rkdouble) then
          !!mpi_reprec = MPI_REAL8
          mpi_reprec = MPI_DOUBLE_PRECISION
       end if

       inprec = kind(val_int)
       if (inprec == ikdouble) then
          mpi_inprec = MPI_INTEGER
       else if (inprec == iksingle) then
          mpi_inprec = MPI_INTEGER2
       end if
       
       mpi_chprec = MPI_CHARACTER

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
       call MPI_GET_PROCESSOR_NAME(procname,namelen,ierr)
       
       call get_procs_info
       
    end subroutine parallel_init

    subroutine parallel_finalize
       implicit none
       integer :: ierr

       call MPI_FINALIZE(ierr)

    end subroutine parallel_finalize

    subroutine get_procs_info
       implicit none
       character(len=MPI_MAX_PROCESSOR_NAME) :: pcnames(0:numprocs-1)
       character(len=MPI_MAX_PROCESSOR_NAME) :: pcname_recv
       integer :: status(MPI_STATUS_SIZE)
       integer :: sender,np,ierr
       character(5) :: str_id
    
       if (myid == master) then
          pcnames(master) = procname
          do np=0,numprocs-1
             if (np == master) cycle
             call MPI_RECV(pcname_recv,MPI_MAX_PROCESSOR_NAME,mpi_chprec, &
                           MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
             sender = status(MPI_SOURCE)
             pcnames(sender) = pcname_recv
          end do
       else
          call MPI_SEND(procname, &
                        MPI_MAX_PROCESSOR_NAME,mpi_chprec, &
                        master,0,MPI_COMM_WORLD,ierr)
       end if
    
       if (myid == master) then
          write(str_id,'(i5)') numprocs
          write(*,'(3a)') " $ The communicator starts ",trim(adjustl(str_id))," processes."
          write(*,*) "  ======================================================="
          write(*,*) "               The processors name and pid"
          write(*,*) "  -------------------------------------------------------"
          do np=0,numprocs-1
             if (mod(np,3) == 0) write(*,'(a)', advance='no') "   "
             write(str_id,'(i5)') np
             write(*, '(a18,1x)', advance='no') trim(pcnames(np))//'('//trim(adjustl(str_id))//')'
             if (mod(np+1,3) == 0 .or. np==numprocs-1) write(*,*) ''
          end do
          write(*,*) "  -------------------------------------------------------"
          write(*,*) "  The master pid:",master
          write(*,*) "  ======================================================="
       end if
    
       call synchronize(11)
    end subroutine get_procs_info    
    
    subroutine synchronize(ierrcode)
       implicit none
       integer :: ierr,ierrcode

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       if (ierr /= 0) then
          call error_mpi(ierrcode,"MPI processors can't be synchronized")
       end if

    end subroutine synchronize
    
    subroutine send_message(message)
       implicit none
       character(len=*) :: message

       if (myid == master) then
          write(*,*) "$ "//trim(adjustl(message))
       end if

    end subroutine send_message    

    subroutine error_mpi(ierrcode,errmsg)
       implicit none
       integer :: ierrcode
       character(len=*) :: errmsg
       integer :: ierr

       write(*,*) "# "//trim(adjustl(errmsg))//":",ierrcode
       call MPI_ABORT(MPI_COMM_WORLD, ierrcode, ierr)

    end subroutine error_mpi   
    
    subroutine distrib_parameter
       use global_variables
       implicit none
       integer,parameter :: packsize=10000
       integer(8) :: packbuf(packsize)
       integer :: position,m,ierr
    
       if (myid == master) then
          position = 0
          
          !!namelist /inflow/moo,reynolds,attack,sideslip,tref,twall,pref,rref,vref,ndim,height
          call MPI_PACK(moo     ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(reynolds,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(attack  ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(sideslip,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(tref    ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(twall   ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(pref    ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(rref    ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(vref    ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(ndim    ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(height  ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          
          !!namelist /control/nstart,nmethod,nbgmax,ndisk,newton,nomax,nforce,nwerror,method,nplot
          call MPI_PACK(nstart ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nmethod,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nbgmax ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(ndisk  ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(newton ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nomax  ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nforce ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nwerror,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(method ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nplot  ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /step/ntmst,cfl,timedt,timedt_rate,timedt_turb,ndualtst,dtdts,nsubstmx,tolsub
          call MPI_PACK(ntmst ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(cfl   ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(timedt,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(timedt_rate,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(timedt_turb,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(ndualtst,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(dtdts   ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nsubstmx,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(tolsub  ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /flowtype/nvis,nchem,ntmodel
          call MPI_PACK(nvis   ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nchem  ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(ntmodel,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /turbulence/nlamtur,nameturb,nreset,nwallfun,ntrans,xtrans,nwtmax,kmaxlim,ncmpcor,ndes
          call MPI_PACK(nlamtur ,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nameturb,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nreset  ,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nwallfun,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(ntrans  ,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(xtrans  ,1  ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nwtmax  ,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(kmaxlim ,1  ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(ncmpcor ,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(ndes    ,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          
          !!namelist /technic/nlhs,nscheme,nlimiter,efix,csrv,nsmooth,nflux
          call MPI_PACK(nlhs    ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nscheme ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nlimiter,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(efix    ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(csrv    ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nsmooth ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nflux   ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /interplate/xk,xb,c_k2,c_k4
          call MPI_PACK(xk  ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(xb  ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(c_k2,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(c_k4,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /filename/flowname,gridname,bcname,forcename,errhis,tecname
          call MPI_PACK(flowname ,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(gridname ,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(bcname   ,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(forcename,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(errhis   ,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(tecname  ,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /force/nwholefield,sref,lfref,lref,xref,yref,zref
          call MPI_PACK(nwholefield,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(sref       ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(lfref      ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(lref       ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(xref       ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(yref       ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(zref       ,1,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /chemical/gasmodel,nchem_source,nchem_rad
          call MPI_PACK(gasmodel    ,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nchem_source,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(nchem_rad   ,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

          !!namelist /gcl_cic/gcl,cic1,cic2,cic3,cic4,cic5
          call MPI_PACK(gcl ,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(cic1,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(cic2,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(cic3,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(cic4,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(cic5,1,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          
          !!namelist /connect_0/connect_point,conpointname,connect_order,dis_tao
          call MPI_PACK(connect_point,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(conpointname ,120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(connect_order,1  ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          call MPI_PACK(dis_tao      ,1  ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
       end if
       
       call MPI_BCAST(position,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)
       
       if (myid /= master) then
          position = 0

          !!namelist /inflow/moo,reynolds,attack,sideslip,tref,twall,pref,rref,vref,ndim,height
          call MPI_UNPACK(packbuf,packsize,position,moo     ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,reynolds,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,attack  ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,sideslip,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,tref    ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,twall   ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,pref    ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,rref    ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,vref    ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,ndim    ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,height  ,1,mpi_reprec,MPI_COMM_WORLD,ierr)

          !!namelist /control/nstart,nmethod,nbgmax,ndisk,newton,nomax,nforce,nwerror,method,nplot
          call MPI_UNPACK(packbuf,packsize,position,nstart ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nmethod,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nbgmax ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,ndisk  ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,newton ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nomax  ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nforce ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nwerror,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,method ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nplot  ,1,mpi_inprec,MPI_COMM_WORLD,ierr)

          !!namelist /step/ntmst,cfl,timedt,timedt_rate,timedt_turb,ndualtst,dtdts,nsubstmx,tolsub
          call MPI_UNPACK(packbuf,packsize,position,ntmst ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,cfl   ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,timedt,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,timedt_rate,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,timedt_turb,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,ndualtst,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,dtdts   ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nsubstmx,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,tolsub  ,1,mpi_reprec,MPI_COMM_WORLD,ierr)

          !!namelist /flowtype/nvis,nchem,ntmodel
          call MPI_UNPACK(packbuf,packsize,position,nvis   ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nchem  ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,ntmodel,1,mpi_inprec,MPI_COMM_WORLD,ierr)

          !!namelist /turbulence/nlamtur,nameturb,nreset,nwallfun,ntrans,xtrans,nwtmax,kmaxlim,ncmpcor,ndes
          call MPI_UNPACK(packbuf,packsize,position,nlamtur ,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nameturb,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nreset  ,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nwallfun,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,ntrans  ,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,xtrans  ,1  ,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nwtmax  ,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,kmaxlim ,1  ,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,ncmpcor ,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,ndes    ,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)

          !!namelist /technic/nlhs,nscheme,nlimiter,efix,csrv,nsmooth,nflux
          call MPI_UNPACK(packbuf,packsize,position,nlhs    ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nscheme ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nlimiter,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,efix    ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,csrv    ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nsmooth ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nflux   ,1,mpi_inprec,MPI_COMM_WORLD,ierr)

          !!namelist /interplate/xk,xb,c_k2,c_k4
          call MPI_UNPACK(packbuf,packsize,position,xk  ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,xb  ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,c_k2,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,c_k4,1,mpi_reprec,MPI_COMM_WORLD,ierr)

          !!namelist /filename/flowname,gridname,bcname,forcename,errhis,tecname
          call MPI_UNPACK(packbuf,packsize,position,flowname ,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,gridname ,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,bcname   ,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,forcename,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,errhis   ,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,tecname  ,120,mpi_chprec,MPI_COMM_WORLD,ierr)

          !!namelist /force/nwholefield,sref,lfref,lref,xref,yref,zref
          call MPI_UNPACK(packbuf,packsize,position,nwholefield,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,sref       ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,lfref      ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,lref       ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,xref       ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,yref       ,1,mpi_reprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,zref       ,1,mpi_reprec,MPI_COMM_WORLD,ierr)

          !!namelist /chemical/gasmodel,nchem_source,nchem_rad
          call MPI_UNPACK(packbuf,packsize,position,gasmodel    ,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nchem_source,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,nchem_rad   ,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)

          !!namelist /gcl_cic/gcl,cic1,cic2,cic3,cic4,cic5
          call MPI_UNPACK(packbuf,packsize,position,gcl ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,cic1,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,cic2,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,cic3,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,cic4,1,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,cic5,1,mpi_inprec,MPI_COMM_WORLD,ierr)

          !!namelist /connect_0/connect_point,conpointname,connect_order,dis_tao
          call MPI_UNPACK(packbuf,packsize,position,connect_point,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,conpointname ,120,mpi_chprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,connect_order,1  ,mpi_inprec,MPI_COMM_WORLD,ierr)
          call MPI_UNPACK(packbuf,packsize,position,dis_tao      ,1  ,mpi_reprec,MPI_COMM_WORLD,ierr)
       end if
       
       !! gasmodel
       if( gasmodel == 'air.dat') then
          if (myid == master) then
             position = 0
             
             call MPI_PACK(ns     ,1     ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             call MPI_PACK(nf     ,1     ,mpi_inprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             
             do m=1,ns
                call MPI_PACK(varname(m),120,mpi_chprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             end do
             call MPI_PACK(cn_init,ns    ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             call MPI_PACK(ws     ,ns    ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             
             call MPI_PACK(gama   ,1     ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             call MPI_PACK(prl    ,1     ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             call MPI_PACK(prt    ,1     ,mpi_reprec,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          end if
          
          call MPI_BCAST(position,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)
          
          if (myid /= master) then
             position = 0
       
             !!namelist /inflow/moo,reynolds,attack,sideslip,tref,twall,pref,rref,vref,ndim,height
             call MPI_UNPACK(packbuf,packsize,position,ns     ,1     ,mpi_inprec,MPI_COMM_WORLD,ierr)
             call MPI_UNPACK(packbuf,packsize,position,nf     ,1     ,mpi_inprec,MPI_COMM_WORLD,ierr)
             
             allocate( varname(ns) )
             allocate( cn_init(ns) )
             allocate( ws(ns) , ms(ns) , ws1(ns) , ms1(ns) )
             
             do m=1,ns
                call MPI_UNPACK(packbuf,packsize,position,varname(m),120,mpi_chprec,MPI_COMM_WORLD,ierr)
             end do
             call MPI_UNPACK(packbuf,packsize,position,cn_init,ns    ,mpi_reprec,MPI_COMM_WORLD,ierr)
             call MPI_UNPACK(packbuf,packsize,position,ws     ,ns    ,mpi_reprec,MPI_COMM_WORLD,ierr)
             
             call MPI_UNPACK(packbuf,packsize,position,gama   ,1     ,mpi_reprec,MPI_COMM_WORLD,ierr)
             call MPI_UNPACK(packbuf,packsize,position,prl    ,1     ,mpi_reprec,MPI_COMM_WORLD,ierr)
             call MPI_UNPACK(packbuf,packsize,position,prt    ,1     ,mpi_reprec,MPI_COMM_WORLD,ierr)
          end if
       else
          call error_mpi(21,"The gas model is not supportted")
       end if
                       
       call send_message("Distribute the parameters successfully!")
       
    end subroutine distrib_parameter
    
    subroutine read_bc_parallel
       use global_variables
       implicit none
       integer :: flow_solver_id,number_of_blocks,nb,nbt,bctype
       integer :: ndif,nr,nrmax,ntemp,js1,js2,ks1,ks2,ls1,ls2
       integer :: m,n,i,j,k,imax,jmax,kmax
       integer :: co,op(2),t_coor(2),ibcwin
       integer :: s_nd,s_lr,s_fix,idelt,jdelt,kdelt
       integer :: t_nd,t_lr,t_fix,s_t_dirction(3,3)
       integer :: s_index,t_index,s_sign(3),t_sign(3),st_sign(3)
       integer,dimension(1:3) :: s_st,s_ed,t_st,t_ed
       character(len=120) :: blockname
       integer :: pid,nprocs,ntms,ierr
       
       open(101,file=bcname,status='old')
       
       read(101,*)flow_solver_id
       if (flow_solver_id /= mbsmarker) then
          call error_mpi(flow_solver_id,"the file BC is incorrect in parallel computing")
       end if
                   
       read(101, *) nprocs
       if (nprocs > 0 .and. mod(nprocs,numprocs) == 0) then
          ntms = nprocs/numprocs
       else
          ntms = -1
       end if
       
       read(101,*,iostat=ierr)number_of_blocks          
       write(*,*)'读边界文件信息'
       if ( number_of_blocks /= nblocks ) then
          write(*,*)'边界条件文件与网格文件中分块数不一致!!'
          stop
       endif
       
       allocate( mb_bc(nblocks) )
       allocate( mb_pids(nblocks) )
       
       
       do nb=1,nblocks
          read(101,*,iostat=ierr) pid
          if (ntms > 0) then
             pid = (pid-1)/ntms + 1
          end if
          mb_pids(nb) = pid
       
          read(101,*)imax,jmax,kmax         !��ÿ�������ά��
          ndif = abs(imax - mb_dim(nb,1)) + abs(jmax - mb_dim(nb,2)) + abs(kmax - mb_dim(nb,3))
          
          if ( ndif /= 0 ) then
             write(*,*)'第',nb,'块边界条件文件的网格维数与网格文件中网格维数不一致!'
             stop
          endif
          read(101,*)blockname              !��ÿ��Ŀ���
          read(101,*)nrmax                  !��ÿ�������(�߽�)��
          mb_bc(nb)%nregions = nrmax
          
          allocate( mb_bc(nb)%bc(nrmax) )
          do nr=1,nrmax
             read(101,*)s_st(1),s_ed(1),s_st(2),s_ed(2),s_st(3),s_ed(3),bctype
             mb_bc(nb)%bc(nr)%bctype = bctype
             mb_bc(nb)%bc(nr)%nbs = nb     !������е���࣬����������д�ϰ�
             mb_bc(nb)%bc(nr)%s_st = s_st
             mb_bc(nb)%bc(nr)%s_ed = s_ed
             if ( bctype < 0 ) then !�����ƴ����
                read(101,*)t_st(1),t_ed(1),t_st(2),t_ed(2),t_st(3),t_ed(3),nbt,ibcwin
                mb_bc(nb)%bc(nr)%nbt = nbt
                mb_bc(nb)%bc(nr)%ibcwin = ibcwin
                mb_bc(nb)%bc(nr)%t_st = t_st
                mb_bc(nb)%bc(nr)%t_ed = t_ed
             endif
          enddo
       enddo
     
       close(101)
       write(*,*)'finished reading bc info'

    end subroutine read_bc_parallel
   
    
    subroutine distrib_grid
       use global_variables
       implicit none
       integer :: nb,pid,packsize,ierr
       integer :: status(MPI_STATUS_SIZE)
       
       ! broadcast the total number of blocks include each processor
       call MPI_BCAST(nblocks,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
       
       ! broadcast the max number of grid points in three directions
       call MPI_BCAST(nmax,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
       
       if (myid /= master) then
          allocate( mb_x(nblocks), mb_y(nblocks), mb_z(nblocks) )
          allocate( mb_dim(nblocks,3) )
          allocate( mb_pids(nblocks) )
       end if
       call MPI_BCAST(mb_dim ,nblocks*3,mpi_inprec,master,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(mb_pids,nblocks  ,mpi_inprec,master,MPI_COMM_WORLD,ierr)
       
       do nb=1,nblocks
          pid = mb_pids(nb)
          if (pid-1 /= master) then
             packsize = mb_dim(nb,1)*mb_dim(nb,2)*mb_dim(nb,3)

             ! send the coor of mesh to each processor
             ! except master processor(myid=0)
             if (myid == master) then
                call MPI_SEND(mb_x(nb)%a3d,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,ierr)
                deallocate(mb_x(nb)%a3d,stat=ierr)
                
                call MPI_SEND(mb_y(nb)%a3d,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,ierr)
                deallocate(mb_y(nb)%a3d,stat=ierr)
                
                call MPI_SEND(mb_z(nb)%a3d,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,ierr)
                deallocate(mb_z(nb)%a3d,stat=ierr)
             end if

             ! receive the coor of mesh of processor(pid-1/=0)
             ! from master processor(myid=0)
             if (myid == pid-1) then
                allocate( mb_x(nb)%a3d(mb_dim(nb,1),mb_dim(nb,2),mb_dim(nb,3)) ) 
                call MPI_RECV(mb_x(nb)%a3d,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
                              
                allocate( mb_y(nb)%a3d(mb_dim(nb,1),mb_dim(nb,2),mb_dim(nb,3)) ) 
                call MPI_RECV(mb_y(nb)%a3d,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
                              
                allocate( mb_z(nb)%a3d(mb_dim(nb,1),mb_dim(nb,2),mb_dim(nb,3)) ) 
                call MPI_RECV(mb_z(nb)%a3d,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
             end if
          end if
          
          call synchronize(31)
       end do
       
       call send_message("Distribute the mesh successfully!")
        
    end subroutine distrib_grid
    
    subroutine distrib_bc
       use global_variables
       implicit none
       integer,parameter :: packsize=10000
       integer(8) :: packbuf(packsize)
       integer :: position,nb,nr,ierr
    
       if (myid /= master) then
          allocate(mb_bc(nblocks),stat=ierr)
       end if

       ! distribute the info. of blocks<*>
       do nb=1,nblocks
          if (myid == master) then
             position = 0
             call MPI_PACK(mb_bc(nb)%nregions,1,mpi_inprec, &
                           packbuf,packsize,position,MPI_COMM_WORLD,ierr)
          end if
    
          ! broadcast the info. of blocks to each processor
          ! from master processor(myid=0)
          call MPI_BCAST(position,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)
    
          if (myid /= master) then
             position = 0
             
             call MPI_UNPACK(packbuf,packsize,position, &
                             mb_bc(nb)%nregions,1,mpi_inprec,MPI_COMM_WORLD,ierr)
             allocate(mb_bc(nb)%bc(mb_bc(nb)%nregions),stat=ierr)
          end if

          call synchronize(41)
       end do

       do nb=1,nblocks
          do nr=1,mb_bc(nb)%nregions
             if (myid == master) then
                position = 0
                
                call MPI_PACK(mb_bc(nb)%bc(nr)%bctype,1,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(mb_bc(nb)%bc(nr)%nbs   ,1,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(mb_bc(nb)%bc(nr)%s_st  ,3,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(mb_bc(nb)%bc(nr)%s_ed  ,3,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                              
                call MPI_PACK(mb_bc(nb)%bc(nr)%nbt   ,1,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(mb_bc(nb)%bc(nr)%ibcwin,1,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(mb_bc(nb)%bc(nr)%t_st  ,3,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(mb_bc(nb)%bc(nr)%t_ed  ,3,mpi_inprec, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
             end if
    
             ! broadcast the info. of BCs to each processor
             ! from master processor(myid=0)
             call MPI_BCAST(position,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
             call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)
    
             if (myid /= master) then
                position = 0
                
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%bctype,1,mpi_inprec,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%nbs   ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%s_st  ,3,mpi_inprec,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%s_ed  ,3,mpi_inprec,MPI_COMM_WORLD,ierr)
                                
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%nbt   ,1,mpi_inprec,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%ibcwin,1,mpi_inprec,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%t_st  ,3,mpi_inprec,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                mb_bc(nb)%bc(nr)%t_ed  ,3,mpi_inprec,MPI_COMM_WORLD,ierr)
             end if
    
             call synchronize(51)
          end do
       end do
       
       call analyze_bc
       
       call get_local_info
       
       call alloc_bc_comm_array
    
       call send_message("Distribute the topology successfully!")
    
    end subroutine distrib_bc
    
    subroutine get_local_info
       use global_variables
       implicit none
       integer :: nb,pid,pnb,ierr
    
       pnblocks = 0
       do nb=1,nblocks
          pid = mb_pids(nb)
          if (myid == pid-1) then
             pnblocks = pnblocks + 1
          end if
       end do
    
       allocate(pnbindexs(pnblocks),stat=ierr)
       
       pnb = 0
       do nb=1,nblocks
          pid = mb_pids(nb)
          if (myid == pid-1) then
             pnb = pnb + 1
             pnbindexs(pnb) = nb
          end if
       end do
    
    end subroutine get_local_info    
    
    subroutine analyze_bc
       use global_variables
       implicit none
       integer :: flow_solver_id,number_of_blocks,nb,nbt,bctype
       integer :: ndif,nr,nrmax,ntemp,js1,js2,ks1,ks2,ls1,ls2
       integer :: m,n,i,j,k,imax,jmax,kmax
       integer :: co,op(2),t_coor(2),ibcwin
       integer :: s_nd,s_lr,s_fix,idelt,jdelt,kdelt
       integer :: t_nd,t_lr,t_fix,s_t_dirction(3,3)
       integer :: s_index,t_index,s_sign(3),t_sign(3),st_sign(3)
       integer,dimension(1:3) :: s_st,s_ed,t_st,t_ed
       character(len=120) :: blockname
       integer :: pid,nprocs,ntms,ierr
       
       do nb=1,nblocks
          do nr=1,mb_bc(nb)%nregions
             bctype = mb_bc(nb)%bc(nr)%bctype
             s_st = mb_bc(nb)%bc(nr)%s_st
             s_ed = mb_bc(nb)%bc(nr)%s_ed
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
                nbt = mb_bc(nb)%bc(nr)%nbt
                t_st = mb_bc(nb)%bc(nr)%t_st
                t_ed = mb_bc(nb)%bc(nr)%t_ed
                
                call analyze_bc_connect(nb,nr,s_nd,s_st,s_ed,t_st,t_ed,nbt)
                
                do m=1,3
                   mb_bc(nb)%bc(nr)%t_st(m) = min( abs(t_st(m)),abs(t_ed(m)) )
                   mb_bc(nb)%bc(nr)%t_ed(m) = max( abs(t_st(m)),abs(t_ed(m)) )
                end do
                
!!                if (myid == master) then
!!                print '(1x,8I4)',(mb_bc(nb)%bc(nr)%s_st(m),mb_bc(nb)%bc(nr)%s_ed(m),m=1,3),mb_bc(nb)%bc(nr)%s_nd,mb_bc(nb)%bc(nr)%s_lr
!!                print '(1x,8I4)',(mb_bc(nb)%bc(nr)%t_st(m),mb_bc(nb)%bc(nr)%t_ed(m),m=1,3),mb_bc(nb)%bc(nr)%t_nd,mb_bc(nb)%bc(nr)%t_lr
!!                print *,"--------------------------"
!!                end if
             endif
          enddo
       enddo
       
       call synchronize(61)
    
    end subroutine analyze_bc
    
    subroutine analyze_bc_connect(nb,nr,s_nd,s_st,s_ed,t_st,t_ed,nbt)                                         
       use global_variables                                                                        
       implicit none                                                                               
       integer :: ndif,nr,nrmax,ntemp,js1,js2,ks1,ks2,ls1,ls2                                      
       integer :: m,n,i,j,k,imax,jmax,kmax,nbt,nb                                                  
       integer :: co,op(2),t_coor(2)                                                               
       integer :: s_nd,s_lr,s_fix,idelt,jdelt,kdelt                                                
       integer :: t_nd,t_lr,t_fix,s_t_dirction(3,3)                                                
       integer :: s_index,t_index,s_sign(3),t_sign(3),st_sign(3)                                   
       integer,dimension(1:3) :: s_st,s_ed,t_st,t_ed                                               
                                                                                                   
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
                                                                                                   
10     continue                                                                                    
                                                                                                   
       ! ȷ���ű�������1 ����-1 �������ű��Ϊ��                                                    
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
    
    end subroutine analyze_bc_connect                                                                  


    subroutine alloc_bc_comm_array
       use global_variables                                                                        
       implicit none
       integer :: nb,pnb,ierr
       integer :: nr,nregions,ibctype
       integer :: nbt,nrt,id_src,id_des,idir2
       integer,parameter :: mcyc(5) = (/1,2,3,1,2/)
       integer :: ibeg(3),iend(3),idir,inrout,m
       integer :: nvtur
       
       if (nlamtur >= 0 ) then
          nvtur = nlamtur + 1
       else
          nvtur = 0
       end if
       nvtot = 5 + nvtur
       
       do nb=1,nblocks
          id_src = mb_pids(nb)
          nregions = mb_bc(nb)%nregions
          do nr=1,nregions
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if (ibctype < 0) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
                if (id_src /= id_des) then
                   nrt = mb_bc(nb)%bc(nr)%ibcwin
                   if (myid == id_src-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
                      
!!                      allocate( mb_bc(nb)%bc(nr)%vol(ibeg(1):iend(1),ibeg(2):iend(2), &
!!                                                     ibeg(3):iend(3)),stat=ierr)  
!!                      allocate( mb_bc(nb)%bc(nr)%sxyz(ibeg(1):iend(1),ibeg(2):iend(2), &   !!2010
!!                                                      ibeg(3):iend(3),12),stat=ierr)  
                      allocate( mb_bc(nb)%bc(nr)%dis(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                     ibeg(3):iend(3)),stat=ierr)  

                      
                      allocate( mb_bc(nb)%bc(nr)%dqv(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                     ibeg(3):iend(3),5),stat=ierr)  
                                                     
                      allocate( mb_bc(nb)%bc(nr)%qpv_t(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                       ibeg(3):iend(3),6),stat=ierr) 

                      allocate( mb_bc(nb)%bc(nr)%dqv_t(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                     ibeg(3):iend(3),5),stat=ierr)  
                                                     

                      ibeg(idir) = ibeg(idir) - 5*max(inrout,0)
                      iend(idir) = iend(idir) - 5*min(inrout,0)
                      allocate( mb_bc(nb)%bc(nr)%vol(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                     ibeg(3):iend(3)),stat=ierr)  
                      allocate( mb_bc(nb)%bc(nr)%sxyz(ibeg(1):iend(1),ibeg(2):iend(2), &   !!2010
                                                      ibeg(3):iend(3),12),stat=ierr)  
                    
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
!!                      ibeg(idir) = ibeg(idir) - 2*max(inrout,0) !!+ min(inrout,0)
!!                      iend(idir) = iend(idir) - 2*min(inrout,0) !!+ max(inrout,0)
                      ibeg(idir) = ibeg(idir) - 4*max(inrout,0)
                      iend(idir) = iend(idir) - 4*min(inrout,0)
                      
                      do m=1,2               !! it is doubtful.���ܵĻ��������Ȳ�ֵ�ٴ��ݡ�
                         idir2 = mcyc(idir+m)
                         if ( ibeg(idir2) > 1 ) then
                            ibeg(idir2) = ibeg(idir2) - 1
                         end if
                         if ( iend(idir2) < mb_dim(nb,idir2) ) then
                            iend(idir2) = iend(idir2) + 1
                         end if
                      end do

                      allocate( mb_bc(nb)%bc(nr)%qpv(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                     ibeg(3):iend(3),nvtot),stat=ierr) 
                   end if                               
                               
                   if (myid == id_des-1) then                                  
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
                      
!!                      allocate( mb_bc(nbt)%bc(nrt)%volpack(ibeg(1):iend(1),ibeg(2):iend(2), &
!!                                                           ibeg(3):iend(3)),stat=ierr)
!!                      allocate( mb_bc(nbt)%bc(nrt)%sxyzpack(ibeg(1):iend(1),ibeg(2):iend(2), &   !!2010
!!                                                            ibeg(3):iend(3),12),stat=ierr)
                      allocate( mb_bc(nbt)%bc(nrt)%dispack(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                           ibeg(3):iend(3)),stat=ierr)
                                                         
                      allocate( mb_bc(nbt)%bc(nrt)%dqvpack(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                           ibeg(3):iend(3),5),stat=ierr)
                                                         
                      allocate( mb_bc(nbt)%bc(nrt)%qpvpack_t(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                             ibeg(3):iend(3),6),stat=ierr)
                                                         
                      allocate( mb_bc(nbt)%bc(nrt)%dqvpack_t(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                             ibeg(3):iend(3),5),stat=ierr)
                                                           
 
                      ibeg(idir) = ibeg(idir) - 5*max(inrout,0) 
                      iend(idir) = iend(idir) - 5*min(inrout,0)
                      allocate( mb_bc(nbt)%bc(nrt)%volpack(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                           ibeg(3):iend(3)),stat=ierr)
                      allocate( mb_bc(nbt)%bc(nrt)%sxyzpack(ibeg(1):iend(1),ibeg(2):iend(2), &   !!2010
                                                            ibeg(3):iend(3),12),stat=ierr)

                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
!!                      ibeg(idir) = ibeg(idir) - 2*max(inrout,0) 
!!                      iend(idir) = iend(idir) - 2*min(inrout,0)
                      ibeg(idir) = ibeg(idir) - 4*max(inrout,0) 
                      iend(idir) = iend(idir) - 4*min(inrout,0)
                      
                      do m=1,2
                         idir2 = mcyc(idir+m)
                         if ( ibeg(idir2) > 1 ) then
                            ibeg(idir2) = ibeg(idir2) - 1
                         end if
                         if ( iend(idir2) < mb_dim(nb,idir2) ) then
                            iend(idir2) = iend(idir2) + 1
                         end if
                      end do
                      
                      allocate( mb_bc(nbt)%bc(nrt)%qpvpack(ibeg(1):iend(1),ibeg(2):iend(2), &
                                                           ibeg(3):iend(3),nvtot),stat=ierr)
                   end if
                end if
             end if
          end do
       end do
       
    end subroutine alloc_bc_comm_array
    
    subroutine exchange_bc_geometric
       use global_variables                                                                        
       implicit none
       integer :: nb,pnb,i,j,k,m,ierr
       integer :: nr,nregions,ibctype
       integer :: iseq,nbt,nrt,packsize
       integer :: id_src,id_des,tag_seq
       integer :: ibeg(3),iend(3),idir,inrout
       integer :: ist,jst,kst,s_lr3d(3)
       integer :: status(MPI_STATUS_SIZE)
       real :: volp,dis1
       
       iseq = 0
       do nb=1,nblocks
          id_src = mb_pids(nb)
          nregions = mb_bc(nb)%nregions
          do nr=1,nregions
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if (ibctype < 0) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
                if (id_src /= id_des) then
                   iseq = iseq + 1
                   tag_seq = iseq
                   nrt = mb_bc(nb)%bc(nr)%ibcwin
                   
                   if (myid == id_src-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr

                      s_lr3d = mb_bc(nb)%bc(nr)%s_lr3d
                      
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1) 
                         ist = i - s_lr3d(1)
                         jst = j - s_lr3d(2)
                         kst = k - s_lr3d(3)
			             dis1 = sqrt( (mb_x(nb)%a3d(i , j , k ) - mb_x(nb)%a3d(ist,jst,kst))**2 + &
						              (mb_y(nb)%a3d(i , j , k ) - mb_y(nb)%a3d(ist,jst,kst))**2 + &
						              (mb_z(nb)%a3d(i , j , k ) - mb_z(nb)%a3d(ist,jst,kst))**2 )
                         mb_bc(nb)%bc(nr)%dis(i,j,k) = dis1
                      end do
                      end do
                      end do
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)
                      call MPI_SEND(mb_bc(nb)%bc(nr)%dis,packsize,mpi_reprec, &
                                    id_des-1,tag_seq+200,MPI_COMM_WORLD,ierr)
                                    
                      ibeg(idir) = ibeg(idir) - 5*max(inrout,0)
                      iend(idir) = iend(idir) - 5*min(inrout,0)
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1) 
                         volp = mb_vol(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%vol(i,j,k) = volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,1) = mb_kcx(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,2) = mb_kcy(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,3) = mb_kcz(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,4) = mb_etx(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,5) = mb_ety(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,6) = mb_etz(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,7) = mb_ctx(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,8) = mb_cty(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,9) = mb_ctz(nb)%a3d(i,j,k) !!/ volp
                         
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,10) = mb_kct(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,11) = mb_ett(nb)%a3d(i,j,k) !!/ volp
                         mb_bc(nb)%bc(nr)%sxyz(i,j,k,12) = mb_ctt(nb)%a3d(i,j,k) !!/ volp
                      end do
                      end do
                      end do
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)
                      call MPI_SEND(mb_bc(nb)%bc(nr)%vol,packsize,mpi_reprec, &
                                    id_des-1,tag_seq,MPI_COMM_WORLD,ierr)
                      call MPI_SEND(mb_bc(nb)%bc(nr)%sxyz,packsize*12,mpi_reprec, &
                                    id_des-1,tag_seq+100,MPI_COMM_WORLD,ierr)
                                    
                      deallocate( mb_bc(nb)%bc(nr)%sxyz,stat=ierr)  
                      deallocate( mb_bc(nb)%bc(nr)%vol,stat=ierr)  
                   end if
                   
                   if (myid == id_des-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)
                      call MPI_RECV(mb_bc(nbt)%bc(nrt)%dispack,packsize,mpi_reprec, &
                                    id_src-1,tag_seq+200,MPI_COMM_WORLD,status,ierr)
                      
                      ibeg(idir) = ibeg(idir) - 5*max(inrout,0)
                      iend(idir) = iend(idir) - 5*min(inrout,0)
                      packsize = product(iend(:)-ibeg(:)+1,1)
                      call MPI_RECV(mb_bc(nbt)%bc(nrt)%volpack,packsize,mpi_reprec, &
                                    id_src-1,tag_seq,MPI_COMM_WORLD,status,ierr)
                      call MPI_RECV(mb_bc(nbt)%bc(nrt)%sxyzpack,packsize*12,mpi_reprec, &
                                    id_src-1,tag_seq+100,MPI_COMM_WORLD,status,ierr)
                   end if
                   
                   call synchronize(71)
                   
                end if
             end if
          end do
       end do
       
       call send_message("Exchange the BC geometric successfully!")
       
    end subroutine exchange_bc_geometric      

    subroutine exchange_bc
       use global_variables                                                                        
       implicit none
       integer :: nb,pnb,i,j,k,m,ierr
       integer :: nr,nregions,ibctype
       integer :: iseq,nbt,nrt,packsize
       integer :: id_src,id_des,tag_seq
       integer :: ibeg(3),iend(3),idir,inrout,idir2
       integer,parameter :: mcyc(5) = (/1,2,3,1,2/)
       integer :: status(MPI_STATUS_SIZE)
       
       iseq = 0
       !!do pnb=1,pnblocks
       !!   nb = pnbindexs(pnb)
       do nb=1,nblocks
          id_src = mb_pids(nb)
          nregions = mb_bc(nb)%nregions
          do nr=1,nregions
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if (ibctype < 0) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
                if (id_src /= id_des) then
                   iseq = iseq + 1
                   tag_seq = iseq
                   nrt = mb_bc(nb)%bc(nr)%ibcwin
                   
                   if (myid == id_src-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
!!                      ibeg(idir) = ibeg(idir) - 2*max(inrout,0)
!!                      iend(idir) = iend(idir) - 2*min(inrout,0)
                      ibeg(idir) = ibeg(idir) - 4*max(inrout,0)
                      iend(idir) = iend(idir) - 4*min(inrout,0)
                      do m=1,2
                         idir2 = mcyc(idir+m)
                         if ( ibeg(idir2) > 1 ) then
                            ibeg(idir2) = ibeg(idir2) - 1
                         end if
                         if ( iend(idir2) < mb_dim(nb,idir2) ) then
                            iend(idir2) = iend(idir2) + 1
                         end if
                      end do
                      
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1)
                         mb_bc(nb)%bc(nr)%qpv(i,j,k,1) = mb_r(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv(i,j,k,2) = mb_u(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv(i,j,k,3) = mb_v(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv(i,j,k,4) = mb_w(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv(i,j,k,5) = mb_p(nb)%a3d(i,j,k)
                      end do
                      end do
                      end do
                      
                      if (nlamtur >= 0 ) then
                         do k=ibeg(3),iend(3)
                         do j=ibeg(2),iend(2)
                         do i=ibeg(1),iend(1)
                            mb_bc(nb)%bc(nr)%qpv(i,j,k,6) = mb_vist(nb)%a3d(i,j,k)
                            do m=7,nvtot
                               mb_bc(nb)%bc(nr)%qpv(i,j,k,m) = mb_qke(nb)%a4d(i,j,k,m-6)
                            end do
                         end do
                         end do
                         end do
                      end if
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*nvtot
                      
                      call MPI_SEND(mb_bc(nb)%bc(nr)%qpv,packsize,mpi_reprec, &
                                    id_des-1,tag_seq,MPI_COMM_WORLD,ierr)
                   end if
                   
                   if (myid == id_des-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
!!                      ibeg(idir) = ibeg(idir) - 2*max(inrout,0)
!!                      iend(idir) = iend(idir) - 2*min(inrout,0)
                      ibeg(idir) = ibeg(idir) - 4*max(inrout,0)
                      iend(idir) = iend(idir) - 4*min(inrout,0)
                      do m=1,2
                         idir2 = mcyc(idir+m)
                         if ( ibeg(idir2) > 1 ) then
                            ibeg(idir2) = ibeg(idir2) - 1
                         end if
                         if ( iend(idir2) < mb_dim(nb,idir2) ) then
                            iend(idir2) = iend(idir2) + 1
                         end if
                      end do
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*nvtot
                      
                      call MPI_RECV(mb_bc(nbt)%bc(nrt)%qpvpack,packsize,mpi_reprec, &
                                    id_src-1,tag_seq,MPI_COMM_WORLD,status,ierr)
                   end if
                   
!!                   call synchronize(81)
                   
                end if
             end if
          end do
       end do

!!       call send_message("Exchange the BC datas successfully!")
       
    end subroutine exchange_bc
    
    subroutine boundary_n1_other(nb,nr,bctype) !�Խӱ߽�����
    !*TGH. ���ʱ���Ȱ�ԭʼ���������ȡ��Ӧ�ص�����ֵ����������mb_qkeҲȡ��Ӧ�ص����ϵ�ֵ��
    !*TGH. Ȼ��call dif_average ��������call BC_face_dif_turbulence�����߽��ϵ�ֵ���ڱ߽������ֵ�ļ�ƽ��
    !*TGH. �Ѿ��޸ĵ��ʺ������߽����������ǻ����ʺϻ�ѧ��Ӧ�����
    
       use global_variables
       implicit none
       integer :: nb,nr,bctype,m,n
       integer :: nbs,s_nd,s_fix,s_lr
       integer :: nbt,t_nd,t_fix,t_lr
       integer :: is,js,ks,i,j,k,nt,ntarg
       integer :: it,jt,kt,ist,jst,kst
       integer :: s_st(3),s_ed(3)
       integer :: nv_s,nv_t,nfile_par_s,nfile_par_t,num_var_s,num_var_t,nchem_s,nchem_t
       real    :: prim_s1(nl),q_1(nl),zf1,ttt1
    
       do m=1,3
          s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
          s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
       enddo
       s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
       s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
       s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
       nbs   = mb_bc(nb)%bc(nr)%nbs             !���
    
       nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
       t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
       t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
       t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)
    
       if(cic1  /= 1)then !*TGH. һ��Խ�ʱ���ȸ�����ϵ�ֵ���ּ���߽��ϵ�ֵ
          do i = s_st(1),s_ed(1)
          do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
               is = i + mb_bc(nb)%bc(nr)%s_lr3d(1)
               js = j + mb_bc(nb)%bc(nr)%s_lr3d(2)
               ks = k + mb_bc(nb)%bc(nr)%s_lr3d(3)
               ist = i - mb_bc(nb)%bc(nr)%s_lr3d(1)
               jst = j - mb_bc(nb)%bc(nr)%s_lr3d(2)
               kst = k - mb_bc(nb)%bc(nr)%s_lr3d(3)
               it = mb_bc(nb)%bc(nr)%image(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(1)
               jt = mb_bc(nb)%bc(nr)%jmage(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(2)
               kt = mb_bc(nb)%bc(nr)%kmage(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(3)
!!               mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)
!!               mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
!!               mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
!!               mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
!!               mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)

               mb_r(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,1)
               mb_u(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,2)
               mb_v(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,3)
               mb_w(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,4)
               mb_p(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,5)
               
               mb_r(nbs)%a3d(i,j,k) = 0.5*(mb_r(nbs)%a3d(is,js,ks) + mb_r(nbs)%a3d(ist,jst,kst))
               mb_u(nbs)%a3d(i,j,k) = 0.5*(mb_u(nbs)%a3d(is,js,ks) + mb_u(nbs)%a3d(ist,jst,kst))
               mb_v(nbs)%a3d(i,j,k) = 0.5*(mb_v(nbs)%a3d(is,js,ks) + mb_v(nbs)%a3d(ist,jst,kst))
               mb_w(nbs)%a3d(i,j,k) = 0.5*(mb_w(nbs)%a3d(is,js,ks) + mb_w(nbs)%a3d(ist,jst,kst))
               mb_p(nbs)%a3d(i,j,k) = 0.5*(mb_p(nbs)%a3d(is,js,ks) + mb_p(nbs)%a3d(ist,jst,kst))
               

               if (nlamtur >= 0 ) then
                  mb_vist(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,6)
                  
                  mb_vist(nbs)%a3d(i,j,k) = 0.5*(mb_vist(nbs)%a3d(is,js,ks) + mb_vist(nbs)%a3d(ist,jst,kst))
                  
                  do m=7,nvtot
                     n = m-6
                     mb_qke(nbs)%a4d(is,js,ks,n) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,m)
                     
                     mb_qke(nbs)%a4d(i,j,k,n) = 0.5*(mb_qke(nbs)%a4d(is,js,ks,n) + mb_qke(nbs)%a4d(ist,jst,kst,n))
                  end do
               end if
               
!!             call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC
          enddo
          enddo
          enddo
          
!!          call dif_average(nb,nr,bctype,s_st,s_ed) !* TGH. �߽��ϵ�ֵ���������ƽ����ͬʱ���д�������ģʽ�߽� *!
    
       else !*TGH. �����Խ�ʱ��ֻ������ϵ�ֵ��������߽��ϵ�ֵ��������ģ������
    
          do i = s_st(1),s_ed(1)
          do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
               is = i + mb_bc(nb)%bc(nr)%s_lr3d(1)
               js = j + mb_bc(nb)%bc(nr)%s_lr3d(2)
               ks = k + mb_bc(nb)%bc(nr)%s_lr3d(3)
               ist = i - mb_bc(nb)%bc(nr)%s_lr3d(1)
               jst = j - mb_bc(nb)%bc(nr)%s_lr3d(2)
               kst = k - mb_bc(nb)%bc(nr)%s_lr3d(3)
               it = mb_bc(nb)%bc(nr)%image(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(1)
               jt = mb_bc(nb)%bc(nr)%jmage(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(2)
               kt = mb_bc(nb)%bc(nr)%kmage(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(3)
!!               mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)
!!               mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
!!               mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
!!               mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
!!               mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)
               
               mb_r(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,1)
               mb_u(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,2)
               mb_v(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,3)
               mb_w(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,4)
               mb_p(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,5)
               
!!               mb_r(nbs)%a3d(i,j,k) = 0.5*(mb_r(nbs)%a3d(is,js,ks) + mb_r(nbs)%a3d(ist,jst,kst))
!!               mb_u(nbs)%a3d(i,j,k) = 0.5*(mb_u(nbs)%a3d(is,js,ks) + mb_u(nbs)%a3d(ist,jst,kst))
!!               mb_v(nbs)%a3d(i,j,k) = 0.5*(mb_v(nbs)%a3d(is,js,ks) + mb_v(nbs)%a3d(ist,jst,kst))
!!               mb_w(nbs)%a3d(i,j,k) = 0.5*(mb_w(nbs)%a3d(is,js,ks) + mb_w(nbs)%a3d(ist,jst,kst))
!!               mb_p(nbs)%a3d(i,j,k) = 0.5*(mb_p(nbs)%a3d(is,js,ks) + mb_p(nbs)%a3d(ist,jst,kst))
               

               if (nlamtur >= 0 ) then
                  mb_vist(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,6)
                  
                  mb_vist(nbs)%a3d(i,j,k) = 0.5*(mb_vist(nbs)%a3d(is,js,ks) + mb_vist(nbs)%a3d(ist,jst,kst))
                  
                  do m=7,nvtot
                     n = m-6
                     mb_qke(nbs)%a4d(is,js,ks,n) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,m)
                     
                     mb_qke(nbs)%a4d(i,j,k,n) = 0.5*(mb_qke(nbs)%a4d(is,js,ks,n) + mb_qke(nbs)%a4d(ist,jst,kst,n))
                  end do
               end if
                  
!!               call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC
!!               call BC_face_dif_turbulence(nb,nb,is,js,ks,ist,jst,kst,i,j,k)
    
          enddo
          enddo
          enddo
    
       endif
    
    end subroutine boundary_n1_other

    subroutine boundary_n1_parallel(nb,nr,bctype)
       use global_variables
       implicit none
       integer :: nb,nr,bctype
       integer :: nbt,id_src,id_des
       
       nbt = mb_bc(nb)%bc(nr)%nbt
       id_src = mb_pids(nb)
       id_des = mb_pids(nbt)
       if (id_src == id_des) then
          call boundary_n1(nb,nr,bctype)    
       else
          call boundary_n1_other(nb,nr,bctype)    
       end if
       
    end subroutine boundary_n1_parallel


    subroutine boundary_n1_vir1_other(nb,nr,bctype) !�Խӱ߽�����
    !_____________________________________________________________________!
    !* 1�׾��� ������ ��� �Խ����ϵ��������ϵ�ԭʼ����
    !* ͨ��ԭʼ�����ڵѿ�������ϵ�еĵ�����������ϵ�ֵ
    !* ����mb_qke���ü򵥶Խӵķ���ֱ��ȡ���ڿ��ϵ�ֵ��
    !*     Ȼ��call dif_average ��������call BC_face_dif_turbulence����
    !*     �߽��ϵ�ֵ���ڱ߽������ֵ�ļ�ƽ��
    !* �����Խ�ʱֻ������ϵ�ֵ��һ��Խ�ʱ��Ҫ��Խ����ϵ�ֵ�����ǻ����ʺϻ�ѧ��Ӧ�����
    !* ��ƣ�ëö��
    !* ���ԣ�Ϳ����
    !_____________________________________________________________________!
    	use define_precision_mod
        use global_variables,only : mb_r,mb_u,mb_v,mb_w,mb_p,mb_bc,method,cic1,mb_vol &
    	                          , mb_kcx,mb_kcy,mb_kcz,mb_etx,mb_ety,mb_etz,mb_ctx,mb_cty,mb_ctz &
    							  , mb_dim,nl,nmax,mb_x,mb_y,mb_z,dis_tao,nlamtur,mb_vist,mb_qke
        implicit none
        integer :: nb,nr,bctype,m,n
        integer :: nbs,s_nd,s_fix,s_lr
        integer :: nbt,t_nd,t_fix,t_lr
        integer :: is,js,ks,i,j,k,nt,ntarg
        integer :: it,jt,kt,ist,jst,kst,it0,jt0,kt0
        integer :: s_st(3),s_ed(3),s_lr3d(3)
        real(prec) :: prim_s1(nl),q_1(nl),zf1,ttt1
    	real(prec) :: d1i,d1j,d1k,d1x,d1y,d1z
    	real(prec) :: volp,nxyz(9),nijk(9),druvwp(15),druvwpdxyz(15)
    	real(prec) :: tk1,tk2,tk3
    	real(prec) :: rm,um,vm,wm,pm
    	real(prec) :: dis1,dis2,dis12 ! �Խ�����������
    	integer :: nit,njt,nkt,nis,njs,nks
    	integer :: nerr,key_scm
    
    !	dis_tao = 5.0 !���ȳ���dis_taoʱ��������㷨��õĲ���
        key_scm = 1
    
        do m=1,3
           s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
           s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
        enddo
        s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
        s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
        s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
        nbs   = mb_bc(nb)%bc(nr)%nbs             !���
    
        nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
        t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
        t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
        t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)
    
        nis = mb_dim(nb,1)
        njs = mb_dim(nb,2)
        nks = mb_dim(nb,3)
    
        nit = mb_dim(nbt,1)
        njt = mb_dim(nbt,2)
        nkt = mb_dim(nbt,3)
    
    	s_lr3d(1) = mb_bc(nb)%bc(nr)%s_lr3d(1)
    	s_lr3d(2) = mb_bc(nb)%bc(nr)%s_lr3d(2)
    	s_lr3d(3) = mb_bc(nb)%bc(nr)%s_lr3d(3)
    	
    	   nerr = 0
           do i = s_st(1),s_ed(1)
           do j = s_st(2),s_ed(2)
           do k = s_st(3),s_ed(3)
    
    
                is = i + s_lr3d(1)
                js = j + s_lr3d(2)
                ks = k + s_lr3d(3)
                ist = i - s_lr3d(1)
                jst = j - s_lr3d(2)
                kst = k - s_lr3d(3)
    
    			it0 = mb_bc(nb)%bc(nr)%image(i,j,k)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,k)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,k)
                it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)
                jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)
                kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)
    
    			dis1 = sqrt( (mb_x(nbs)%a3d(i , j , k ) - mb_x(nbs)%a3d(ist,jst,kst))**2 + &
    						 (mb_y(nbs)%a3d(i , j , k ) - mb_y(nbs)%a3d(ist,jst,kst))**2 + &
    						 (mb_z(nbs)%a3d(i , j , k ) - mb_z(nbs)%a3d(ist,jst,kst))**2 )
    
!!    			dis2 = sqrt( (mb_x(nbt)%a3d(it0,jt0,kt0) - mb_x(nbt)%a3d(it ,jt ,kt ))**2 + &
!!    						 (mb_y(nbt)%a3d(it0,jt0,kt0) - mb_y(nbt)%a3d(it ,jt ,kt ))**2 + &
!!    						 (mb_z(nbt)%a3d(it0,jt0,kt0) - mb_z(nbt)%a3d(it ,jt ,kt ))**2 )
    			dis2 = mb_bc(nb)%bc(nr)%dispack(it0,jt0,kt0)
    
    			dis12 = dis1/dis2
    
                !----------------------------------------------------------------------!
                !* ���¿�ʼ�����Ӧ���ڶԽ�����ԭʼ������ֱ������ϵ�е�һ�׵���        !
!!    			volp    = mb_vol(nbt)%a3d(it0,jt0,kt0)
!!    			nxyz(1) = mb_kcx(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(2) = mb_kcy(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(3) = mb_kcz(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(4) = mb_etx(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(5) = mb_ety(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(6) = mb_etz(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(7) = mb_ctx(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(8) = mb_cty(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(9) = mb_ctz(nbt)%a3d(it0,jt0,kt0) / volp
    			volp    = mb_bc(nb)%bc(nr)%volpack(it0,jt0,kt0)
    			nxyz(1:9) = mb_bc(nb)%bc(nr)%sxyzpack(it0,jt0,kt0,1:9) / volp
    
    
    			if(it0 == 1) then
    
    			  if(key_scm == 1)then
!!    			  druvwp(1) =  mb_r(nbt)%a3d(2,jt0,kt0) &
!!    			              -mb_r(nbt)%a3d(1,jt0,kt0)
!!    			  druvwp(2) =  mb_u(nbt)%a3d(2,jt0,kt0) &
!!    			              -mb_u(nbt)%a3d(1,jt0,kt0)
!!    			  druvwp(3) =  mb_v(nbt)%a3d(2,jt0,kt0) &
!!    			              -mb_v(nbt)%a3d(1,jt0,kt0)
!!    			  druvwp(4) =  mb_w(nbt)%a3d(2,jt0,kt0) &
!!    			              -mb_w(nbt)%a3d(1,jt0,kt0)
!!    			  druvwp(5) =  mb_p(nbt)%a3d(2,jt0,kt0) &
!!    			              -mb_p(nbt)%a3d(1,jt0,kt0)
    			  druvwp(1:5) =  mb_bc(nb)%bc(nr)%qpvpack(2,jt0,kt0,1:5) &
    			                -mb_bc(nb)%bc(nr)%qpvpack(1,jt0,kt0,1:5)
    			  else
    
!!    			  druvwp(1) = -1.5_prec*mb_r(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_r(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_r(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(2) = -1.5_prec*mb_u(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_u(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_u(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(3) = -1.5_prec*mb_v(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_v(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_v(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(4) = -1.5_prec*mb_w(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_w(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_w(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(5) = -1.5_prec*mb_p(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_p(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_p(nbt)%a3d(3,jt0,kt0)
    			  druvwp(1:5) = -1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(1,jt0,kt0,1:5) &
    			                +2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(2,jt0,kt0,1:5) &
    			                -0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(3,jt0,kt0,1:5)
    			  endif
    
    			elseif(it0 == nit)then
    
    			  if(key_scm == 1)then
    
!!    			  druvwp(1) =  mb_r(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -mb_r(nbt)%a3d(it0-1,jt0,kt0)
!!    			  druvwp(2) =  mb_u(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -mb_u(nbt)%a3d(it0-1,jt0,kt0)
!!    			  druvwp(3) =  mb_v(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -mb_v(nbt)%a3d(it0-1,jt0,kt0)
!!    			  druvwp(4) =  mb_w(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -mb_w(nbt)%a3d(it0-1,jt0,kt0)
!!    			  druvwp(5) =  mb_p(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -mb_p(nbt)%a3d(it0-1,jt0,kt0)
    			  druvwp(1:5) =  mb_bc(nb)%bc(nr)%qpvpack(it0  ,jt0,kt0,1:5) &
    			                -mb_bc(nb)%bc(nr)%qpvpack(it0-1,jt0,kt0,1:5)
    			  else
    
!!    			  druvwp(1) =  1.5_prec*mb_r(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_r(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_r(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(2) =  1.5_prec*mb_u(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_u(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_u(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(3) =  1.5_prec*mb_v(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_v(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_v(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(4) =  1.5_prec*mb_w(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_w(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_w(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(5) =  1.5_prec*mb_p(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_p(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_p(nbt)%a3d(it0-2,jt0,kt0)
    			  druvwp(1:5) =  1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0  ,jt0,kt0,1:5) &
    			                -2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0-1,jt0,kt0,1:5) &
    			                +0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0-2,jt0,kt0,1:5)
    			  endif
    
    			else
!!    			  druvwp(1) =  0.5_prec*( mb_r(nbt)%a3d(it0+1,jt0,kt0)-mb_r(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(2) =  0.5_prec*( mb_u(nbt)%a3d(it0+1,jt0,kt0)-mb_u(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(3) =  0.5_prec*( mb_v(nbt)%a3d(it0+1,jt0,kt0)-mb_v(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(4) =  0.5_prec*( mb_w(nbt)%a3d(it0+1,jt0,kt0)-mb_w(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(5) =  0.5_prec*( mb_p(nbt)%a3d(it0+1,jt0,kt0)-mb_p(nbt)%a3d(it0-1,jt0,kt0) )
    			  druvwp(1:5) =  0.5_prec*( mb_bc(nb)%bc(nr)%qpvpack(it0+1,jt0,kt0,1:5) &
    			                           -mb_bc(nb)%bc(nr)%qpvpack(it0-1,jt0,kt0,1:5) )
    			endif
    
    			if(jt0 == 1)then
    
      			  if(key_scm == 1)then
    
!!    			  druvwp(6) =  mb_r(nbt)%a3d(it0,2,kt0) &
!!    			              -mb_r(nbt)%a3d(it0,1,kt0)
!!    			  druvwp(7) =  mb_u(nbt)%a3d(it0,2,kt0) &
!!    			              -mb_u(nbt)%a3d(it0,1,kt0)
!!    			  druvwp(8) =  mb_v(nbt)%a3d(it0,2,kt0) &
!!    			              -mb_v(nbt)%a3d(it0,1,kt0)
!!    			  druvwp(9) =  mb_w(nbt)%a3d(it0,2,kt0) &
!!    			              -mb_w(nbt)%a3d(it0,1,kt0)
!!    			  druvwp(10)=  mb_p(nbt)%a3d(it0,2,kt0) &
!!    			              -mb_p(nbt)%a3d(it0,1,kt0)
    			  druvwp(6:10) =  mb_bc(nb)%bc(nr)%qpvpack(it0,2,kt0,1:5) &
    			                 -mb_bc(nb)%bc(nr)%qpvpack(it0,1,kt0,1:5)
    
    			  else
    
!!    			  druvwp(6) = -1.5_prec*mb_r(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_r(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_r(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(7) = -1.5_prec*mb_u(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_u(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_u(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(8) = -1.5_prec*mb_v(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_v(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_v(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(9) = -1.5_prec*mb_w(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_w(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_w(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(10)= -1.5_prec*mb_p(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_p(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_p(nbt)%a3d(it0,3,kt0)
    			  druvwp(6:10) = -1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,1,kt0,1:5) &
    			                 +2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,2,kt0,1:5) &
    			                 -0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,3,kt0,1:5)
    
    			  endif
    
    			elseif(jt0 == njt)then
    
    			  if(key_scm == 1)then
!!    			  druvwp(6) =  mb_r(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -mb_r(nbt)%a3d(it0,jt0-1,kt0)
!!    			  druvwp(7) =  mb_u(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -mb_u(nbt)%a3d(it0,jt0-1,kt0)
!!    			  druvwp(8) =  mb_v(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -mb_v(nbt)%a3d(it0,jt0-1,kt0)
!!    			  druvwp(9) =  mb_w(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -mb_w(nbt)%a3d(it0,jt0-1,kt0)
!!    			  druvwp(10)=  mb_p(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -mb_p(nbt)%a3d(it0,jt0-1,kt0)
    			  druvwp(6:10) =  mb_bc(nb)%bc(nr)%qpvpack(it0,jt0  ,kt0,1:5) &
    			                 -mb_bc(nb)%bc(nr)%qpvpack(it0,jt0-1,kt0,1:5)
    
    			  else
    
!!    			  druvwp(6) =  1.5_prec*mb_r(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(7) =  1.5_prec*mb_u(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(8) =  1.5_prec*mb_v(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(9) =  1.5_prec*mb_w(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(10)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0-2,kt0)
    			  druvwp(6:10) =  1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0  ,kt0,1:5) &
    			                 -2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0-1,kt0,1:5) &
    			                 +0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0-2,kt0,1:5)
    			  endif
    
    			else
!!    			  druvwp(6) =  0.5_prec*( mb_r(nbt)%a3d(it0,jt0+1,kt0)-mb_r(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(7) =  0.5_prec*( mb_u(nbt)%a3d(it0,jt0+1,kt0)-mb_u(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(8) =  0.5_prec*( mb_v(nbt)%a3d(it0,jt0+1,kt0)-mb_v(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(9) =  0.5_prec*( mb_w(nbt)%a3d(it0,jt0+1,kt0)-mb_w(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(10)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0+1,kt0)-mb_p(nbt)%a3d(it0,jt0-1,kt0) )
    			  druvwp(6:10)=  0.5_prec*( mb_bc(nb)%bc(nr)%qpvpack(it0,jt0+1,kt0,1:5) &
    			                           -mb_bc(nb)%bc(nr)%qpvpack(it0,jt0-1,kt0,1:5) )
    			endif
    
    			if(kt0 == 1)then
    
    			  if(key_scm == 1)then
    
!!    			  druvwp(11)=  mb_r(nbt)%a3d(it0,jt0,2) &
!!    			              -mb_r(nbt)%a3d(it0,jt0,1)
!!    			  druvwp(12)=  mb_u(nbt)%a3d(it0,jt0,2) &
!!    			              -mb_u(nbt)%a3d(it0,jt0,1)
!!    			  druvwp(13)=  mb_v(nbt)%a3d(it0,jt0,2) &
!!    			              -mb_v(nbt)%a3d(it0,jt0,1)
!!    			  druvwp(14)=  mb_w(nbt)%a3d(it0,jt0,2) &
!!    			              -mb_w(nbt)%a3d(it0,jt0,1)
!!    			  druvwp(15)=  mb_p(nbt)%a3d(it0,jt0,2) &
!!    			              -mb_p(nbt)%a3d(it0,jt0,1)
    			  druvwp(11:15)=  mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,2,1:5) &
    			                 -mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,1,1:5)
    
    			  else
    
!!    			  druvwp(11)= -1.5_prec*mb_r(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_r(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_r(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(12)= -1.5_prec*mb_u(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_u(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_u(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(13)= -1.5_prec*mb_v(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_v(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_v(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(14)= -1.5_prec*mb_w(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_w(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_w(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(15)= -1.5_prec*mb_p(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_p(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_p(nbt)%a3d(it0,jt0,3)
    			  druvwp(11:15)= -1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,1,1:5) &
    			                 +2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,2,1:5) &
    			                 -0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,3,1:5)
    			  endif
    
    			elseif(kt0 == nkt)then
    
    			  if(key_scm == 1)then
    
!!    			  druvwp(11)=  mb_r(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -mb_r(nbt)%a3d(it0,jt0,kt0-1)
!!    			  druvwp(12)=  mb_u(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -mb_u(nbt)%a3d(it0,jt0,kt0-1)
!!    			  druvwp(13)=  mb_v(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -mb_v(nbt)%a3d(it0,jt0,kt0-1)
!!    			  druvwp(14)=  mb_w(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -mb_w(nbt)%a3d(it0,jt0,kt0-1)
!!    			  druvwp(15)=  mb_p(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -mb_p(nbt)%a3d(it0,jt0,kt0-1)
    			  druvwp(11:15)=  mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0  ,1:5) &
    			                 -mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0-1,1:5)
    
    			  else
    
!!    			  druvwp(11)=  1.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(12)=  1.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(13)=  1.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(14)=  1.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(15)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0-2)
    			  druvwp(11:15)=  1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0  ,1:5) &
    			                 -2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0-1,1:5) &
    			                 +0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0-2,1:5)
    
    			  endif
    
    			else
!!    			  druvwp(11)=  0.5_prec*( mb_r(nbt)%a3d(it0,jt0,kt0+1)-mb_r(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(12)=  0.5_prec*( mb_u(nbt)%a3d(it0,jt0,kt0+1)-mb_u(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(13)=  0.5_prec*( mb_v(nbt)%a3d(it0,jt0,kt0+1)-mb_v(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(14)=  0.5_prec*( mb_w(nbt)%a3d(it0,jt0,kt0+1)-mb_w(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(15)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0,kt0+1)-mb_p(nbt)%a3d(it0,jt0,kt0-1) )
    			  druvwp(11:15)=  0.5_prec*( mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0+1,1:5) &
    			                            -mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0-1,1:5) )
    			endif
    
    			call druvwp_xyz(nxyz,druvwp,druvwpdxyz)
                     !* �Ѽ���ϵ��ԭʼ������һ�׵���ת����ֱ��ϵ��һ�׵���
                     !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                     !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk
                     !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz
    
                !------------------------------------------------------------------------------------!
    			!* ���°ѶԽ�����ԭʼ������ֱ������ϵ�е�һ�׵��������ڱ��������м���ϵ�µ�һ�׵���  !
    
    			volp    = mb_vol(nb)%a3d(i,j,k)
    			nxyz(1) = mb_kcx(nb)%a3d(i,j,k) / volp
    			nxyz(2) = mb_kcy(nb)%a3d(i,j,k) / volp
    			nxyz(3) = mb_kcz(nb)%a3d(i,j,k) / volp
    			nxyz(4) = mb_etx(nb)%a3d(i,j,k) / volp
    			nxyz(5) = mb_ety(nb)%a3d(i,j,k) / volp
    			nxyz(6) = mb_etz(nb)%a3d(i,j,k) / volp
    			nxyz(7) = mb_ctx(nb)%a3d(i,j,k) / volp
    			nxyz(8) = mb_cty(nb)%a3d(i,j,k) / volp
    			nxyz(9) = mb_ctz(nb)%a3d(i,j,k) / volp
    
    !			call get_nijk(nb,nis,njs,nks,i,j,k,nijk) !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
    			call nxyz_nijk(nxyz,nijk,volp) !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                                               !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
    			call druvwp_xyz(nijk,druvwpdxyz,druvwp)
                     !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz
    			     !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk
    
               !---------------------------------------------------------------------------------!
                !* ���¸����ݶ��� �����߽��֮��Ĳ���                 ------------------------!
    			rm =  s_lr3d(1)*druvwp(1) + s_lr3d(2)*druvwp(6)  + s_lr3d(3)*druvwp(11)
    			um =  s_lr3d(1)*druvwp(2) + s_lr3d(2)*druvwp(7)  + s_lr3d(3)*druvwp(12)
    			vm =  s_lr3d(1)*druvwp(3) + s_lr3d(2)*druvwp(8)  + s_lr3d(3)*druvwp(13)
    			wm =  s_lr3d(1)*druvwp(4) + s_lr3d(2)*druvwp(9)  + s_lr3d(3)*druvwp(14)
    			pm =  s_lr3d(1)*druvwp(5) + s_lr3d(2)*druvwp(10) + s_lr3d(3)*druvwp(15)
    
    			if(dis12 > dis_tao)then
    				dis12 = dis_tao/dis12
    
    				rm = rm * dis12
    				um = um * dis12
    				vm = vm * dis12
    				wm = wm * dis12
    				pm = pm * dis12
    			endif
    
              !* ���¸��� "����" ������ϵ�ֵ                               ------------------------!
    			rm = mb_r(nb)%a3d(i,j,k) + rm
    			um = mb_u(nb)%a3d(i,j,k) + um
    			vm = mb_v(nb)%a3d(i,j,k) + vm
    			wm = mb_w(nb)%a3d(i,j,k) + wm
    			pm = mb_p(nb)%a3d(i,j,k) + pm
    
901	format(1x,'��',i3,' ��ĵ�',i2,' �߽���õ����������ֵʧ�ܣ��ܶȣ�',e12.4,3i4)
902	format(1x,'��',i3,' ��ĵ�',i2,' �߽���õ����������ֵʧ�ܣ�ѹ����',e12.4,3i4)
    
     			if(pm <= 0._prec .or. rm <= 0._prec) then
    			    nerr = nerr + 1
    !				if(nerr<3)write(*,902)nb,nr,pm,i,j,k
    
!!                    rm = mb_r(nbt)%a3d(it,jt,kt)
!!                    um = mb_u(nbt)%a3d(it,jt,kt)
!!                    vm = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. ����ԭʼ���� �� ��Ӧ������������һ����ֵ
!!                    wm = mb_w(nbt)%a3d(it,jt,kt)
!!                    pm = mb_p(nbt)%a3d(it,jt,kt)
                    rm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,1)
                    um = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,2)
                    vm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,3)  !*TGH. ����ԭʼ���� �� ��Ӧ������������һ����ֵ
                    wm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,4)
                    pm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,5)
    
    			endif
    
                mb_r(nb)%a3d(is,js,ks) = rm
                mb_u(nb)%a3d(is,js,ks) = um
                mb_v(nb)%a3d(is,js,ks) = vm
                mb_w(nb)%a3d(is,js,ks) = wm
                mb_p(nb)%a3d(is,js,ks) = pm
               !* ��㸳ֵ���  ------------------------!
               
    
!!    			call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC
!!    			do m=1,cic1
!!    			    call BC_face_dif_turbulence(nb,nb,is,js,ks,ist,jst,kst,i,j,k)
!!    			enddo
               if (cic1 /= 1) then
                  mb_r(nbs)%a3d(i,j,k) = 0.5*(mb_r(nbs)%a3d(is,js,ks) + mb_r(nbs)%a3d(ist,jst,kst))
                  mb_u(nbs)%a3d(i,j,k) = 0.5*(mb_u(nbs)%a3d(is,js,ks) + mb_u(nbs)%a3d(ist,jst,kst))
                  mb_v(nbs)%a3d(i,j,k) = 0.5*(mb_v(nbs)%a3d(is,js,ks) + mb_v(nbs)%a3d(ist,jst,kst))
                  mb_w(nbs)%a3d(i,j,k) = 0.5*(mb_w(nbs)%a3d(is,js,ks) + mb_w(nbs)%a3d(ist,jst,kst))
                  mb_p(nbs)%a3d(i,j,k) = 0.5*(mb_p(nbs)%a3d(is,js,ks) + mb_p(nbs)%a3d(ist,jst,kst))
               end if
               if (nlamtur >= 0 ) then
                  mb_vist(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,6)
                  
                  mb_vist(nbs)%a3d(i,j,k) = 0.5*(mb_vist(nbs)%a3d(is,js,ks) + mb_vist(nbs)%a3d(ist,jst,kst))
                
                  do m=7,nvtot
                     n = m-6
                     mb_qke(nbs)%a4d(is,js,ks,n) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,m)
                     
                     mb_qke(nbs)%a4d(i,j,k,n) = 0.5*(mb_qke(nbs)%a4d(is,js,ks,n) + mb_qke(nbs)%a4d(ist,jst,kst,n))
                  end do
               end if
           enddo
           enddo
           enddo
    
903	format(1x,'��',i3,' ��ĵ�',i2,'�߽繲', i4 ' ����õ����������ֵʧ��')
    
    	   if(nerr > 0)write(*,903)nb,nr,nerr
    
!!    	   if(cic1 /= 1) call dif_average(nb,nr,bctype,s_st,s_ed) !* TGH. �߽��ϵ�ֵ���������ƽ����ͬʱ���д�������ģʽ�߽� *!
    
    
        return
    end subroutine boundary_n1_vir1_other


    subroutine boundary_n1_vir1_parallel(nb,nr,bctype)
       use global_variables
       implicit none
       integer :: nb,nr,bctype
       integer :: nbt,id_src,id_des
       
       nbt = mb_bc(nb)%bc(nr)%nbt
       id_src = mb_pids(nb)
       id_des = mb_pids(nbt)
       if (id_src == id_des) then
          call boundary_n1_vir1(nb,nr,bctype)    
       else
          call boundary_n1_vir1_other(nb,nr,bctype)    
       end if
       
    end subroutine boundary_n1_vir1_parallel

    subroutine boundary_n1_vir2_other(nb,nr,bctype) !�Խӱ߽�����
    !_____________________________________________________________________!
    !* 2�׾��� ������ ��� �Խ����ϵ��������ϵ�ԭʼ����
    !* ͨ��ԭʼ�����ڵѿ�������ϵ�еĵ�����������ϵ�ֵ
    !* ����mb_qke���ü򵥶Խӵķ���ֱ��ȡ���ڿ��ϵ�ֵ��
    !*     Ȼ��call dif_average ��������call BC_face_dif_turbulence����
    !*     �߽��ϵ�ֵ���ڱ߽������ֵ�ļ�ƽ��
    !* �����Խ�ʱֻ������ϵ�ֵ��һ��Խ�ʱ��Ҫ��Խ����ϵ�ֵ�����ǻ����ʺϻ�ѧ��Ӧ�����
    !* ��ƣ�ëö��
    !* ���ԣ�Ϳ����
    !_____________________________________________________________________!
    	use define_precision_mod
        use global_variables,only : mb_r,mb_u,mb_v,mb_w,mb_p,mb_bc,method,cic1,mb_vol &
    	                          , mb_kcx,mb_kcy,mb_kcz,mb_etx,mb_ety,mb_etz,mb_ctx,mb_cty,mb_ctz &
    							  , mb_dim,nl,nmax,mb_x,mb_y,mb_z,dis_tao,nlamtur,mb_vist,mb_qke
        implicit none
        integer :: nb,nr,bctype,m,n
        integer :: nbs,s_nd,s_fix,s_lr
        integer :: nbt,t_nd,t_fix,t_lr
        integer :: is,js,ks,i,j,k,nt,ntarg
        integer :: it,jt,kt,ist,jst,kst,it0,jt0,kt0
        integer :: s_st(3),s_ed(3),s_lr3d(3)
        real(prec) :: prim_s1(nl),q_1(nl),zf1,ttt1
    	real(prec) :: d1i,d1j,d1k,d1x,d1y,d1z
    	real(prec) :: volp,nxyz(9),nijk(9),druvwp(15),druvwpdxyz(15)
    	real(prec) :: tk1,tk2,tk3
    	real(prec) :: rm,um,vm,wm,pm
    	real(prec) :: dis1,dis2,dis12  ! �Խ�����������
    	integer :: nit,njt,nkt,nis,njs,nks
    	integer :: nerr
    
    !	dis_tao = 2.0 !���ȳ���dis_taoʱ��������㷨��õĲ���
        do m=1,3
           s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
           s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
        enddo
        s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
        s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
        s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
        nbs   = mb_bc(nb)%bc(nr)%nbs             !���
    
        nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
        t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
        t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
        t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)
    
        nis = mb_dim(nb,1)
        njs = mb_dim(nb,2)
        nks = mb_dim(nb,3)
    
        nit = mb_dim(nbt,1)
        njt = mb_dim(nbt,2)
        nkt = mb_dim(nbt,3)
    
    	s_lr3d(1) = mb_bc(nb)%bc(nr)%s_lr3d(1)
    	s_lr3d(2) = mb_bc(nb)%bc(nr)%s_lr3d(2)
    	s_lr3d(3) = mb_bc(nb)%bc(nr)%s_lr3d(3)
    
    	   nerr = 0
           do i = s_st(1),s_ed(1)
           do j = s_st(2),s_ed(2)
           do k = s_st(3),s_ed(3)
    
    
                is = i + s_lr3d(1)
                js = j + s_lr3d(2)
                ks = k + s_lr3d(3)
                ist = i - s_lr3d(1)
                jst = j - s_lr3d(2)
                kst = k - s_lr3d(3)
    
    			it0 = mb_bc(nb)%bc(nr)%image(i,j,k)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,k)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,k)
                it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)
                jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)
                kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)
    
    			dis1 = sqrt( (mb_x(nbs)%a3d(i , j , k ) - mb_x(nbs)%a3d(ist,jst,kst))**2 + &
    						 (mb_y(nbs)%a3d(i , j , k ) - mb_y(nbs)%a3d(ist,jst,kst))**2 + &
    						 (mb_z(nbs)%a3d(i , j , k ) - mb_z(nbs)%a3d(ist,jst,kst))**2 )
    
!!    			dis2 = sqrt( (mb_x(nbt)%a3d(it0,jt0,kt0) - mb_x(nbt)%a3d(it ,jt ,kt ))**2 + &
!!    						 (mb_y(nbt)%a3d(it0,jt0,kt0) - mb_y(nbt)%a3d(it ,jt ,kt ))**2 + &
!!    						 (mb_z(nbt)%a3d(it0,jt0,kt0) - mb_z(nbt)%a3d(it ,jt ,kt ))**2 )
    			dis2 = mb_bc(nb)%bc(nr)%dispack(it0,jt0,kt0)
    
    			dis12 = dis1/dis2
    
                !----------------------------------------------------------------------!
                !* ���¿�ʼ�����Ӧ���ڶԽ�����ԭʼ������ֱ������ϵ�е�һ�׵���        !
!!    			volp    = mb_vol(nbt)%a3d(it0,jt0,kt0)
!!    			nxyz(1) = mb_kcx(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(2) = mb_kcy(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(3) = mb_kcz(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(4) = mb_etx(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(5) = mb_ety(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(6) = mb_etz(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(7) = mb_ctx(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(8) = mb_cty(nbt)%a3d(it0,jt0,kt0) / volp
!!    			nxyz(9) = mb_ctz(nbt)%a3d(it0,jt0,kt0) / volp
    			volp    = mb_bc(nb)%bc(nr)%volpack(it0,jt0,kt0)
    			nxyz(1:9) = mb_bc(nb)%bc(nr)%sxyzpack(it0,jt0,kt0,1:9) / volp
    
    
    			if(it0 == 1) then
!!    			  druvwp(1) = -1.5_prec*mb_r(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_r(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_r(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(2) = -1.5_prec*mb_u(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_u(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_u(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(3) = -1.5_prec*mb_v(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_v(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_v(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(4) = -1.5_prec*mb_w(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_w(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_w(nbt)%a3d(3,jt0,kt0)
!!    			  druvwp(5) = -1.5_prec*mb_p(nbt)%a3d(1,jt0,kt0) &
!!    			              +2.0_prec*mb_p(nbt)%a3d(2,jt0,kt0) &
!!    			              -0.5_prec*mb_p(nbt)%a3d(3,jt0,kt0)
    			  druvwp(1:5) = -1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(1,jt0,kt0,1:5) &
    			                +2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(2,jt0,kt0,1:5) &
    			                -0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(3,jt0,kt0,1:5)
    			              
    			elseif(it0 == nit)then
!!    			  druvwp(1) =  1.5_prec*mb_r(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_r(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_r(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(2) =  1.5_prec*mb_u(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_u(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_u(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(3) =  1.5_prec*mb_v(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_v(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_v(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(4) =  1.5_prec*mb_w(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_w(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_w(nbt)%a3d(it0-2,jt0,kt0)
!!    			  druvwp(5) =  1.5_prec*mb_p(nbt)%a3d(it0  ,jt0,kt0) &
!!    			              -2.0_prec*mb_p(nbt)%a3d(it0-1,jt0,kt0) &
!!    			              +0.5_prec*mb_p(nbt)%a3d(it0-2,jt0,kt0)
    			  druvwp(1:5) =  1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0  ,jt0,kt0,1:5) &
    			                -2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0-1,jt0,kt0,1:5) &
    			                +0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0-2,jt0,kt0,1:5)
    			else
!!    			  druvwp(1) =  0.5_prec*( mb_r(nbt)%a3d(it0+1,jt0,kt0)-mb_r(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(2) =  0.5_prec*( mb_u(nbt)%a3d(it0+1,jt0,kt0)-mb_u(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(3) =  0.5_prec*( mb_v(nbt)%a3d(it0+1,jt0,kt0)-mb_v(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(4) =  0.5_prec*( mb_w(nbt)%a3d(it0+1,jt0,kt0)-mb_w(nbt)%a3d(it0-1,jt0,kt0) )
!!    			  druvwp(5) =  0.5_prec*( mb_p(nbt)%a3d(it0+1,jt0,kt0)-mb_p(nbt)%a3d(it0-1,jt0,kt0) )
    			  druvwp(1:5) =  0.5_prec*( mb_bc(nb)%bc(nr)%qpvpack(it0+1,jt0,kt0,1:5) &
    			                           -mb_bc(nb)%bc(nr)%qpvpack(it0-1,jt0,kt0,1:5) )
    			endif
    
    			if(jt0 == 1)then
!!    			  druvwp(6) = -1.5_prec*mb_r(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_r(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_r(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(7) = -1.5_prec*mb_u(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_u(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_u(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(8) = -1.5_prec*mb_v(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_v(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_v(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(9) = -1.5_prec*mb_w(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_w(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_w(nbt)%a3d(it0,3,kt0)
!!    			  druvwp(10)= -1.5_prec*mb_p(nbt)%a3d(it0,1,kt0) &
!!    			              +2.0_prec*mb_p(nbt)%a3d(it0,2,kt0) &
!!    			              -0.5_prec*mb_p(nbt)%a3d(it0,3,kt0)
    			  druvwp(6:10) = -1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,1,kt0,1:5) &
    			                 +2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,2,kt0,1:5) &
    			                 -0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,3,kt0,1:5)
    			              
    			elseif(jt0 == njt)then
!!    			  druvwp(6) =  1.5_prec*mb_r(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(7) =  1.5_prec*mb_u(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(8) =  1.5_prec*mb_v(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(9) =  1.5_prec*mb_w(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0-2,kt0)
!!    			  druvwp(10)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0  ,kt0) &
!!    			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0-1,kt0) &
!!    			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0-2,kt0)
    			  druvwp(6:10) =  1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0  ,kt0,1:5) &
    			                 -2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0-1,kt0,1:5) &
    			                 +0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0-2,kt0,1:5)
    			else
!!    			  druvwp(6) =  0.5_prec*( mb_r(nbt)%a3d(it0,jt0+1,kt0)-mb_r(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(7) =  0.5_prec*( mb_u(nbt)%a3d(it0,jt0+1,kt0)-mb_u(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(8) =  0.5_prec*( mb_v(nbt)%a3d(it0,jt0+1,kt0)-mb_v(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(9) =  0.5_prec*( mb_w(nbt)%a3d(it0,jt0+1,kt0)-mb_w(nbt)%a3d(it0,jt0-1,kt0) )
!!    			  druvwp(10)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0+1,kt0)-mb_p(nbt)%a3d(it0,jt0-1,kt0) )
    			  druvwp(6:10)=  0.5_prec*( mb_bc(nb)%bc(nr)%qpvpack(it0,jt0+1,kt0,1:5) &
    			                           -mb_bc(nb)%bc(nr)%qpvpack(it0,jt0-1,kt0,1:5) )
    			endif
    
    			if(kt0 == 1)then
!!    			  druvwp(11)= -1.5_prec*mb_r(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_r(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_r(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(12)= -1.5_prec*mb_u(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_u(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_u(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(13)= -1.5_prec*mb_v(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_v(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_v(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(14)= -1.5_prec*mb_w(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_w(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_w(nbt)%a3d(it0,jt0,3)
!!    			  druvwp(15)= -1.5_prec*mb_p(nbt)%a3d(it0,jt0,1) &
!!    			              +2.0_prec*mb_p(nbt)%a3d(it0,jt0,2) &
!!    			              -0.5_prec*mb_p(nbt)%a3d(it0,jt0,3)
    			  druvwp(11:15)= -1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,1,1:5) &
    			                 +2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,2,1:5) &
    			                 -0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,3,1:5)

    			elseif(kt0 == nkt)then
!!    			  druvwp(11)=  1.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_r(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_r(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(12)=  1.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_u(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_u(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(13)=  1.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_v(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_v(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(14)=  1.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_w(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_w(nbt)%a3d(it0,jt0,kt0-2)
!!    			  druvwp(15)=  1.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0  ) &
!!    			              -2.0_prec*mb_p(nbt)%a3d(it0,jt0,kt0-1) &
!!    			              +0.5_prec*mb_p(nbt)%a3d(it0,jt0,kt0-2)
    			  druvwp(11:15)=  1.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0  ,1:5) &
    			                 -2.0_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0-1,1:5) &
    			                 +0.5_prec*mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0-2,1:5)
    			else
!!    			  druvwp(11)=  0.5_prec*( mb_r(nbt)%a3d(it0,jt0,kt0+1)-mb_r(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(12)=  0.5_prec*( mb_u(nbt)%a3d(it0,jt0,kt0+1)-mb_u(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(13)=  0.5_prec*( mb_v(nbt)%a3d(it0,jt0,kt0+1)-mb_v(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(14)=  0.5_prec*( mb_w(nbt)%a3d(it0,jt0,kt0+1)-mb_w(nbt)%a3d(it0,jt0,kt0-1) )
!!    			  druvwp(15)=  0.5_prec*( mb_p(nbt)%a3d(it0,jt0,kt0+1)-mb_p(nbt)%a3d(it0,jt0,kt0-1) )
    			  druvwp(11:15)=  0.5_prec*( mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0+1,1:5) &
    			                            -mb_bc(nb)%bc(nr)%qpvpack(it0,jt0,kt0-1,1:5) )
    			endif
    
    			call druvwp_xyz(nxyz,druvwp,druvwpdxyz)
                     !* �Ѽ���ϵ��ԭʼ������һ�׵���ת����ֱ��ϵ��һ�׵���
                     !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                     !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk
                     !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz
    
                !------------------------------------------------------------------------------------!
    			!* ���°ѶԽ�����ԭʼ������ֱ������ϵ�е�һ�׵��������ڱ��������м���ϵ�µ�һ�׵���  !
    
    			volp    = mb_vol(nb)%a3d(i,j,k)
    			nxyz(1) = mb_kcx(nb)%a3d(i,j,k) / volp
    			nxyz(2) = mb_kcy(nb)%a3d(i,j,k) / volp
    			nxyz(3) = mb_kcz(nb)%a3d(i,j,k) / volp
    			nxyz(4) = mb_etx(nb)%a3d(i,j,k) / volp
    			nxyz(5) = mb_ety(nb)%a3d(i,j,k) / volp
    			nxyz(6) = mb_etz(nb)%a3d(i,j,k) / volp
    			nxyz(7) = mb_ctx(nb)%a3d(i,j,k) / volp
    			nxyz(8) = mb_cty(nb)%a3d(i,j,k) / volp
    			nxyz(9) = mb_ctz(nb)%a3d(i,j,k) / volp
    
    !			call get_nijk(nb,nis,njs,nks,i,j,k,nijk) !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
    			call nxyz_nijk(nxyz,nijk,volp) !* nxyz(9): kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz
                                               !* nijk(9): xkc,xet,xct,ykc,yet,yct,zkc,zet,zct
    			call druvwp_xyz(nijk,druvwpdxyz,druvwp)
                     !* druvwpdxyz: drx,dux,dvx,dwx,dpx,dry,duy,dvy,dwy,dpy,drz,duz,dvz,dwz,dpz
    			     !* druvwp(15): dri,dui,dvi,dwi,dpi,drj,duj,dvj,dwj,dpj,drk,duk,dvk,dwk,dpk
    
               !---------------------------------------------------------------------------------!
                !* ���¸����ݶ��� �����߽��֮��Ĳ���                 ------------------------!
    			rm =  s_lr3d(1)*druvwp(1) + s_lr3d(2)*druvwp(6)  + s_lr3d(3)*druvwp(11)
    			um =  s_lr3d(1)*druvwp(2) + s_lr3d(2)*druvwp(7)  + s_lr3d(3)*druvwp(12)
    			vm =  s_lr3d(1)*druvwp(3) + s_lr3d(2)*druvwp(8)  + s_lr3d(3)*druvwp(13)
    			wm =  s_lr3d(1)*druvwp(4) + s_lr3d(2)*druvwp(9)  + s_lr3d(3)*druvwp(14)
    			pm =  s_lr3d(1)*druvwp(5) + s_lr3d(2)*druvwp(10) + s_lr3d(3)*druvwp(15)
    
    			if(dis12 > dis_tao)then
    				dis12 = dis_tao/dis12
    
    				rm = rm * dis12
    				um = um * dis12
    				vm = vm * dis12
    				wm = wm * dis12
    				pm = pm * dis12
    			endif
    
              !* ���¸��� "����" ������ϵ�ֵ                               ------------------------!
    			rm = mb_r(nb)%a3d(ist,jst,kst) + rm*2.0
    			um = mb_u(nb)%a3d(ist,jst,kst) + um*2.0
    			vm = mb_v(nb)%a3d(ist,jst,kst) + vm*2.0
    			wm = mb_w(nb)%a3d(ist,jst,kst) + wm*2.0
    			pm = mb_p(nb)%a3d(ist,jst,kst) + pm*2.0
    
901	format(1x,'��',i3,' ��ĵ�',i2,' �߽���õ����������ֵʧ�ܣ��ܶȣ�',e12.4,3i4)
902	format(1x,'��',i3,' ��ĵ�',i2,' �߽���õ����������ֵʧ�ܣ�ѹ����',e12.4,3i4)
903	format(1x,'��',i3,' ��ĵ�',i2,'�߽繲', i4 ' ����õ����������ֵʧ��')
    
     			if(pm <= 0._prec .or. rm <= 0._prec) then
    			    nerr = nerr + 1
    !				if(nerr<3)write(*,902)nb,nr,pm,i,j,k
    
!!                    rm = mb_r(nbt)%a3d(it,jt,kt)
!!                    um = mb_u(nbt)%a3d(it,jt,kt)
!!                    vm = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. ����ԭʼ���� �� ��Ӧ������������һ����ֵ
!!                    wm = mb_w(nbt)%a3d(it,jt,kt)
!!                    pm = mb_p(nbt)%a3d(it,jt,kt)
                    rm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,1)
                    um = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,2)
                    vm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,3)  !*TGH. ����ԭʼ���� �� ��Ӧ������������һ����ֵ
                    wm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,4)
                    pm = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,5)
    
    			endif
    
                mb_r(nb)%a3d(is,js,ks) = rm
                mb_u(nb)%a3d(is,js,ks) = um
                mb_v(nb)%a3d(is,js,ks) = vm
                mb_w(nb)%a3d(is,js,ks) = wm
                mb_p(nb)%a3d(is,js,ks) = pm
               !* ��㸳ֵ���  ------------------------!
    
!!    			call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC
!!    
!!    			do m=1,cic1
!!    			    call BC_face_dif_turbulence(nb,nb,is,js,ks,ist,jst,kst,i,j,k)
!!    			enddo
               if (cic1 /= 1) then
                  mb_r(nbs)%a3d(i,j,k) = 0.5*(mb_r(nbs)%a3d(is,js,ks) + mb_r(nbs)%a3d(ist,jst,kst))
                  mb_u(nbs)%a3d(i,j,k) = 0.5*(mb_u(nbs)%a3d(is,js,ks) + mb_u(nbs)%a3d(ist,jst,kst))
                  mb_v(nbs)%a3d(i,j,k) = 0.5*(mb_v(nbs)%a3d(is,js,ks) + mb_v(nbs)%a3d(ist,jst,kst))
                  mb_w(nbs)%a3d(i,j,k) = 0.5*(mb_w(nbs)%a3d(is,js,ks) + mb_w(nbs)%a3d(ist,jst,kst))
                  mb_p(nbs)%a3d(i,j,k) = 0.5*(mb_p(nbs)%a3d(is,js,ks) + mb_p(nbs)%a3d(ist,jst,kst))
               end if
               if (nlamtur >= 0 ) then
                  mb_vist(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,6)
                  
                  mb_vist(nbs)%a3d(i,j,k) = 0.5*(mb_vist(nbs)%a3d(is,js,ks) + mb_vist(nbs)%a3d(ist,jst,kst))
                
                  do m=7,nvtot
                     n = m-6
                     mb_qke(nbs)%a4d(is,js,ks,n) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,m)
                     
                     mb_qke(nbs)%a4d(i,j,k,n) = 0.5*(mb_qke(nbs)%a4d(is,js,ks,n) + mb_qke(nbs)%a4d(ist,jst,kst,n))
                  end do
               end if
    
           enddo
           enddo
           enddo
    
    	   if(nerr > 0)write(*,903)nb,nr,nerr
    
!!    	   if(cic1 /= 1) call dif_average(nb,nr,bctype,s_st,s_ed) !* TGH. �߽��ϵ�ֵ���������ƽ����ͬʱ���д�������ģʽ�߽� *!
    
    
        return
    end subroutine boundary_n1_vir2_other    

    subroutine boundary_n1_vir2_parallel(nb,nr,bctype)
       use global_variables
       implicit none
       integer :: nb,nr,bctype
       integer :: nbt,id_src,id_des
       
       nbt = mb_bc(nb)%bc(nr)%nbt
       id_src = mb_pids(nb)
       id_des = mb_pids(nbt)
       if (id_src == id_des) then
          call boundary_n1_vir2(nb,nr,bctype)    
       else
          call boundary_n1_vir2_other(nb,nr,bctype)    
       end if
       
    end subroutine boundary_n1_vir2_parallel

    subroutine boundary_n1_vir3_other(nb,nr,bctype) !�Խӱ߽�����
    !*TGH. ���ʱ���Ȱ�ԭʼ���������ȡ��Ӧ�ص�����ֵ����������mb_qkeҲȡ��Ӧ�ص����ϵ�ֵ��
    !*TGH. Ȼ��call dif_average ��������call BC_face_dif_turbulence�����߽��ϵ�ֵ���ڱ߽������ֵ�ļ�ƽ��
    !*TGH. �Ѿ��޸ĵ��ʺ������߽����������ǻ����ʺϻ�ѧ��Ӧ�����
   
       use global_variables
       implicit none
       integer :: nb,nr,bctype,m,n
       integer :: nbs,s_nd,s_fix,s_lr
       integer :: nbt,t_nd
!!       integer :: t_fix,t_lr
       integer :: is,js,ks,i,j,k,nt
       integer :: it,jt,kt,it0,jt0,kt0,ist,jst,kst
       integer :: s_st(3),s_ed(3),s_lr3d(3),t_lr3d(3)
       integer :: nv_s,nv_t,nfile_par_s,nfile_par_t,num_var_s,num_var_t,nchem_s,nchem_t
       real    :: prim_s1(nl),q_1(nl),zf1,ttt1
   
       do m=1,3
          s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
          s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
   	      s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)
   	      t_lr3d(m) = mb_bc(nb)%bc(nr)%t_lr3d(m)
       end do
       s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
       s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!!       s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
       nbs   = nb !mb_bc(nb)%bc(nr)%nbs         !��� 
   
       nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
       t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k 
!!       t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!!       t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)
   
   
       do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
       do k = s_st(3),s_ed(3)
          it0 = mb_bc(nb)%bc(nr)%image(i,j,k )
          jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,k )
          kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,k )
   
          is = i + s_lr3d(1)     !��һ���ص����λ��
          js = j + s_lr3d(2)
          ks = k + s_lr3d(3)
   
          it = it0 - t_lr3d(1)
          jt = jt0 - t_lr3d(2)
          kt = kt0 - t_lr3d(3)
   
          ist = i - s_lr3d(1)
          jst = j - s_lr3d(2)
          kst = k - s_lr3d(3)

          if (nlamtur >= 0) then
!!   	         mb_vist(nbs)%a3d(is,js,ks) =  mb_vist(nbt)%a3d(it,jt,kt)
   	         mb_vist(nbs)%a3d(is,js,ks) =  mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,6)
   		     mb_vist(nbs)%a3d(i,j,k)    = (mb_vist(nbs)%a3d(is,js,ks) + mb_vist(nbs)%a3d(ist,jst,kst) ) * 0.5
!!   	         do m=1,nvtot
!!   	            mb_qke(nbs)%a4d(is,js,ks,m)  =  mb_qke(nbt)%a4d(it,jt,kt,m)
!!   		        mb_qke(nbs)%a4d(i,j,k,m)     = (mb_qke(nbs)%a4d(is,js,ks,m) + &
!!   		                                        mb_qke(nbs)%a4d(ist,jst,kst,m) ) * 0.5_prec
   	         do m=7,nvtot
   	            n = m-6
   	            mb_qke(nbs)%a4d(is,js,ks,n)  =  mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,m)
   		        mb_qke(nbs)%a4d(i,j,k,n)     = (mb_qke(nbs)%a4d(is,js,ks,n) + mb_qke(nbs)%a4d(ist,jst,kst,n) ) * 0.5
   	         end do
   	      end if
   
   
!!          mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)  !��һ���ص��㸳ֵ
!!          mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
!!          mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
!!          mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
!!          mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)
          mb_r(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,1)  !��һ���ص��㸳ֵ
          mb_u(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,2)
          mb_v(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,3)
          mb_w(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,4)
          mb_p(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,5)
 
   
          is = is + s_lr3d(1)  !��2���ص����λ��
          js = js + s_lr3d(2)
          ks = ks + s_lr3d(3)
   
          it = it - t_lr3d(1)
          jt = jt - t_lr3d(2)
          kt = kt - t_lr3d(3)
   
!!          mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)   !��2���ص��㸳ֵ
!!          mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
!!          mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
!!          mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
!!          mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)
          mb_r(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,1)   !��2���ص��㸳ֵ
          mb_u(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,2)
          mb_v(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,3)
          mb_w(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,4)
          mb_p(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,5)
   
          is = is + s_lr3d(1)  !��3���ص����λ��
          js = js + s_lr3d(2)
          ks = ks + s_lr3d(3)
   
          it = it - t_lr3d(1)
          jt = jt - t_lr3d(2)
          kt = kt - t_lr3d(3)
   
!!          mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)   !��3���ص��㸳ֵ
!!          mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
!!          mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
!!          mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
!!          mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)
          mb_r(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,1)   !��3���ص��㸳ֵ
          mb_u(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,2)
          mb_v(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,3)
          mb_w(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,4)
          mb_p(nbs)%a3d(is,js,ks) = mb_bc(nb)%bc(nr)%qpvpack(it,jt,kt,5)
       end do
       end do
       end do
        
       return
    end subroutine boundary_n1_vir3_other

    subroutine boundary_n1_vir3_parallel(nb,nr,bctype)
       use global_variables
       implicit none
       integer :: nb,nr,bctype
       integer :: nbt,id_src,id_des
       
       nbt = mb_bc(nb)%bc(nr)%nbt
       id_src = mb_pids(nb)
       id_des = mb_pids(nbt)
       if (id_src == id_des) then
          call boundary_n1_vir3(nb,nr,bctype)    
       else
          call boundary_n1_vir3_other(nb,nr,bctype)    
       end if
       
    end subroutine boundary_n1_vir3_parallel
    
    subroutine build_nppos_list
       use global_variables
       use mod_comp_connects
       implicit none
       integer :: i,j,k,m,n,ierr
       integer :: pid,i_st(numprocs)
       
       allocate(nppos_local(numprocs),stat=ierr)
       allocate(ipp_st_local(numprocs),stat=ierr)
       
       nptt_ex = 0
	   do i=1,nptt
	      if (cdate(i)%npp > 2) then
             nptt_ex = nptt_ex + 1
          end if 
       end do
       
       allocate(cdate_ex(nptt_ex),stat=ierr)
       
       nppos = 0
       nppos_local(:) = 0
       k = 0
	   do i=1,nptt
	      if (cdate(i)%npp > 2) then
	         k = k + 1
	         cdate_ex(k)%npp = cdate(i)%npp
             allocate(cdate_ex(k)%ppos_indexs(cdate_ex(k)%npp),stat=ierr) 
             allocate(cdate_ex(k)%nbijk(4,cdate_ex(k)%npp),stat=ierr) 
          
             nppos = nppos + cdate_ex(k)%npp
          
             do j=1,cdate_ex(k)%npp
                cdate_ex(k)%nbijk(:,j) = cdate(i)%nbijk(j,:)
                pid = mb_pids(cdate_ex(k)%nbijk(1,j))
                nppos_local(pid) = nppos_local(pid) + 1
	         end do
	      end if
       end do
       
       allocate(ppos_nbijks(4,nppos),stat=ierr)  
       allocate(ppos_icdate(nppos),stat=ierr)  
       
       ipp_st_local(1) = 0
       do i=2,numprocs
          ipp_st_local(i) = ipp_st_local(i-1) + nppos_local(i-1)
       end do
       
       i_st(:) = ipp_st_local(:)
       do i=1,nptt_ex
          do j=1,cdate_ex(i)%npp
             pid = mb_pids(cdate_ex(i)%nbijk(1,j))
             i_st(pid) = i_st(pid) + 1
             ppos_nbijks(:,i_st(pid)) = cdate_ex(i)%nbijk(:,j)
             ppos_icdate(i_st(pid)) = i
             cdate_ex(i)%ppos_indexs(j) = i_st(pid)
          end do
       end do
       
       n = ipp_st_local(myid+1) + 1
       m = ipp_st_local(myid+1) + nppos_local(myid+1)
       
       allocate(dq_npp(5,nppos),stat=ierr)
       allocate(dq_npp_local(5,n:m),stat=ierr)
       
       allocate(pv_npp(6,nppos),stat=ierr)
       allocate(pv_npp_local(6,n:m),stat=ierr)
       
       call send_message("Build multiple-points connection information successfully!")
       if (myid == master) then
          write(*,*) "$ multiple-points connection:",nptt_ex
       end if
       
    end subroutine build_nppos_list
    
    subroutine build_nppos_list_ex
       use mod_comp_connects
       implicit none
       integer :: i,j,k,m,n,ierr
       integer :: pid,i_st(numprocs)
       integer :: ncount(2)
       type(comp_pnt_t), pointer :: current
       type(simp_pnt_t), pointer :: simpcur
       
       allocate(nppos_local(numprocs),stat=ierr)
       allocate(ipp_st_local(numprocs),stat=ierr)
       
       call connect_pre_parallel
       ncount = comp_list_count(comp_list)
       
       nptt_ex = ncount(2)
       
       allocate(cdate_ex(nptt_ex),stat=ierr)
       
       nppos = 0
       nppos_local(:) = 0
       k = 0
       current => comp_list
       do while ( associated(current) )
          if (current%nsimp_pnts > 2) then
	         k = k + 1
	         cdate_ex(k)%npp = current%nsimp_pnts
             allocate(cdate_ex(k)%ppos_indexs(cdate_ex(k)%npp),stat=ierr) 
             allocate(cdate_ex(k)%nbijk(4,cdate_ex(k)%npp),stat=ierr) 
          
             nppos = nppos + cdate_ex(k)%npp
             
             j = 0
             simpcur => current%simp_head
             do while ( associated(simpcur) )
                j = j + 1
                cdate_ex(k)%nbijk(:,j) = simpcur%nbijk(:)
                pid = mb_pids(cdate_ex(k)%nbijk(1,j))
                nppos_local(pid) = nppos_local(pid) + 1
                simpcur => simpcur%next
             end do
          end if
          current => current%next
       end do
       call comp_list_destroy( comp_list )
       
       allocate(ppos_nbijks(4,nppos),stat=ierr)  
       allocate(ppos_icdate(nppos),stat=ierr)  
       
       ipp_st_local(1) = 0
       do i=2,numprocs
          ipp_st_local(i) = ipp_st_local(i-1) + nppos_local(i-1)
       end do
       
       i_st(:) = ipp_st_local(:)
       do i=1,nptt_ex
          do j=1,cdate_ex(i)%npp
             pid = mb_pids(cdate_ex(i)%nbijk(1,j))
             i_st(pid) = i_st(pid) + 1
             ppos_nbijks(:,i_st(pid)) = cdate_ex(i)%nbijk(:,j)
             ppos_icdate(i_st(pid)) = i
             cdate_ex(i)%ppos_indexs(j) = i_st(pid)
          end do
       end do
       
       n = ipp_st_local(myid+1) + 1
       m = ipp_st_local(myid+1) + nppos_local(myid+1)
       
       allocate(dq_npp(5,nppos),stat=ierr)
       allocate(dq_npp_local(5,n:m),stat=ierr)
       
       allocate(pv_npp(6,nppos),stat=ierr)
       allocate(pv_npp_local(6,n:m),stat=ierr)
       
       call send_message("Build multiple-points connection information successfully!")
       if (myid == master) then
          write(*,*) "$ multiple-points connection:",nptt_ex
       end if
       
    end subroutine build_nppos_list_ex
        
    subroutine communicate_dq_npp
       use global_variables
       implicit none
       integer :: i,j,n,m,ierr
       integer :: nb,it,jt,kt
       integer :: sendcount,recvcounts(numprocs),displs(numprocs)
       real :: vol_p
       
       call exchange_bc_dq_vol
    
       n = ipp_st_local(myid+1) + 1
       m = ipp_st_local(myid+1) + nppos_local(myid+1)
       
       do i=n,m
          nb = ppos_nbijks(1,i)
          it = ppos_nbijks(2,i)
          jt = ppos_nbijks(3,i)
          kt = ppos_nbijks(4,i)
		  vol_p = mb_vol(nb)%a3d(it,jt,kt)
		  do j=1,5
		     dq_npp_local(j,i) = mb_dq(nb)%a4d(j,it,jt,kt)/vol_p
		  end do
       end do
       
       sendcount = (m - n + 1)*5
       recvcounts = nppos_local(:)*5
       displs(:) = ipp_st_local(:)*5
       
       call MPI_ALLGATHERV(dq_npp_local,sendcount,mpi_reprec, &
                           dq_npp,recvcounts,displs,mpi_reprec,MPI_COMM_WORLD,ierr)
                           

       call boundary_match_dq_2pm
       
       call boundary_match_dq_3pp

!!       call send_message("Communicate dq_npp successfully!")
    
    end subroutine communicate_dq_npp
    
    subroutine boundary_match_dq_2pm
       use global_variables
       implicit none
       integer :: nb,nbt,pnb
       integer :: nr,nrmax,ibctype
       integer :: ibeg(3),iend(3),id_src,id_des
       integer :: i,j,k,it,jt,kt,m
       real :: vol_p

       do pnb=1,pnblocks
          nb = pnbindexs(pnb)
          id_src = mb_pids(nb)

          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if ( ibctype < 0 ) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
             
                ibeg = mb_bc(nb)%bc(nr)%s_st
                iend = mb_bc(nb)%bc(nr)%s_ed
                
                if (id_src == id_des) then
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1)
                      it = mb_bc(nb)%bc(nr)%image(i,j,k )
                      jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )
                      kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )
                      vol_p = mb_vol(nb)%a3d(i,j,k)/mb_vol(nbt)%a3d(it,jt,kt)
                      
                      do m=1,5
                         mb_dq(nb)%a4d(m,i,j,k) = 0.5*(mb_dq(nbt)%a4d(m,it,jt,kt)*vol_p + &
                                                       mb_dq(nb)%a4d(m,i,j,k))
                      end do
                   end do
                   end do
                   end do
                else
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1)
                      it = mb_bc(nb)%bc(nr)%image(i,j,k )
                      jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )
                      kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )
                      vol_p = mb_vol(nb)%a3d(i,j,k)
                      
                      do m=1,5
                         mb_dq(nb)%a4d(m,i,j,k) = 0.5*(mb_bc(nb)%bc(nr)%dqvpack_t(it,jt,kt,m)*vol_p + &
                                                       mb_dq(nb)%a4d(m,i,j,k))
                      end do
                   end do
                   end do
                   end do
                end if
             end if
          end do
       end do
    
    end subroutine boundary_match_dq_2pm
    
    subroutine boundary_match_dq_3pp
    !   ��1��nptt �Խ����ϵ㣻��2��cdateΪ�ṹ������npp��Nbijk(:,:)����Ԫ��                !    
    !   ��3��npp��ʾ���жԽ����ϵ�����������NbijkΪλ����Ϣ                              !
    !   ��4��Nbijk(n,1)��Nb��Nbijk(n,2)��i��Nbijk(n,3)��i��Nbijk(n,4)��k��nΪ1��npp        !
    !    (5) �ýṹ�ķ�ʽ���������λ�õĽǶȼ�¼�Խ����ϵ��λ��
    
       use global_variables,only:mb_vol,mb_dq,cdate,nptt,nl
       use define_precision_mod
       implicit none
       integer :: nb,i,j,k,it,jt,kt,npp,m,n,m1,ic,pid
       real(prec):: dqq(1:5),vol_p
  
       n = ipp_st_local(myid+1) + 1
       m1 = ipp_st_local(myid+1) + nppos_local(myid+1)
       
       do ic=n,m1
          i = ppos_icdate(ic)
          npp=cdate_ex(i)%npp
          do m=1,5
             dqq(m)=0.0_prec
          end do
          do j=1,npp
             k = cdate_ex(i)%ppos_indexs(j)
             do m=1,5
                dqq(m)=dqq(m)+dq_npp(m,k)
             end do
          end do
    
          do m=1,5
             dqq(m)=dqq(m)/(npp*1.0_prec)
          end do
          
          nb = ppos_nbijks(1,ic)
          it = ppos_nbijks(2,ic)
          jt = ppos_nbijks(3,ic)
          kt = ppos_nbijks(4,ic)
		  vol_p = mb_vol(nb)%a3d(it,jt,kt)
          do m=1,5
             mb_dq(nb)%a4d(m,it,jt,kt)=dqq(m)*vol_p
          end do
       end do

    end subroutine boundary_match_dq_3pp
    
    subroutine communicate_pv_npp
       use global_variables
       implicit none
       integer :: i,j,n,m,ierr
       integer :: nb,it,jt,kt
       integer :: sendcount,recvcounts(numprocs),displs(numprocs)
       real :: vol_p
       
       call exchange_bc_pv
    
       n = ipp_st_local(myid+1) + 1
       m = ipp_st_local(myid+1) + nppos_local(myid+1)
       
       do i=n,m
          nb = ppos_nbijks(1,i)
          it = ppos_nbijks(2,i)
          jt = ppos_nbijks(3,i)
          kt = ppos_nbijks(4,i)
		  vol_p = 1.0  !!mb_vol(nb)%a3d(it,jt,kt)
		  pv_npp_local(1,i) = mb_r(nb)%a3d(it,jt,kt)/vol_p
		  pv_npp_local(2,i) = mb_u(nb)%a3d(it,jt,kt)/vol_p
		  pv_npp_local(3,i) = mb_v(nb)%a3d(it,jt,kt)/vol_p
		  pv_npp_local(4,i) = mb_w(nb)%a3d(it,jt,kt)/vol_p
		  pv_npp_local(5,i) = mb_p(nb)%a3d(it,jt,kt)/vol_p
		  pv_npp_local(6,i) = mb_t(nb)%a3d(it,jt,kt)/vol_p
		  
!!          write(30+myid,'(4(1x,i4),6(1x,e12.5))')ppos_nbijks(:,i),pv_npp_local(:,i)
       end do
       
       sendcount = (m - n + 1)*6
       recvcounts = nppos_local(:)*6
       displs(:) = ipp_st_local(:)*6
       
       call MPI_ALLGATHERV(pv_npp_local,sendcount,mpi_reprec, &
                           pv_npp,recvcounts,displs,mpi_reprec,MPI_COMM_WORLD,ierr)
                           
!!       if (myid == master) then
!!          do i=1,nppos
!!             write(41,'(4(1x,i4),6(1x,e12.5))')ppos_nbijks(:,i),pv_npp(:,i)
!!          end do
!!          stop
!!       end if
                           
       call boundary_match_pv_2pm
       
       call boundary_match_pv_3pp
                           
!!       call send_message("Communicate pv_npp successfully!")
    
    end subroutine communicate_pv_npp
    
    subroutine boundary_match_pv_2pm
       use global_variables
       implicit none
       integer :: nb,nbt,pnb
       integer :: nr,nrmax,ibctype
       integer :: ibeg(3),iend(3),id_src,id_des
       integer :: i,j,k,it,jt,kt,m

       do pnb=1,pnblocks
          nb = pnbindexs(pnb)
          id_src = mb_pids(nb)

          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if ( ibctype < 0 ) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
             
                ibeg = mb_bc(nb)%bc(nr)%s_st
                iend = mb_bc(nb)%bc(nr)%s_ed
                
                if (id_src == id_des) then
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1)
                      it = mb_bc(nb)%bc(nr)%image(i,j,k )
                      jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )
                      kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )
                      
                      mb_r(nb)%a3d(i,j,k) = 0.5*(mb_r(nbt)%a3d(it,jt,kt) + mb_r(nb)%a3d(i,j,k))
                      mb_u(nb)%a3d(i,j,k) = 0.5*(mb_u(nbt)%a3d(it,jt,kt) + mb_u(nb)%a3d(i,j,k))
                      mb_v(nb)%a3d(i,j,k) = 0.5*(mb_v(nbt)%a3d(it,jt,kt) + mb_v(nb)%a3d(i,j,k))
                      mb_w(nb)%a3d(i,j,k) = 0.5*(mb_w(nbt)%a3d(it,jt,kt) + mb_w(nb)%a3d(i,j,k))
                      mb_p(nb)%a3d(i,j,k) = 0.5*(mb_p(nbt)%a3d(it,jt,kt) + mb_p(nb)%a3d(i,j,k))
                      mb_t(nb)%a3d(i,j,k) = 0.5*(mb_t(nbt)%a3d(it,jt,kt) + mb_t(nb)%a3d(i,j,k))
                   end do
                   end do
                   end do
                else
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1)
                      it = mb_bc(nb)%bc(nr)%image(i,j,k )
                      jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )
                      kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )
                   
                      mb_r(nb)%a3d(i,j,k) = 0.5*(mb_bc(nb)%bc(nr)%qpvpack_t(it,jt,kt,1) + mb_r(nb)%a3d(i,j,k))
                      mb_u(nb)%a3d(i,j,k) = 0.5*(mb_bc(nb)%bc(nr)%qpvpack_t(it,jt,kt,2) + mb_u(nb)%a3d(i,j,k))
                      mb_v(nb)%a3d(i,j,k) = 0.5*(mb_bc(nb)%bc(nr)%qpvpack_t(it,jt,kt,3) + mb_v(nb)%a3d(i,j,k))
                      mb_w(nb)%a3d(i,j,k) = 0.5*(mb_bc(nb)%bc(nr)%qpvpack_t(it,jt,kt,4) + mb_w(nb)%a3d(i,j,k))
                      mb_p(nb)%a3d(i,j,k) = 0.5*(mb_bc(nb)%bc(nr)%qpvpack_t(it,jt,kt,5) + mb_p(nb)%a3d(i,j,k))
                      mb_t(nb)%a3d(i,j,k) = 0.5*(mb_bc(nb)%bc(nr)%qpvpack_t(it,jt,kt,6) + mb_t(nb)%a3d(i,j,k))
                   end do
                   end do
                   end do
                end if
             end if
          end do
       end do
           
    end subroutine boundary_match_pv_2pm
    
    subroutine boundary_match_pv_3pp
    !   ��1��nptt �Խ����ϵ㣻��2��cdateΪ�ṹ������npp��Nbijk(:,:)����Ԫ��                !    
    !   ��3��npp��ʾ���жԽ����ϵ�����������NbijkΪλ����Ϣ                              !
    !   ��4��Nbijk(n,1)��Nb��Nbijk(n,2)��i��Nbijk(n,3)��i��Nbijk(n,4)��k��nΪ1��npp        !
    !    (5) �ýṹ�ķ�ʽ���������λ�õĽǶȼ�¼�Խ����ϵ��λ��
    
       use global_variables,only:mb_vol,mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,cdate,nptt,nl
       use define_precision_mod
       implicit none
       integer :: nb,i,j,k,it,jt,kt,npp,m,n,m1,ic,pid
       real(prec):: dqq(1:6),vol_p
  
       n = ipp_st_local(myid+1) + 1
       m1 = ipp_st_local(myid+1) + nppos_local(myid+1)
       
       do ic=n,m1
          i = ppos_icdate(ic)
          npp=cdate_ex(i)%npp
          do m=1,6
             dqq(m)=0.0_prec
          end do
          do j=1,npp
             k = cdate_ex(i)%ppos_indexs(j)
             do m=1,6
                dqq(m)=dqq(m)+pv_npp(m,k)
             end do
          end do
    
          do m=1,6
             dqq(m)=dqq(m)/(npp*1.0_prec)
          end do
          
          nb = ppos_nbijks(1,ic)
          it = ppos_nbijks(2,ic)
          jt = ppos_nbijks(3,ic)
          kt = ppos_nbijks(4,ic)
		  vol_p = 1.0   !!mb_vol(nb)%a3d(it,jt,kt)
          
          mb_r(nb)%a3d(it,jt,kt)=dqq(1)*vol_p
          mb_u(nb)%a3d(it,jt,kt)=dqq(2)*vol_p
          mb_v(nb)%a3d(it,jt,kt)=dqq(3)*vol_p
          mb_w(nb)%a3d(it,jt,kt)=dqq(4)*vol_p
          mb_p(nb)%a3d(it,jt,kt)=dqq(5)*vol_p
          mb_t(nb)%a3d(it,jt,kt)=dqq(6)*vol_p
       end do

    end subroutine boundary_match_pv_3pp    
    
    
    subroutine exchange_bc_pv
       use global_variables                                                                        
       implicit none
       integer :: nb,pnb,i,j,k,m,ierr
       integer :: nr,nregions,ibctype
       integer :: iseq,nbt,nrt,packsize
       integer :: id_src,id_des,tag_seq
       integer :: ibeg(3),iend(3),idir,inrout
       integer :: status(MPI_STATUS_SIZE)
       
       iseq = 0
       !!do pnb=1,pnblocks
       !!   nb = pnbindexs(pnb)
       do nb=1,nblocks
          id_src = mb_pids(nb)
          nregions = mb_bc(nb)%nregions
          do nr=1,nregions
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if (ibctype < 0) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
                if (id_src /= id_des) then
                   iseq = iseq + 1
                   tag_seq = iseq
                   nrt = mb_bc(nb)%bc(nr)%ibcwin
                   
                   if (myid == id_src-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
                      
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1)
                         mb_bc(nb)%bc(nr)%qpv_t(i,j,k,1) = mb_r(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv_t(i,j,k,2) = mb_u(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv_t(i,j,k,3) = mb_v(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv_t(i,j,k,4) = mb_w(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv_t(i,j,k,5) = mb_p(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%qpv_t(i,j,k,6) = mb_t(nb)%a3d(i,j,k)
                      end do
                      end do
                      end do
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*6
                      
                      call MPI_SEND(mb_bc(nb)%bc(nr)%qpv_t,packsize,mpi_reprec, &
                                    id_des-1,tag_seq,MPI_COMM_WORLD,ierr)
                   end if
                   
                   if (myid == id_des-1) then
                      ibeg = mb_bc(nbt)%bc(nrt)%t_st
                      iend = mb_bc(nbt)%bc(nrt)%t_ed
                      idir = mb_bc(nbt)%bc(nrt)%t_nd
                      inrout = mb_bc(nbt)%bc(nrt)%t_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*6
                      
                      call MPI_RECV(mb_bc(nbt)%bc(nrt)%qpvpack_t,packsize,mpi_reprec, &
                                    id_src-1,tag_seq,MPI_COMM_WORLD,status,ierr)
                   end if
                   
!!                   call synchronize(81)
                   
                end if
             end if
          end do
       end do

!!       call send_message("Exchange the BC pv successfully!")
       
    end subroutine exchange_bc_pv
        
    subroutine exchange_bc_dq_vol
       use global_variables                                                                        
       implicit none
       integer :: nb,pnb,i,j,k,m,ierr
       integer :: nr,nregions,ibctype
       integer :: iseq,nbt,nrt,packsize
       integer :: id_src,id_des,tag_seq
       integer :: ibeg(3),iend(3),idir,inrout
       integer :: status(MPI_STATUS_SIZE)
       real :: vol_p
       
       iseq = 0
       !!do pnb=1,pnblocks
       !!   nb = pnbindexs(pnb)
       do nb=1,nblocks
          id_src = mb_pids(nb)
          nregions = mb_bc(nb)%nregions
          do nr=1,nregions
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if (ibctype < 0) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
                if (id_src /= id_des) then
                   iseq = iseq + 1
                   tag_seq = iseq
                   nrt = mb_bc(nb)%bc(nr)%ibcwin
                   
                   if (myid == id_src-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
                      
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1) 
                         vol_p = mb_vol(nb)%a3d(i,j,k)
                         mb_bc(nb)%bc(nr)%dqv_t(i,j,k,1) = mb_dq(nb)%a4d(1,i,j,k)/vol_p
                         mb_bc(nb)%bc(nr)%dqv_t(i,j,k,2) = mb_dq(nb)%a4d(2,i,j,k)/vol_p
                         mb_bc(nb)%bc(nr)%dqv_t(i,j,k,3) = mb_dq(nb)%a4d(3,i,j,k)/vol_p
                         mb_bc(nb)%bc(nr)%dqv_t(i,j,k,4) = mb_dq(nb)%a4d(4,i,j,k)/vol_p
                         mb_bc(nb)%bc(nr)%dqv_t(i,j,k,5) = mb_dq(nb)%a4d(5,i,j,k)/vol_p
                      end do
                      end do
                      end do
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*5
                      
                      call MPI_SEND(mb_bc(nb)%bc(nr)%dqv_t,packsize,mpi_reprec, &
                                    id_des-1,tag_seq,MPI_COMM_WORLD,ierr)
                   end if
                   
                   if (myid == id_des-1) then
                      ibeg = mb_bc(nbt)%bc(nrt)%t_st
                      iend = mb_bc(nbt)%bc(nrt)%t_ed
                      idir = mb_bc(nbt)%bc(nrt)%t_nd
                      inrout = mb_bc(nbt)%bc(nrt)%t_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*5
                      
                      call MPI_RECV(mb_bc(nbt)%bc(nrt)%dqvpack_t,packsize,mpi_reprec, &
                                    id_src-1,tag_seq,MPI_COMM_WORLD,status,ierr)
                   end if
                   
!!                   call synchronize(91)
                   
                end if
             end if
          end do
       end do
       
!!       call send_message("Exchange the BC dq_vol successfully!")
       
    end subroutine exchange_bc_dq_vol   
    
    subroutine exchange_bc_dq
       use global_variables                                                                        
       implicit none
       integer :: nb,pnb,i,j,k,m,ierr
       integer :: nr,nregions,ibctype
       integer :: iseq,nbt,nrt,packsize
       integer :: id_src,id_des,tag_seq
       integer :: ibeg(3),iend(3),idir,inrout
       integer :: status(MPI_STATUS_SIZE)
       
       iseq = 0
       !!do pnb=1,pnblocks
       !!   nb = pnbindexs(pnb)
       do nb=1,nblocks
          id_src = mb_pids(nb)
          nregions = mb_bc(nb)%nregions
          do nr=1,nregions
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if (ibctype < 0) then
                nbt = mb_bc(nb)%bc(nr)%nbt
                id_des = mb_pids(nbt)
                if (id_src /= id_des) then
                   iseq = iseq + 1
                   tag_seq = iseq
                   nrt = mb_bc(nb)%bc(nr)%ibcwin
                   
                   if (myid == id_src-1) then
                      ibeg = mb_bc(nb)%bc(nr)%s_st
                      iend = mb_bc(nb)%bc(nr)%s_ed
                      idir = mb_bc(nb)%bc(nr)%s_nd
                      inrout = mb_bc(nb)%bc(nr)%s_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
                      
                      do k=ibeg(3),iend(3)
                      do j=ibeg(2),iend(2)
                      do i=ibeg(1),iend(1) 
                         mb_bc(nb)%bc(nr)%dqv(i,j,k,1) = mb_dq(nb)%a4d(1,i,j,k)
                         mb_bc(nb)%bc(nr)%dqv(i,j,k,2) = mb_dq(nb)%a4d(2,i,j,k)
                         mb_bc(nb)%bc(nr)%dqv(i,j,k,3) = mb_dq(nb)%a4d(3,i,j,k)
                         mb_bc(nb)%bc(nr)%dqv(i,j,k,4) = mb_dq(nb)%a4d(4,i,j,k)
                         mb_bc(nb)%bc(nr)%dqv(i,j,k,5) = mb_dq(nb)%a4d(5,i,j,k)
                      end do
                      end do
                      end do
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*5
                      
                      call MPI_SEND(mb_bc(nb)%bc(nr)%dqv,packsize,mpi_reprec, &
                                    id_des-1,tag_seq,MPI_COMM_WORLD,ierr)
                   end if
                   
                   if (myid == id_des-1) then
                      ibeg = mb_bc(nbt)%bc(nrt)%t_st
                      iend = mb_bc(nbt)%bc(nrt)%t_ed
                      idir = mb_bc(nbt)%bc(nrt)%t_nd
                      inrout = mb_bc(nbt)%bc(nrt)%t_lr
!!                      ibeg(idir) = ibeg(idir) - inrout
!!                      iend(idir) = iend(idir) - inrout
                      
                      packsize = product(iend(:)-ibeg(:)+1,1)*5
                      
                      call MPI_RECV(mb_bc(nbt)%bc(nrt)%dqvpack,packsize,mpi_reprec, &
                                    id_src-1,tag_seq,MPI_COMM_WORLD,status,ierr)
                   end if
                   
!!                   call synchronize(91)
                   
                end if
             end if
          end do
       end do
       
!!       call send_message("Exchange the BC dq successfully!")
       
    end subroutine exchange_bc_dq
            
    subroutine bc_n1_cic_pre_other(nb,nr,bctype) 
    !*TGH. �Խӱ߽����� ǰ���������������Խӱ߽� ֻ�������޲��
    !*TGH. ��Ŀ����ڶԽ����ϵ�dq���ֵ�Դ��������
       use global_variables,only : mb_bc,mb_dq,nblocks,method,nl !!,mb_r,mb_u,mb_v,mb_w,mb_p
       implicit none
       integer :: nb,nr,bctype,m,n,n1
       integer :: nbs,s_nd,s_fix,s_lr
       integer :: nbt,t_nd,t_fix,t_lr
       integer :: is,js,ks,i,j,k,nt,ntarg
       integer :: it0,jt0,kt0 !,it,jt,kt
       integer :: s_st(3),s_ed(3)
       
       do m=1,3
          s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
          s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
       enddo
       s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
       s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
       s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
       nbs   = mb_bc(nb)%bc(nr)%nbs             !��� 
       
       nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
       t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k 
       t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
       t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)

       !*TGH. ���� ע��:��
       !*TGH. Դ�������ԭʼ���� ȡ Ŀ����ϵĶ�Ӧֵ������������һ����ֵ����border_all���Ѿ���ֵ��
       !*TGH. Դ�������dq ȡ Ŀ����϶Խ����ϵ�ֵ ��Ϊ�˱��ں��������߽�ʹ�ã�
       do i = s_st(1),s_ed(1)
          do j = s_st(2),s_ed(2)
             do k = s_st(3),s_ed(3)
                do n=1,1
                   is = i + mb_bc(nb)%bc(nr)%s_lr3d(1)*n
                   js = j + mb_bc(nb)%bc(nr)%s_lr3d(2)*n  !*TGH. Դ��������һ�����
                   ks = k + mb_bc(nb)%bc(nr)%s_lr3d(3)*n
                   it0= mb_bc(nb)%bc(nr)%image(i,j,k )
                   jt0= mb_bc(nb)%bc(nr)%jmage(i,j,k )   !*TGH. Ŀ���棬���ʱ�ڶԽ����ϣ����ʱ��������һ��
                   kt0= mb_bc(nb)%bc(nr)%kmage(i,j,k )
!                   nt = n - 1 + method
!                   it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)*nt
!                   jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)*nt !*TGH. Ŀ���棬���ʱ��������һ�㣬���ʱ���䣨����������
!                   kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)*nt
!                   mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)
!                   mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)  !*TGH. Դ��������һ��
!                   mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. Ŀ�����������һ��
!                   mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)  !*TGH. ��Դ������ȡĿ�������Ӧ���ֵ
!                   mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)
       
                   do n1=1,5
!!	   				mb_dq(nbs)%a4d(n1,is,js,ks) = mb_dq(nbt)%a4d(n1,it0,jt0,kt0)
	   				mb_dq(nbs)%a4d(n1,is,js,ks) = mb_bc(nb)%bc(nr)%dqvpack(it0,jt0,kt0,n1)
	   			!*TGH. ע�⣬��Դ������ֵ����¼Ŀ������ڶԽ����ϵ�RHS
	   			!*TGH. ���ǵ���Ҫʹ��ʱ�����Բ�����Ŀ������Ϣ��ֱ�Ӳ�������ϵ�ֵ����
	   			enddo
	   													
                enddo
             enddo
          enddo
       enddo
    end subroutine bc_n1_cic_pre_other    
    
    subroutine bc_n1_cic_pre_parallel(nb,nr,bctype)
       use global_variables
       implicit none
       integer :: nb,nr,bctype
       integer :: nbt,id_src,id_des
       
       nbt = mb_bc(nb)%bc(nr)%nbt
       id_src = mb_pids(nb)
       id_des = mb_pids(nbt)
       if (id_src == id_des) then
          call bc_n1_cic_pre(nb,nr,bctype)    
       else
          call bc_n1_cic_pre_other(nb,nr,bctype)    
       end if
       
    end subroutine bc_n1_cic_pre_parallel    
    
    
    subroutine dif_cic_new_other(nb,nr) !�ж�����ֵʱ���ǵ���������ʱΪ��Сֵ�����
    !*TGH. �����Խ��޸�rhs����dq��������ϵ�ֵ��¼�Խ����϶�Ӧ��Ŀ����rhs
    !*TGH. ��Ҫ�õ��Խ����Լ�����ϵ�rhs1�����������õ��Խ����ϵ�dq��
    !*tgh. �����Խ���Ȼ�ü�ƽ��
    !*TGH. ʹ��ǰ����Recast������������
    !* ����ϵ������Ϊeps���ڡ�subroutine Aps_dq_new�� ��
    !* ��ƣ�Ϳ����                             ----------------------------------!
    !* ���ԣ�Ϳ����  2009.3                     ----------------------------------!

    use define_precision_mod
    use global_const,only:rmin_limit,pmin_limit,rmax_limit,pmax_limit
    use global_variables ,only : mb_bc,mb_r,mb_u,mb_v,mb_w,mb_p,mb_q,mb_dq, &
                                  mb_kcx,mb_kcy,mb_kcz,mb_kct,&
                                  mb_etx,mb_ety,mb_etz,mb_ett,&
                                  mb_ctx,mb_cty,mb_ctz,mb_ctt,&
                                  kcx,kcy,kcz,kct,etx,ety,etz,ett,ctx,cty,ctz,ctt,&
                                  mb_vol,mb_dim,nl,gama,r,u,v,w,p,dq,q,mb_dtdt,sml_sss


    implicit none
    integer :: s_st(3),s_ed(3),s_lr3d(3),nb,nr,nsurf !,bctype
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: i,j,k,is,js,ks,it,jt,kt,m,n1,n2
    real(prec) :: nx,ny,nz,nt,nvols,nvolt
    real(prec) :: ARRAY_P1(5,5),ARRAY_P(5,5),ARRAY_Asp(5,5),ARRAY_Asn(5,5),ary_c(5,5)
    real(prec) :: rm,um,vm,wm,pm,mm,cm,tm,cgm,ct,qbar1(5)
    real(prec) :: l1,l4,l5
    real(prec) :: qsp(5),qsn(5),qspn(5)
    integer :: is1,js1,ks1,is2,js2,ks2
    integer :: s01,s02
    real(prec) :: dqq(nl)
    real(prec) :: nxyz(4)

    do m=1,3
       s_st(m)   = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m)   = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)
    enddo

    nsurf  = mb_bc(nb)%bc(nr)%s_nd         

    nbs   = mb_bc(nb)%bc(nr)%nbs             !��� 
    s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
    s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
    s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)

    nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
    t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k 
    t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
    t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)

             it = mb_bc(nb)%bc(nr)%image(i,j,k )
             jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )    !*TGH. �Խ��棨Ŀ��飩�ϵ�����
             kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )

             is1 = i + mb_bc(nb)%bc(nr)%s_lr3d(1)
             js1 = j + mb_bc(nb)%bc(nr)%s_lr3d(2)  !*TGH. �Խ��棨Դ��������һ�������
             ks1 = k + mb_bc(nb)%bc(nr)%s_lr3d(3)

             is2 = i - mb_bc(nb)%bc(nr)%s_lr3d(1)
             js2 = j - mb_bc(nb)%bc(nr)%s_lr3d(2)  !*TGH. �Խ��棨Դ����������һ�������
             ks2 = k - mb_bc(nb)%bc(nr)%s_lr3d(3)

             if ( nsurf == 1 ) then    !*tgh. ���ݶԽ����������ϵķ���ʸ��
                nx = kcx(i,j,k) ! mb_kcx(nb)%a3d(i,j,k)
                ny = kcy(i,j,k) ! mb_kcy(nb)%a3d(i,j,k)
                nz = kcz(i,j,k) ! mb_kcz(nb)%a3d(i,j,k)  !*tgh. Դ���ϵķ���ʸ��
                nt = kct(i,j,k) ! mb_kct(nb)%a3d(i,j,k)
             elseif ( nsurf == 2 ) then
                nx = etx(i,j,k) ! mb_etx(nb)%a3d(i,j,k)
                ny = ety(i,j,k) ! mb_ety(nb)%a3d(i,j,k)
                nz = etz(i,j,k) ! mb_etz(nb)%a3d(i,j,k)
                nt = ett(i,j,k) ! mb_ett(nb)%a3d(i,j,k)
             else
                nx = ctx(i,j,k) ! mb_ctx(nb)%a3d(i,j,k) 
                ny = cty(i,j,k) ! mb_cty(nb)%a3d(i,j,k) 
                nz = ctz(i,j,k) ! mb_ctz(nb)%a3d(i,j,k) 
                nt = ctt(i,j,k) ! mb_ctt(nb)%a3d(i,j,k) 
             endif

             nvols = mb_vol(nbs)%a3d(i,j,k)
!!             nvolt = mb_vol(nbt)%a3d(it,jt,kt)
             nvolt = mb_bc(nb)%bc(nr)%volpack(it,jt,kt)


             rm = r(i,j,k)  !mb_r(nb)%a3d(i,j,k)  !*tgh. Դ���ϵ�ԭʼ����
             um = u(i,j,k)  !mb_u(nb)%a3d(i,j,k)
             vm = v(i,j,k)  !mb_v(nb)%a3d(i,j,k)
             wm = w(i,j,k)  !mb_w(nb)%a3d(i,j,k)
             pm = p(i,j,k)  !mb_p(nb)%a3d(i,j,k)

             qbar1(1) = rm
             qbar1(2) = um
             qbar1(3) = vm
             qbar1(4) = wm
             qbar1(5) = pm
!_______

             !������
             cm  = sqrt(gama*pm/rm)
             ct  = nx*um + ny*vm + nz*wm + nt  !*tgh. ����ٶ�
             cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)   !Ҫ���ǵ�cgmΪ������

             l1 = ct  !+ 1.e-30
             l4 = ct + cm * cgm !+ 1.e-30
             l5 = ct - cm * cgm !+ 1.e-30

             l1 = l1 * s_lr  
             l4 = l4 * s_lr   !*tgh. ���ݶԽӷ����ж�����ֵ����Դ�飩
             l5 = l5 * s_lr 

!*********************************
             !_______
             if(abs(l4)>abs(l5))then    
             s01 = 0.0
             s02 = 1.0  !*TGH. ����ٶȴ�����������������ٶ�������ͬ����
             else
             s01 = 1.0  !*TGH. ����ٶ�С��������
             s02 = 0.0
             endif
             !________

             ARRAY_Asp = 0._prec !1.0d-30        !__A~ !*TGH. ע�⣬ȡ��ֵ���������
             ARRAY_Asn = 0._prec !1.0d-30        !__A+
             
!______���濪ʼ��A~*RHS  A+*RHS

             qsp = 0.0
             qsn = 0.0
         
             nxyz(1) = nx
             nxyz(2) = ny
             nxyz(3) = nz
             nxyz(4) = nt

             do n1 = 1,5
!               dqq(n1) = mb_dq(nb)%a4d(n1,i,j,k)    !*TGH. Դ���ϵĶԽ���
                dqq(n1) = dq(n1,i,j,k)    !*TGH. Դ���ϵĶԽ���
             enddo
             call Aps_dq_new(qbar1,nxyz,s_lr,dqq,qsp,1)  !*TGH. ��A~��A+������ʾ��Ȼ��ֱ�Ӷ�RHS����ת��
             !TGH. ����Ϊ ԭʼ���������������Խӷ���rhs�����������ϵ�rhs��

                do n1 = 1,5
                   dqq(n1) = dq(n1,is1,js1,ks1)*nvols/nvolt  !*TGH. Դ��������һ���RHS
                enddo
              
                call Aps_dq_new(qbar1,nxyz,s_lr,dqq,qsn,-1)
!____         
              
!_____        
                do n1 = 1,5
!                   mb_dq(nb)%a4d(n1,i,j,k) = (qsp(n1)+qsn(n1))    !*TGH. ??? ʹ���������RHS?��dq?
                    dq(n1,i,j,k) = (qsp(n1)+qsn(n1))    !*TGH. ??? ʹ���������RHS?��dq?
                enddo
              
                !*TGH. ����ϵ�dq��ʹ�ú���������
                do n1 = 1,5
                   dq(n1,is1,js1,ks1) = 0._prec  !*TGH. Դ��������һ���dq
                enddo
                !*TGH. end ��ʹ�ú�����ϵ�dq��������


           !*TGH. ���µ��������ж������Խ��Ƿ����
           if(.false.)then
                mb_q(nb)%a4d(1,i,j,k) = rm
                mb_q(nb)%a4d(2,i,j,k) = rm*um
                mb_q(nb)%a4d(3,i,j,k) = rm*vm
                mb_q(nb)%a4d(4,i,j,k) = rm*wm
                mb_q(nb)%a4d(5,i,j,k) = pm/(gama - 1.d0) + 0.5d0*rm*( um*um + vm*vm + wm*wm )
         
                do n1 = 1,5
                    qspn(n1) = mb_q(nb)%a4d(n1,i,j,k) - dq(n1,i,j,k)*mb_dtdt(nb)%a3d(i,j,k) 
                enddo
         
         
                rm = qspn(1)          
                um = qspn(2)/qspn(1)
                vm = qspn(3)/qspn(1)
                wm = qspn(4)/qspn(1)
                pm = (qspn(5) - 0.5*rm*(um*um + vm*vm + wm*wm))*(gama -1.0 )
         
         
                if(rm<rmin_limit .or. pm<pmin_limit .or. rm>rmax_limit .or. pm>pmax_limit )then
                   write(*,*)'error?     rm = ',rm,'pm= ',pm
                   write(*,*)'nb= ',nb,i,j,k
                   rm = mb_r(nb)%a3d(i,j,k) != 0.5*( r(is,js,ks) + r(it,jt,kt) )
                   um = mb_u(nb)%a3d(i,j,k) != 0.5*( u(is,js,ks) + u(it,jt,kt) )
                   vm = mb_v(nb)%a3d(i,j,k) != 0.5*( v(is,js,ks) + v(it,jt,kt) )
                   wm = mb_w(nb)%a3d(i,j,k) != 0.5*( w(is,js,ks) + w(it,jt,kt) )
                   pm = mb_p(nb)%a3d(i,j,k) != 0.5*( p(is,js,ks) + p(it,jt,kt) )
                endif
           endif
           !*TGH. end �ж������Խ��Ƿ����
         
           call BC_face_dif_turbulence(nb,nb,is1,js1,ks1,is2,js2,ks2,i,j,k) !*tgh. �����Խ���Ȼ�ü�ƽ��
          enddo
       enddo
    enddo
    return
    end subroutine dif_cic_new_other    
    
    subroutine dif_cic_new_parallel(nb,nr)
       use global_variables
       implicit none
       integer :: nb,nr
       integer :: nbt,id_src,id_des
       
       nbt = mb_bc(nb)%bc(nr)%nbt
       id_src = mb_pids(nb)
       id_des = mb_pids(nbt)
       if (id_src == id_des) then
          call dif_cic_new(nb,nr)    
       else
          call dif_cic_new_other(nb,nr)    
       end if
       
    end subroutine dif_cic_new_parallel       
    
    subroutine read_flow_file_parallel
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
       real :: tmp,prim(nl),q_q(nl)
       real,pointer :: qpv_tmp(:,:,:,:)
	
	   if (myid == master) then
       open(1,file=flowname,status='unknown',form='unformatted')
       read(1)nout,wholetime,tmp,tmp,tmp,tmp,tmp,tmp
       end if
       
       call MPI_BCAST(nout     ,1,mpi_inprec,master,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(wholetime,1,mpi_reprec,master,MPI_COMM_WORLD,ierr)
          
       do nb=1,nblocks
          pid = mb_pids(nb)
       
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          if (master == pid-1) then
	         if (myid == master) then
                read(1)(((mb_r(nb)%a3d(i,j,k),mb_u(nb)%a3d(i,j,k),mb_v(nb)%a3d(i,j,k), &
                          mb_w(nb)%a3d(i,j,k),mb_p(nb)%a3d(i,j,k),mb_t(nb)%a3d(i,j,k), &
                          i=1,idim),j=1,jdim),k=1,kdim)
                          
                if (ndualtst > 0) then
                   do k=1,kdim
                   do j=1,jdim
                   do i=1,idim
                      prim(1) = mb_r(nb)%a3d(i,j,k)
                      prim(2) = mb_u(nb)%a3d(i,j,k)
                      prim(3) = mb_v(nb)%a3d(i,j,k)
                      prim(4) = mb_w(nb)%a3d(i,j,k)
                      prim(5) = mb_p(nb)%a3d(i,j,k)
                      call prim_to_q(prim,q_q,gama)
                      mb_q(nb)%a4d(1,i,j,k) = q_q(1)
                      mb_q(nb)%a4d(2,i,j,k) = q_q(2)
                      mb_q(nb)%a4d(3,i,j,k) = q_q(3)
                      mb_q(nb)%a4d(4,i,j,k) = q_q(4)
                      mb_q(nb)%a4d(5,i,j,k) = q_q(5)
                      
                      do m=1,nl
                         mb_qnc(nb)%a4d(m,i,j,k) = mb_q(nb)%a4d(m,i,j,k)
                      end do
                   
                   enddo
                   enddo
                   enddo
                end if         
             
             end if
             
          else
	         if (myid == master) then
                allocate(qpv_tmp(6,idim,jdim,kdim),stat=ierr)
                read(1)((((qpv_tmp(m,i,j,k),m=1,6),i=1,idim),j=1,jdim),k=1,kdim)
                
                packsize = idim*jdim*kdim*6
                      
                call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,ierr)
                
                deallocate(qpv_tmp,stat=ierr)
             endif
             
	         if (myid == pid-1) then
                allocate(qpv_tmp(6,idim,jdim,kdim),stat=ierr)
                
                packsize = idim*jdim*kdim*6
                
                call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
                              
                do k=1,kdim
                do j=1,jdim
                do i=1,idim
                   mb_r(nb)%a3d(i,j,k) = qpv_tmp(1,i,j,k)
                   mb_u(nb)%a3d(i,j,k) = qpv_tmp(2,i,j,k)
                   mb_v(nb)%a3d(i,j,k) = qpv_tmp(3,i,j,k)
                   mb_w(nb)%a3d(i,j,k) = qpv_tmp(4,i,j,k)
                   mb_p(nb)%a3d(i,j,k) = qpv_tmp(5,i,j,k)
                   mb_t(nb)%a3d(i,j,k) = qpv_tmp(6,i,j,k)
                end do
                end do
                end do
                
                deallocate(qpv_tmp,stat=ierr)
                
                if (ndualtst > 0) then
                   do k=1,kdim
                   do j=1,jdim
                   do i=1,idim
                      prim(1) = mb_r(nb)%a3d(i,j,k)
                      prim(2) = mb_u(nb)%a3d(i,j,k)
                      prim(3) = mb_v(nb)%a3d(i,j,k)
                      prim(4) = mb_w(nb)%a3d(i,j,k)
                      prim(5) = mb_p(nb)%a3d(i,j,k)
                      call prim_to_q(prim,q_q,gama)
                      mb_q(nb)%a4d(1,i,j,k) = q_q(1)
                      mb_q(nb)%a4d(2,i,j,k) = q_q(2)
                      mb_q(nb)%a4d(3,i,j,k) = q_q(3)
                      mb_q(nb)%a4d(4,i,j,k) = q_q(4)
                      mb_q(nb)%a4d(5,i,j,k) = q_q(5)
                      
                      do m=1,nl
                         mb_qnc(nb)%a4d(m,i,j,k) = mb_q(nb)%a4d(m,i,j,k)
                      end do
                   
                   enddo
                   enddo
                   enddo
                end if
             end if
          end if
       end do
       
	   if (myid == master) then
       close(1)
       end if
       
       return
    end subroutine read_flow_file_parallel     

    subroutine read_flow_dual_parallel
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,pnb,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
       
       if (ndualtst == 2) then
	
	      if (myid == master) then
          open(1,file=flowname,status='unknown',form='unformatted')
          read(1)ttdts
          end if
          
          call MPI_BCAST(ttdts,1,mpi_reprec,master,MPI_COMM_WORLD,ierr)
             
          do nb=1,nblocks
             pid = mb_pids(nb)
          
             idim = mb_dim(nb,1) 
             jdim = mb_dim(nb,2) 
             kdim = mb_dim(nb,3) 
             
             if (master == pid-1) then
	            if (myid == master) then
                   read(1)((((mb_qmc(nb)%a4d(m,i,j,k),m=1,nl),i=1,idim),j=1,jdim),k=1,kdim)
                end if
             else
	            if (myid == master) then
                   allocate(mb_qmc(nb)%a4d(nl,idim,jdim,kdim),stat=ierr)
                   read(1)((((mb_qmc(nb)%a4d(m,i,j,k),m=1,nl),i=1,idim),j=1,jdim),k=1,kdim)
                   
                   packsize = idim*jdim*kdim*nl
                         
                   call MPI_SEND(mb_qmc(nb)%a4d,packsize,mpi_reprec, &
                                 pid-1,nb,MPI_COMM_WORLD,ierr)
                   
                   deallocate(mb_qmc(nb)%a4d,stat=ierr)
                endif
                
	            if (myid == pid-1) then
                   packsize = idim*jdim*kdim*nl
                   
                   call MPI_RECV(mb_qmc(nb)%a4d,packsize,mpi_reprec, &
                                 master,nb,MPI_COMM_WORLD,status,ierr)
                end if
             end if
          end do
       
	      if (myid == master) then
          close(1)
          end if
       
       else
       
          do pnb=1,pnblocks
             nb = pnbindexs(pnb)
       
             idim = mb_dim(nb,1) 
             jdim = mb_dim(nb,2) 
             kdim = mb_dim(nb,3) 
             
             do k=1,kdim
             do j=1,jdim
             do i=1,idim
                do m=1,nl
                   mb_qmc(nb)%a4d(m,i,j,k) = mb_qnc(nb)%a4d(m,i,j,k)
                end do
             end do
             end do
             end do
          end do
       
       end if
       
       return
    end subroutine read_flow_dual_parallel     


    subroutine write_flow_file_parallel
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
       real,pointer :: qpv_tmp(:,:,:,:)
	
	   if (myid == master) then
       open(1,file=flowname,status='unknown',form='unformatted')
       write(1)nout,wholetime,roo,uoo,voo,woo,poo,moo
       end if
       
       do nb=1,nblocks
          pid = mb_pids(nb)
       
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          if (master == pid-1) then
	         if (myid == master) then
                write(1)(((mb_r(nb)%a3d(i,j,k),mb_u(nb)%a3d(i,j,k),mb_v(nb)%a3d(i,j,k), &
                           mb_w(nb)%a3d(i,j,k),mb_p(nb)%a3d(i,j,k),mb_t(nb)%a3d(i,j,k), &
                           i=1,idim),j=1,jdim),k=1,kdim)
             end if
         
          else
	         if (myid == master) then
                allocate(qpv_tmp(6,idim,jdim,kdim),stat=ierr)
                
                packsize = idim*jdim*kdim*6
                      
                call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,status,ierr)
                              
                write(1)((((qpv_tmp(m,i,j,k),m=1,6),i=1,idim),j=1,jdim),k=1,kdim)
                
                deallocate(qpv_tmp,stat=ierr)
             endif
             
	         if (myid == pid-1) then
                allocate(qpv_tmp(6,idim,jdim,kdim),stat=ierr)
                
                do k=1,kdim
                do j=1,jdim
                do i=1,idim
                   qpv_tmp(1,i,j,k) = mb_r(nb)%a3d(i,j,k)
                   qpv_tmp(2,i,j,k) = mb_u(nb)%a3d(i,j,k)
                   qpv_tmp(3,i,j,k) = mb_v(nb)%a3d(i,j,k)
                   qpv_tmp(4,i,j,k) = mb_w(nb)%a3d(i,j,k)
                   qpv_tmp(5,i,j,k) = mb_p(nb)%a3d(i,j,k)
                   qpv_tmp(6,i,j,k) = mb_t(nb)%a3d(i,j,k)
                end do
                end do
                end do
               
                packsize = idim*jdim*kdim*6
                
                call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,ierr)
                
                deallocate(qpv_tmp,stat=ierr)
             end if
          end if
       end do
       
	   if (myid == master) then
       close(1)
       end if
       
       return
    end subroutine write_flow_file_parallel  

    subroutine write_flow_dual_parallel
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
	
	   if (myid == master) then
       open(1,file=trim(flowname)//".dts",status='unknown',form='unformatted')
       write(1)ttdts
       end if
       
       do nb=1,nblocks
          pid = mb_pids(nb)
       
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          if (master == pid-1) then
	         if (myid == master) then
                write(1)((((mb_qmc(nb)%a4d(m,i,j,k),m=1,nl),i=1,idim),j=1,jdim),k=1,kdim)
             end if
         
          else
	         if (myid == master) then
                allocate(mb_qmc(nb)%a4d(nl,idim,jdim,kdim),stat=ierr)
                
                packsize = idim*jdim*kdim*nl
                      
                call MPI_RECV(mb_qmc(nb)%a4d,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,status,ierr)
                              
                write(1)((((mb_qmc(nb)%a4d(m,i,j,k),m=1,nl),i=1,idim),j=1,jdim),k=1,kdim)
                
                deallocate(mb_qmc(nb)%a4d,stat=ierr)
             endif
             
	         if (myid == pid-1) then
                packsize = idim*jdim*kdim*nl
                
                call MPI_SEND(mb_qmc(nb)%a4d,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,ierr)
                
             end if
          end if
       end do
       
	   if (myid == master) then
       close(1)
       end if
       
       return
    end subroutine write_flow_dual_parallel  

    subroutine connect_pre_parallel
       use global_variables
       use mod_comp_connects
       implicit none
       integer :: nb,nbt,iused,ierr
       integer :: nr,nrt,nrmax,ibctype
       integer :: ibeg(3),iend(3),idir,idir2,idir3
       integer :: i,j,k,it,jt,kt,nbijk(4),nbijkt(4)
       integer,parameter :: mcyc(5) = (/1,2,3,1,2/)
       integer,pointer :: flags(:,:,:)
       
       do nb=1,nblocks
          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             mb_bc(nb)%bc(nr)%iused = 0
          end do
       end do
       
       do nb=1,nblocks
          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if ( ibctype < 0 ) then
                iused = mb_bc(nb)%bc(nr)%iused
                if (iused == 0) then
                   mb_bc(nb)%bc(nr)%iused = 1
                   nbt = mb_bc(nb)%bc(nr)%nbt
                   nrt = mb_bc(nb)%bc(nr)%ibcwin
                   mb_bc(nbt)%bc(nrt)%iused = -1
                end if
             end if
          end do
       end do
       
       ncomp_pnts = 0
       comp_list => null()
       do nb=1,nblocks
          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             ibctype = mb_bc(nb)%bc(nr)%bctype
             if ( ibctype < 0 ) then
                iused = mb_bc(nb)%bc(nr)%iused
                if (iused == 1) then
                   nbt = mb_bc(nb)%bc(nr)%nbt
                
                   ibeg = mb_bc(nb)%bc(nr)%s_st
                   iend = mb_bc(nb)%bc(nr)%s_ed
                   idir = mb_bc(nb)%bc(nr)%s_nd
                   
                   allocate(flags(ibeg(1):iend(1),ibeg(2):iend(2),ibeg(3):iend(3)),stat=ierr)
                   
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1) 
                      flags(i,j,k) = 1
                   end do
                   end do
                   end do
                   
                   idir2 = mcyc(idir+1)
                   idir3 = mcyc(idir+2)
                   ibeg(idir2) = ibeg(idir2) + 1
                   iend(idir2) = iend(idir2) - 1
                   ibeg(idir3) = ibeg(idir3) + 1
                   iend(idir3) = iend(idir3) - 1
                   
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1) 
                      flags(i,j,k) = 0
                   end do
                   end do
                   end do
                   
                   ibeg = mb_bc(nb)%bc(nr)%s_st
                   iend = mb_bc(nb)%bc(nr)%s_ed
                  
                   do k=ibeg(3),iend(3)
                   do j=ibeg(2),iend(2)
                   do i=ibeg(1),iend(1) 
                      if (flags(i,j,k) == 1) then
                         it = mb_bc(nb)%bc(nr)%image(i,j,k )
                         jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )
                         kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )
                         nbijk(1) = nb
                         nbijk(2) = i
                         nbijk(3) = j
                         nbijk(4) = k
                         nbijkt(1) = nbt
                         nbijkt(2) = it
                         nbijkt(3) = jt
                         nbijkt(4) = kt
                         
                         call comp_list_insert( comp_list,nbijk,nbijkt,ncomp_pnts)
                      end if
                   end do
                   end do
                   end do
                   
                   deallocate(flags,stat=ierr)
                end if
             end if
          end do
       end do
       
!!       if (myid == master) then
!!          if ( ncomp_pnts>0 ) call comp_list_print(comp_list,ncomp_pnts)
!!       end if

    end subroutine connect_pre_parallel
    
    subroutine calc_wall_dists_parallel
       use global_variables
       implicit none
       integer :: nwallcells
       real,pointer,dimension(:,:) :: coord_wallpnts
       integer :: nwcs_local(numprocs),nw_st_local(numprocs)
       integer :: nwcs_bc,nb,nr,nrmax,ibctype
       integer :: ibeg(3),iend(3),pnb,pid
       integer :: i,j,k,idim,jdim,kdim,m,n,ierr
       integer :: sendcount,recvcounts(numprocs),displs(numprocs)
       real :: xwall,ywall,zwall,xc,yc,zc,dx,dy,dz,dist

       nwallcells = 0
       nwcs_local = 0
       do nb=1,nblocks
          pid = mb_pids(nb)
          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             ibctype = mb_bc(nb)%bc(nr)%bctype
             select case(ibctype)
             case(2,20,21)
                ibeg = mb_bc(nb)%bc(nr)%s_st
                iend = mb_bc(nb)%bc(nr)%s_ed
                nwcs_bc = product(iend(:)-ibeg(:)+1,1)
                nwallcells = nwallcells + nwcs_bc
                nwcs_local(pid) = nwcs_local(pid) + nwcs_bc
             end select
          end do
       end do

       allocate(coord_wallpnts(3,nwallcells),stat=ierr)
       
       nw_st_local(1) = 0
       do i=2,numprocs
          nw_st_local(i) = nw_st_local(i-1) + nwcs_local(i-1)
       end do

       m = nw_st_local(myid+1)
       do pnb=1,pnblocks
          nb = pnbindexs(pnb)
          pid = mb_pids(nb)
          nrmax = mb_bc(nb)%nregions
          do nr = 1,nrmax
             ibctype = mb_bc(nb)%bc(nr)%bctype
             select case(ibctype)
             case(2,20,21)
                ibeg = mb_bc(nb)%bc(nr)%s_st
                iend = mb_bc(nb)%bc(nr)%s_ed
                
                do k=ibeg(3),iend(3)
                do j=ibeg(2),iend(2)
                do i=ibeg(1),iend(1) 
                   m = m + 1
                   coord_wallpnts(1,m) = mb_x(nb)%a3d(i,j,k)
                   coord_wallpnts(2,m) = mb_y(nb)%a3d(i,j,k)
                   coord_wallpnts(3,m) = mb_z(nb)%a3d(i,j,k)
                end do
                end do
                end do
             end select
          end do
       end do
       
       sendcount = nwcs_local(myid+1)*3
       recvcounts = nwcs_local(:)*3
       displs(:) = nw_st_local(:)*3
       
       n = nw_st_local(myid+1) + 1
       
       call MPI_ALLGATHERV(coord_wallpnts(1,n),sendcount,mpi_reprec, &
                           coord_wallpnts,recvcounts,displs,mpi_reprec,MPI_COMM_WORLD,ierr)


       do pnb=1,pnblocks
          nb = pnbindexs(pnb)
          
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          do k=1,kdim
          do j=1,jdim
          do i=1,idim
             mb_distance(nb)%a3d(i,j,k) = 1.0e30
          end do
          end do
          end do
       end do
       
       if (myid == master) then
          write(*,*) "  ======================================================="
          write(*,*) "           The distance to the nearest wall-point"
          write(*,*) "  -------------------------------------------------------"
          write(*,*) "   Total number of wall cells:",nwallcells
          write(*,*) "  ======================================================="
       end if
                                  
       do m=1,nwallcells
          xwall = coord_wallpnts(1,m)
          ywall = coord_wallpnts(2,m)
          zwall = coord_wallpnts(3,m)
       
          do pnb=1,pnblocks
             nb = pnbindexs(pnb)
             
             idim = mb_dim(nb,1) 
             jdim = mb_dim(nb,2) 
             kdim = mb_dim(nb,3) 
             
             do k=1,kdim
             do j=1,jdim
             do i=1,idim
                xc = mb_x(nb)%a3d(i,j,k)
                yc = mb_y(nb)%a3d(i,j,k)
                zc = mb_z(nb)%a3d(i,j,k)
                dx = xc - xwall
                dy = yc - ywall
                dz = zc - zwall
                dist = dx*dx + dy*dy + dz*dz
                
                mb_distance(nb)%a3d(i,j,k) = min(mb_distance(nb)%a3d(i,j,k),dist)
             end do
             end do
             end do
          end do
          
          if (myid == master) then
             if (mod(m,100) == 0) then
                write(*,*) m," of",nwallcells
             end if
          end if
       end do
              
       do pnb=1,pnblocks
          nb = pnbindexs(pnb)
          
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          do k=1,kdim
          do j=1,jdim
          do i=1,idim
             mb_distance(nb)%a3d(i,j,k) = sqrt(mb_distance(nb)%a3d(i,j,k)) + small
          end do
          end do
          end do
       end do
              
       deallocate(coord_wallpnts,stat=ierr)
       
       call write_wall_dists_parallel
       
       call send_message("Compute wall distance successfully!")
       
    end subroutine calc_wall_dists_parallel

    subroutine write_wall_dists_parallel  
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
       real,pointer :: dst_tmp(:,:,:)

	   if (myid == master) then
       open(1,file=trim(gridname)//".wdst",status='unknown',form='unformatted')
       end if
       
       do nb=1,nblocks
          pid = mb_pids(nb)
       
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          if (master == pid-1) then
	         if (myid == master) then
                write(1)(((mb_distance(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
             end if
         
          else
	         if (myid == master) then
                allocate(dst_tmp(idim,jdim,kdim),stat=ierr)
                
                packsize = idim*jdim*kdim
                      
                call MPI_RECV(dst_tmp,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,status,ierr)
                              
                write(1)(((dst_tmp(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
                
                deallocate(dst_tmp,stat=ierr)
             endif
             
	         if (myid == pid-1) then
                allocate(dst_tmp(idim,jdim,kdim),stat=ierr)
                
                do k=1,kdim
                do j=1,jdim
                do i=1,idim
                   dst_tmp(i,j,k) = mb_distance(nb)%a3d(i,j,k)
                end do
                end do
                end do
               
                packsize = idim*jdim*kdim
                
                call MPI_SEND(dst_tmp,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,ierr)
                
                deallocate(dst_tmp,stat=ierr)
             end if
          end if
       end do
       
	   if (myid == master) then
       close(1)
       end if
       
       return
    end subroutine write_wall_dists_parallel

    subroutine read_wall_dists_parallel
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
       real,pointer :: dst_tmp(:,:,:)
       
	   if (myid == master) then
       open(1,file=trim(gridname)//".wdst",status='old',form='unformatted')
       end if
       
       do nb=1,nblocks
          pid = mb_pids(nb)
       
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          if (master == pid-1) then
	         if (myid == master) then
                read(1)(((mb_distance(nb)%a3d(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
             end if
         
          else
	         if (myid == master) then
                allocate(dst_tmp(idim,jdim,kdim),stat=ierr)
                read(1)(((dst_tmp(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
                
                packsize = idim*jdim*kdim
                      
                call MPI_SEND(dst_tmp,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,ierr)
                
                deallocate(dst_tmp,stat=ierr)
             endif
             
	         if (myid == pid-1) then
                allocate(dst_tmp(idim,jdim,kdim),stat=ierr)
                
                packsize = idim*jdim*kdim
                
                call MPI_RECV(dst_tmp,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
                              
                do k=1,kdim
                do j=1,jdim
                do i=1,idim
                   mb_distance(nb)%a3d(i,j,k) = dst_tmp(i,j,k)
                end do
                end do
                end do
                
                deallocate(dst_tmp,stat=ierr)
             end if
          end if
       end do
       
	   if (myid == master) then
       close(1)
       end if
       
       call send_message("Read wall distance successfully!")
       
       return
    end subroutine read_wall_dists_parallel    

    subroutine read_turbulence_parallel
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
	   integer :: nvtur
       real,pointer :: qpv_tmp(:,:,:,:)
       
       nvtur = nlamtur + 1

	   if (myid == master) then
       open(1,file=trim(flowname)//".turb",status='old',form='unformatted')
       end if
       
       do nb=1,nblocks
          pid = mb_pids(nb)
       
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          if (master == pid-1) then
	         if (myid == master) then
                read(1)((((mb_qke(nb)%a4d(i,j,k,m),m=1,nlamtur),mb_vist(nb)%a3d(i,j,k), &
                           i=1,idim),j=1,jdim),k=1,kdim)
             end if
         
          else
	         if (myid == master) then
                allocate(qpv_tmp(nvtur,idim,jdim,kdim),stat=ierr)
                read(1)((((qpv_tmp(m,i,j,k),m=1,nvtur),i=1,idim),j=1,jdim),k=1,kdim)
                
                packsize = idim*jdim*kdim*nvtur
                      
                call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,ierr)
                
                deallocate(qpv_tmp,stat=ierr)
             endif
             
	         if (myid == pid-1) then
                allocate(qpv_tmp(nvtur,idim,jdim,kdim),stat=ierr)
                
                packsize = idim*jdim*kdim*nvtur
                
                call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
                              
                do k=1,kdim
                do j=1,jdim
                do i=1,idim
                   do m=1,nlamtur
                      mb_qke(nb)%a4d(i,j,k,m) = qpv_tmp(m,i,j,k)
                   end do
                   mb_vist(nb)%a3d(i,j,k) = qpv_tmp(nvtur,i,j,k)
                end do
                end do
                end do
                
                deallocate(qpv_tmp,stat=ierr)
             end if
          end if
       end do
       
	   if (myid == master) then
       close(1)
       end if
       
       return
    end subroutine read_turbulence_parallel     

    subroutine write_turbulence_parallel  
       use global_variables
       implicit none
       integer :: nb,i,j,k,m,pid,ierr
       integer :: idim,jdim,kdim,packsize
       integer :: status(MPI_STATUS_SIZE)
	   integer :: nvtur
       real,pointer :: qpv_tmp(:,:,:,:)
       
       nvtur = nlamtur + 1

	   if (myid == master) then
       open(1,file=trim(flowname)//".turb",status='unknown',form='unformatted')
       end if
       
       do nb=1,nblocks
          pid = mb_pids(nb)
       
          idim = mb_dim(nb,1) 
          jdim = mb_dim(nb,2) 
          kdim = mb_dim(nb,3) 
          
          if (master == pid-1) then
	         if (myid == master) then
                write(1)((((mb_qke(nb)%a4d(i,j,k,m),m=1,nlamtur),mb_vist(nb)%a3d(i,j,k), &
                           i=1,idim),j=1,jdim),k=1,kdim)
             end if
         
          else
	         if (myid == master) then
                allocate(qpv_tmp(nvtur,idim,jdim,kdim),stat=ierr)
                
                packsize = idim*jdim*kdim*nvtur
                      
                call MPI_RECV(qpv_tmp,packsize,mpi_reprec, &
                              pid-1,nb,MPI_COMM_WORLD,status,ierr)
                              
                write(1)((((qpv_tmp(m,i,j,k),m=1,nvtur),i=1,idim),j=1,jdim),k=1,kdim)
                
                deallocate(qpv_tmp,stat=ierr)
             endif
             
	         if (myid == pid-1) then
                allocate(qpv_tmp(nvtur,idim,jdim,kdim),stat=ierr)
                
                do k=1,kdim
                do j=1,jdim
                do i=1,idim
                   do m=1,nlamtur
                      qpv_tmp(m,i,j,k) = mb_qke(nb)%a4d(i,j,k,m)
                   end do
                   qpv_tmp(nvtur,i,j,k) = mb_vist(nb)%a3d(i,j,k)
                end do
                end do
                end do
               
                packsize = idim*jdim*kdim*nvtur
                
                call MPI_SEND(qpv_tmp,packsize,mpi_reprec, &
                              master,nb,MPI_COMM_WORLD,ierr)
                
                deallocate(qpv_tmp,stat=ierr)
             end if
          end if
       end do
       
	   if (myid == master) then
       close(1)
       end if
       
       return
    end subroutine write_turbulence_parallel
        
#endif
end module mod_parallels


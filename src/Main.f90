!_____________________________________________________________________!
program main
    use global_variables
    use global_const, only : method,nscheme
#ifdef PARALLEL
    use mod_parallels,only : master,myid,parallel_init,parallel_finalize, &
                             exchange_bc,communicate_dq_npp,exchange_bc_dq, &
                             communicate_pv_npp
#endif
    implicit none
    integer :: nb,nbg,ngs,ngstep,ngzone,ngz,m
    real    :: CPU_BEGIN,CPU_END
    real :: ajdu,bjdu
    real :: TA_MML(2)

    call parallel_init

!   ����Ԥ��������
!-------------------
    nijk2d  = 4    !*tgh ��ʱ�ӣ�<=nijk2d----��ά
    nijk2nd = 4    !*tgh ��ʱ�ӣ�<=nijk2nd----2��
!-------------------

    call pre_processing

    if ( nstart == 0 .or. nstart == 1 .or.  nstart == 3 ) then
       nout = 0
       wholetime = 0.0
       ajdu = attack*180.0/3.1415926
       bjdu = sideslip*180.0/3.1415926

    if (myid == master) then

       open(1,file=forcename,status='unknown')
       rewind(1)
       write(1,10)'     **************************************************************************'
       write(1,10)'     ****  Convegent History of Aerodynamics  ****'
       write(1,20)'     ****  ������    ��',moo   ,'  ��ŵ����',reynolds
       write(1,20)'     ****  ����      ��',ajdu,'  �⻬�ǣ�',bjdu
       write(1,20)'     ****  �����¶�  ��',tref,'  �����¶ȣ�',twall
       write(1,20)'     ****  �ο����  ��',sref ,'  �ο����ȣ�',lfref
       write(1,10)' ��������  ������  ������  ������  ��ת���� ƫ������ �������� ���� ���� ѹ��'
       write(1,10)'     **************************************************************************'
       write(1,30)'   variables=STEP,CPUTIME,cfx,cfy,cfz,cmx,cmr,cmz,cd,cl,cpc'
       close(1)

       open(1,file=errhis,status='unknown')
       rewind(1)
       write(1,10)'      *************************************************************************'
       write(1,10)'      ****  Convegent History of ERROR    ****'
       write(1,20)'      ****  ������    ��',moo   ,'  ��ŵ����',reynolds
       write(1,20)'      ****  ����      ��',ajdu,'  �⻬�ǣ�',bjdu
       write(1,20)'      ****  �����¶�  ��',tref,'  �����¶ȣ�',twall
       write(1,20)'      ****  �ο����  ��',sref ,'  �ο����ȣ�',lfref
       write(1,10)'  ��������  cfl��  ƽ���в�  ���в�  ���  ni nj nk nv(�������)  '
       write(1,10)'      *************************************************************************'
       write(1,30)'   variables=STEP,CPUTIME,cfl,AVER,MAX,NB,ni,nj,nk,nv'
       close(1)
10     format("#",a)
20     format("#",2(a,f14.6,1x))
30     format(a)
    end if

    endif

    nsubstmx = 1

if (myid == master) then
    write(*,*)'-----������ģ�飬�������㿪ʼ-----'
end if

    CALL CPU_TIME ( CPU_BEGIN )        !  CPU starts to operate

    nbg   =1
    ngzone=1
    ngstep=nomax

!*tgh modify ---begin ************************
    do ngs=1,nomax  !  nomax-�����ܲ���
       nout = nout + 1
       nstopsub = 0
       do nsubstep=1,nsubstmx
          phydtime=10000._prec

          call exchange_bc

          call boundary_all

          call initial_q_dq_c_t  !�����غ�����������������ٺ��¶�
          
	      call timestep_tgh   !*TGH. ����ʱ�䲽��*!

          if(nvis == 1)then
                call r_h_s_vis
          endif

          if( nvis ==1 )then !*tgh. ��ճʱճ������ɢ���ƽ��
            call communicate_dq_npp
          endif

          call r_h_s_invis    !* ��ճ����ɢ

        call communicate_dq_npp

        call l_h_s_tgh

        call communicate_pv_npp

        call residual_curve(nbg,ngzone)

       end do  ! *** ˫ʱ�䲽���� ***

       phytime=phytime+phydtime


!*tgh modify ---end ************************

       !   ���ú�������
       CALL CPU_TIME ( CPU_END )   !  CPU Stop in each block

       CPU_TIME_USED = CPU_END - CPU_BEGIN

       if ( nout/nwerror*nwerror == nout ) call saveres2file(nbg,ngzone)

       call post_processing(nout,nbg,ngzone,ngs,ngstep)

    enddo

#ifdef PARALLEL
    call parallel_finalize
#endif

end program main

!_____________________________________________________________________!
subroutine spectrum
    use global_variables
    implicit none

    call spectinv
    if ( nvis == 1  ) then
       call compute_visl
       if ( csrv > 0.0 ) call spectvis
    endif

    return
end subroutine spectrum
!_____________________________________________________________________!

subroutine get_c(nchem,nchem_source)
    implicit none
    integer :: nchem,nchem_source
          call get_c_t
    return
end subroutine get_c
!_____________________________________________________________________!
subroutine compute_visl
    use global_variables
    implicit none
       call compute_visl_ns
    return
end subroutine compute_visl

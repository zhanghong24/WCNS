subroutine initial_q_dq_c_t
   use global_variables, only : nblocks,nvis
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

		  call recast_field(nb)

		  if(nvis == 1)then
             call get_c_t_only !*TGH. �õ����ٺ��¶�
             call compute_visl_ns  !! Added by Fox Liu.
		  else
             call get_c_only !*TGH. �õ�����
		  endif

	      call get_q_only !�õ��غ����
          call set_dq_to_0


    enddo


end subroutine initial_q_dq_c_t

!=====================================================================
subroutine initial_t !��ʼ���ڳ����¶�
   use global_variables, only : nblocks,nvis
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
		  call recast_field(nb)

		  if(nvis == 1)then
             call get_t_ini
		  endif

    enddo

end subroutine initial_t

!=====================================================================
!=====================================================================
subroutine recast_gcl(nb) !*tgh. ӭ���ʽ+�����غ���ʱʹ��
    use global_variables,only : mb_kcx_l,mb_kcy_l,mb_kcz_l,mb_kct_l,&
                                mb_etx_l,mb_ety_l,mb_etz_l,mb_ett_l,&
                                mb_ctx_l,mb_cty_l,mb_ctz_l,mb_ctt_l,&
                                mb_vol_l,mb_volt_l,                 &
                                mb_kcx_r,mb_kcy_r,mb_kcz_r,mb_kct_r,&
                                mb_etx_r,mb_ety_r,mb_etz_r,mb_ett_r,&
                                mb_ctx_r,mb_cty_r,mb_ctz_r,mb_ctt_r,&
                                mb_vol_r,mb_volt_r,                 &
								vol_l ,volt_l,              &
                                kcx_l ,kcy_l ,kcz_l ,kct_l ,&
                                etx_l ,ety_l ,etz_l ,ett_l ,&
                                ctx_l ,cty_l ,ctz_l ,ctt_l ,&
								vol_r ,volt_r,              &
                                kcx_r ,kcy_r ,kcz_r ,kct_r ,&
                                etx_r ,ety_r ,etz_r ,ett_r ,&
                                ctx_r ,cty_r ,ctz_r ,ctt_r ,&
								gcl   ,nscheme
    implicit none
    integer :: nb

    vol_l => mb_vol_l(nb)%a3d
 !   volt_l => mb_volt_l(nb)%a3d   !*tgh. ע�⣺�Ƕ�������ʱ��ʹ��

    kcx_l => mb_kcx_l(nb)%a3d
    kcy_l => mb_kcy_l(nb)%a3d
    kcz_l => mb_kcz_l(nb)%a3d
    kct_l => mb_kct_l(nb)%a3d

    etx_l => mb_etx_l(nb)%a3d
    ety_l => mb_ety_l(nb)%a3d
    etz_l => mb_etz_l(nb)%a3d
    ett_l => mb_ett_l(nb)%a3d

    ctx_l => mb_ctx_l(nb)%a3d
    cty_l => mb_cty_l(nb)%a3d
    ctz_l => mb_ctz_l(nb)%a3d
    ctt_l => mb_ctt_l(nb)%a3d

    vol_r => mb_vol_r(nb)%a3d
 !   volt_r => mb_volt_r(nb)%a3d   !*tgh. ע�⣺�Ƕ�������ʱ��ʹ��

    kcx_r => mb_kcx_r(nb)%a3d
    kcy_r => mb_kcy_r(nb)%a3d
    kcz_r => mb_kcz_r(nb)%a3d
    kct_r => mb_kct_r(nb)%a3d

    etx_r => mb_etx_r(nb)%a3d
    ety_r => mb_ety_r(nb)%a3d
    etz_r => mb_etz_r(nb)%a3d
    ett_r => mb_ett_r(nb)%a3d

    ctx_r => mb_ctx_r(nb)%a3d
    cty_r => mb_cty_r(nb)%a3d
    ctz_r => mb_ctz_r(nb)%a3d
    ctt_r => mb_ctt_r(nb)%a3d

    return
end subroutine recast_gcl

!=============================================================================!
!=============================================================================!
subroutine boundary_all !*tgh. ����߽紦��
    use define_precision_mod
    use global_variables,only:nblocks
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
    integer :: nb,pnb

!!#ifdef PARALLEL
!!    do pnb=1,pnblocks
!!       nb = pnbindexs(pnb)
!!#else
!!    do nb=1,nblocks
!!#endif
!!!       call recast_method(nb)  !*TGH. ע��˳���ܸı�
!!        call recast_grid(nb)
!!        call recast_field(nb)
!!	    call boundary(nb)
!!	enddo



!!#ifdef PARALLEL
!!    do pnb=1,pnblocks
!!       nb = pnbindexs(pnb)
!!#else
!!    do nb=1,nblocks
!!#endif
!!!       call recast_method(nb)  !*TGH. ע��˳���ܸı�
!!        call recast_grid(nb)
!!        call recast_field(nb) 
!!	    call boundary_a(nb)
!!
!!	enddo
!!
!!#ifdef PARALLEL
!!    do pnb=1,pnblocks
!!       nb = pnbindexs(pnb)
!!#else
!!    do nb=1,nblocks
!!#endif
!!!       call recast_method(nb)  !*TGH. ע��˳���ܸı�
!!        call recast_grid(nb)
!!        call recast_field(nb) 
!!	    call boundary_b(nb)
!!
!!	enddo
	
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
!       call recast_method(nb)  !*TGH. ע��˳���ܸı�
        call recast_grid(nb)
        call recast_field(nb) 
	    call boundary_sequence(nb)

	enddo
	
	
end subroutine boundary_all

!=============================================================================!
!=============================================================================!

subroutine bc_connect_vis
!*TGH. ���ʱ����ճ����(dq)�ڶԽӱ߽�����ƽ��
!*TGH.  ���ǻ����ʺϻ�ѧ��Ӧ�����

	use define_precision_mod
    use global_variables, only : mb_bc,mb_dq,nblocks,nl,mb_vol
    implicit none
    integer :: nb,nr,nrmax,bctype,m

    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks
    integer :: it,jt,kt
    integer :: s_st(3),s_ed(3)
	real(prec) :: vols,volt

	do nb=1,nblocks

       nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����

	   do nr=1,nrmax
		  bctype=mb_bc(nb)%bc(nr)%bctype

	     if( bctype < 0 )then !*�����Խӱ߽�

             do m=1,3
                s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
                s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
             enddo

              nbs   = mb_bc(nb)%bc(nr)%nbs             !Դ���

		       IF(NBS /= NB)THEN
			      write(*,*)'==============�߽�Խӿ��ܳ���============'
		      ENDIF

              nbt   = mb_bc(nb)%bc(nr)%nbt     !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������

              do is = s_st(1),s_ed(1)
              do js = s_st(2),s_ed(2)
              do ks = s_st(3),s_ed(3)
                it = mb_bc(nb)%bc(nr)%image(is,js,ks )
                jt = mb_bc(nb)%bc(nr)%jmage(is,js,ks )
                kt = mb_bc(nb)%bc(nr)%kmage(is,js,ks )

				vols = mb_vol(nbs)%a3d(is,js,ks)
				volt = mb_vol(nbt)%a3d(it,jt,kt)

				do m=2,nl  !*TGH. ��һ��ճ�Է���Ϊ0
                  mb_dq(nbs)%a4d(m,is,js,ks) = 0.5_prec*( mb_dq(nbs)%a4d(m,is,js,ks)/vols &
				                                     &  + mb_dq(nbt)%a4d(m,it,jt,kt)/volt )
				  mb_dq(nbt)%a4d(m,it,jt,kt) =	mb_dq(nbs)%a4d(m,is,js,ks)*volt
				  mb_dq(nbs)%a4d(m,is,js,ks) =	mb_dq(nbs)%a4d(m,is,js,ks)*vols
                enddo
              enddo
              enddo
              enddo

		  endif  !*END. �����Խӱ߽�

	   enddo !* nr=1,nrmax

	enddo !* nb=1,nblocks

    return

end subroutine bc_connect_vis

!=============================================================================!
!=============================================================================!

subroutine bc_connect_cic
	use define_precision_mod
    use global_variables, only : mb_bc,mb_dq,nblocks,nl,method, &
	                           & cic1,cic2,cic3,cic4,cic5
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs,bc_n1_cic_pre_parallel,dif_cic_new_parallel
#endif
    implicit none

    integer :: nb,nr,nrmax,bctype,m,pnb

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif

        call recast_grid(nb)
        call recast_field(nb)

       nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����

	   do nr=1,nrmax
          bctype = mb_bc(nb)%bc(nr)%bctype
          if( bctype < 0 ) then !�Խӱ߽�
#ifdef PARALLEL
              call bc_n1_cic_pre_parallel(nb,nr,bctype)
#else
		      call bc_n1_cic_pre(nb,nr,bctype) !*TGH. �����Խ�Ԥ����
		                  !*TGH. ��Ŀ����Ϊ�˰ѶԽӵ�Ŀ����ϵ�rhs��dq���浽Դ��������
#endif
         endif

	   enddo !* nr=1,nrmax

	enddo !*nb = 1,nblocks
!*tgh. �����ǶԽ�ǰ����
!*tgh. �����������Խ�
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif

        call recast_grid(nb)
        call recast_field(nb)

       nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����

	   do nr=1,nrmax
          bctype = mb_bc(nb)%bc(nr)%bctype
          if( bctype < 0 ) then
!!             call dif_cic_tgh(nb,nr)  !*tgh. �����Խ�
#ifdef PARALLEL
             call dif_cic_new_parallel(nb,nr)
#else
             call dif_cic_new(nb,nr)  !*tgh. �����Խ�
#endif

          end if
!*TGH. �����Խ��޸�dq����-rhs����������ϵ�ֵ��¼�Խ����϶�Ӧ��Ŀ����dq
	   enddo !* nr=1,nrmax

	enddo !*nb = 1,nblocks
!*tgh. �����������Խ�
!*tgh. �����ǰѷǶԽ����ϵ�dq����

!* �Ѿ��ĵ�LHS�У�����������Recastһ������
!	do nb=1,nblocks
!        call recast_field(nb)
!        nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
!
!	    do nr=1,nrmax
!           bctype = mb_bc(nb)%bc(nr)%bctype
!	          if( bctype > 0)  then  !�ǶԽӱ߽磬���������桢Զ����
!                 call boundary_dq_to_0(nb,nr,bctype)
!			  ENDIF
!		ENDDO
!	ENDDO


end subroutine bc_connect_cic

!=============================================================================!
!=============================================================================!

subroutine bc_n1_cic_pre(nb,nr,bctype)
!*TGH. �Խӱ߽����� ǰ���������������Խӱ߽� ֻ�������޲��
!*TGH. ��Ŀ����ڶԽ����ϵ�dq���ֵ�Դ��������
    use global_variables,only : mb_bc,mb_dq,nblocks,method,nl !&
!	                          &,mb_r,mb_u,mb_v,mb_w,mb_p
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
!                nt = n - 1 + method
!                it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)*nt
!                jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)*nt !*TGH. Ŀ���棬���ʱ��������һ�㣬���ʱ���䣨����������
!                kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)*nt
!                mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)
!                mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)  !*TGH. Դ��������һ��
!                mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. Ŀ�����������һ��
!                mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)  !*TGH. ��Դ������ȡĿ�������Ӧ���ֵ
!                mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)

                do n1=1,5
					mb_dq(nbs)%a4d(n1,is,js,ks) = mb_dq(nbt)%a4d(n1,it0,jt0,kt0)
				!*TGH. ע�⣬��Դ������ֵ����¼Ŀ������ڶԽ����ϵ�RHS
				!*TGH. ���ǵ���Ҫʹ��ʱ�����Բ�����Ŀ������Ϣ��ֱ�Ӳ�������ϵ�ֵ����
				enddo

             enddo
          enddo
       enddo
    enddo
end subroutine bc_n1_cic_pre
!=============================================================================!
subroutine bc_n1_cic_post(nb,nr,bctype)
!*TGH. �����Խӱ߽�ʱ������ֻ�������޲��
!*TGH. Ŀ���ǰѶԽ����ϵ�dq����ƽ��
	use define_precision_mod
    use global_variables,only : mb_bc,mb_dq,nblocks,method,nl
    implicit none
    integer :: nb,nr,bctype,m,n,n1
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: i,j,k,nt,ntarg
    integer :: it,jt,kt
    integer :: s_st(3),s_ed(3)

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
    enddo
!    s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
!    s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!    s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
    nbs   = mb_bc(nb)%bc(nr)%nbs             !���

    nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
!    t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
!    t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!    t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             do n=1,1
                it= mb_bc(nb)%bc(nr)%image(i,j,k )
                jt= mb_bc(nb)%bc(nr)%jmage(i,j,k )   !*TGH. Ŀ���棬���ʱ�ڶԽ����ϣ����ʱ��������һ��
                kt= mb_bc(nb)%bc(nr)%kmage(i,j,k )

                do n1=1,5
					mb_dq(nbs)%a4d(n1,i,j,k) = 0.5_prec*( mb_dq(nbt)%a4d(n1,it,jt,kt) + &
					                                    & mb_dq(nbs)%a4d(n1,i,j,k)    )

				enddo

             enddo
          enddo
       enddo
    enddo
end subroutine bc_n1_cic_post
!=============================================================================!
!=============================================================================!
subroutine dif_cic_tgh(nb,nr) !�����Խ�����
!*TGH. �����Խ��޸�rhs����dq��������ϵ�ֵ��¼�Խ����϶�Ӧ��Ŀ����rhs
!*TGH. ��Ҫ�õ��Խ����Լ�����ϵ�rhs1�����������õ��Խ����ϵ�dq��
!*tgh. �����Խ���Ȼ�ü�ƽ��
!*TGH. ʹ��ǰ����Recast������������

	use define_precision_mod
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
	real  :: nx,ny,nz,nt,nvols,nvolt
	real  :: ARRAY_P1(5,5),ARRAY_P(5,5),ARRAY_Asp(5,5),ARRAY_Asn(5,5),ary_c(5,5)
    real  :: rm,um,vm,wm,pm,mm,cm,tm,cgm,ct,qbar1(5)
    real    :: l1,l4,l5
    real  :: qsp(5),qsn(5),qspn(5)
    integer :: is1,js1,ks1,is2,js2,ks2
    integer :: s01,s02
!	  real  :: nxt,nyt,nzt,ntt !,nxt1,nyt1,nzt1,ntt1
	  real  :: dqq(nl)
	  real  :: nxyz(4)

!	  integer :: LWARN
!	integer :: it,jt,kt,it1,jt1,kt1
!    real  :: pstar
!    real  :: rm1,um1,vm1,wm1,pm1,l1t,l4t,l5t,f(nl)


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

!             it1 = mb_bc(nb)%bc(nr)%image(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(1)
!             jt1 = mb_bc(nb)%bc(nr)%jmage(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(2) !*TGH. �Խ��棨Ŀ��飩��������һ�������
!             kt1 = mb_bc(nb)%bc(nr)%kmage(i,j,k ) - mb_bc(nb)%bc(nr)%t_lr3d(3)

             if ( nsurf == 1 ) then    !*tgh. ���ݶԽ����������ϵķ���ʸ��
                nx = kcx(i,j,k) ! mb_kcx(nb)%a3d(i,j,k)
                ny = kcy(i,j,k) ! mb_kcy(nb)%a3d(i,j,k)
                nz = kcz(i,j,k) ! mb_kcz(nb)%a3d(i,j,k)  !*tgh. Դ���ϵķ���ʸ��
                nt = kct(i,j,k) ! mb_kct(nb)%a3d(i,j,k)

!                nxt1 = mb_kcx(nbt)%a3d(it1,jt1,kt1)
!                nyt1 = mb_kcy(nbt)%a3d(it1,jt1,kt1)
!                nzt1 = mb_kcz(nbt)%a3d(it1,jt1,kt1) !*tgh. Ŀ�����������һ���ķ���ʸ��
!                ntt1 = mb_kct(nbt)%a3d(it1,jt1,kt1)
             elseif ( nsurf == 2 ) then
                nx = etx(i,j,k) ! mb_etx(nb)%a3d(i,j,k)
                ny = ety(i,j,k) ! mb_ety(nb)%a3d(i,j,k)
                nz = etz(i,j,k) ! mb_etz(nb)%a3d(i,j,k)
                nt = ett(i,j,k) ! mb_ett(nb)%a3d(i,j,k)

!                nxt1 = mb_etx(nbt)%a3d(it1,jt1,kt1)
!                nyt1 = mb_ety(nbt)%a3d(it1,jt1,kt1)
!                nzt1 = mb_etz(nbt)%a3d(it1,jt1,kt1)
!                ntt1 = mb_ett(nbt)%a3d(it1,jt1,kt1)
             else
                nx = ctx(i,j,k) ! mb_ctx(nb)%a3d(i,j,k)
                ny = cty(i,j,k) ! mb_cty(nb)%a3d(i,j,k)
                nz = ctz(i,j,k) ! mb_ctz(nb)%a3d(i,j,k)
                nt = ctt(i,j,k) ! mb_ctt(nb)%a3d(i,j,k)

!                nxt1 = mb_ctx(nbt)%a3d(it1,jt1,kt1)
!                nyt1 = mb_cty(nbt)%a3d(it1,jt1,kt1)
!                nzt1 = mb_ctz(nbt)%a3d(it1,jt1,kt1)
!                ntt1 = mb_ctt(nbt)%a3d(it1,jt1,kt1)
             endif

			 nvols = mb_vol(nbs)%a3d(i,j,k)
			 nvolt = mb_vol(nbt)%a3d(it,jt,kt)


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
!			    dqq(n1) = mb_dq(nb)%a4d(n1,i,j,k)    !*TGH. Դ���ϵĶԽ���
			    dqq(n1) = dq(n1,i,j,k)    !*TGH. Դ���ϵĶԽ���
			 enddo
			 call  Aps_dq(qbar1,nxyz,s_lr,dqq,qsp,1)  !*TGH. ��A~��A+������ʾ��Ȼ��ֱ�Ӷ�RHS����ת��
!TGH. ����Ϊ ԭʼ���������������Խӷ���rhs�����������ϵ�rhs��

			 do n1 = 1,5
			    dqq(n1) = dq(n1,is1,js1,ks1)*nvols/nvolt  !*TGH. Դ��������һ���RHS
			 enddo

			 call  Aps_dq(qbar1,nxyz,s_lr,dqq,qsn,-1)
!____

!_____
			 do n1 = 1,5
!			     mb_dq(nb)%a4d(n1,i,j,k) = (qsp(n1)+qsn(n1))    !*TGH. ??? ʹ���������RHS?��dq?
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


			 if(rm<1.e-20 .or. pm<1.e-20 .or. rm>1.e10 .or. pm>1.e10 )then
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
end subroutine dif_cic_tgh



subroutine  post_processing_tgh
    use define_precision_mod
    use global_variables,only:nblocks,nout,nomax
    implicit none
    integer :: ngzone,nbg,ngs,ngstep


	ngzone = 1
	nbg    = 1
	ngs    = nout
	ngstep = nomax


	call post_processing(nout,nbg,ngzone,ngs,ngstep)
	     !* TGH. nout-ÿ���Ѿ��������ܲ���, nbg-��ǰ���, ngzone-������, ����ʱngzone=1 *!
	     !* TGH. ngs-ÿ�鵱ǰ��������,ngstep-ÿ���趨�ĵ����ܲ��� *!



end subroutine  post_processing_tgh

!=============================================================================!
!=============================================================================!

subroutine BC_farfield_turbulence_tgh(nb,is,js,ks,it,jt,kt,lr3d,s_nd)
!*TGH. ����ģ�ͱ߽�
!*TGH. ��ڲ���Զ��ֵ
!*TGH. ���ڲ�����һ���ֵ

  use define_precision_mod
  use global_variables,only : mb_qke,nlamtur,nameturb,koo,goo,muoo,omgoo &
                            &,mb_u,mb_v,mb_w,      mb_kcx,mb_kcy,mb_kcz  &
                            &,mb_etx,mb_ety,mb_etz,mb_ctx,mb_cty,mb_ctz

  implicit none

	integer :: lr3d(3),s_nd
	integer :: nb,is,js,ks,it,jt,kt,m
	real :: ut,tem_lr

	tem_lr=real(lr3d(s_nd))
	
	if(s_nd ==1 )then
	   ut = tem_lr* ( mb_kcx(nb)%a3d(is,js,ks) * mb_u(nb)%a3d(is,js,ks) + &
	                & mb_kcy(nb)%a3d(is,js,ks) * mb_v(nb)%a3d(is,js,ks) + &
	                & mb_kcz(nb)%a3d(is,js,ks) * mb_w(nb)%a3d(is,js,ks) )
	elseif(s_nd == 2)then
	   ut = tem_lr* ( mb_etx(nb)%a3d(is,js,ks) * mb_u(nb)%a3d(is,js,ks) + &
	            & mb_ety(nb)%a3d(is,js,ks) * mb_v(nb)%a3d(is,js,ks) + &
	            & mb_etz(nb)%a3d(is,js,ks) * mb_w(nb)%a3d(is,js,ks) )
	elseif(s_nd ==3)then
	   ut = tem_lr* ( mb_ctx(nb)%a3d(is,js,ks) * mb_u(nb)%a3d(is,js,ks) + &
	                & mb_cty(nb)%a3d(is,js,ks) * mb_v(nb)%a3d(is,js,ks) + &
	                & mb_ctz(nb)%a3d(is,js,ks) * mb_w(nb)%a3d(is,js,ks) )
	else
	   write(*,*)"�ڴ���������Զ���߽�ʱ����"
	   write(*,*)"����������Զ���ķ���s_nd����"
	   stop
	endif
	
	do m=0,1
	  is = is + lr3d(1) * m
	  js = js + lr3d(2) * m
	  ks = ks + lr3d(3) * m

	  if(ut < 0._prec )then !*tgh. ������
	  
	
	   if(nlamtur>0)then
		 if(nameturb=='SA')then
			mb_qke(nb)%a4d(is,js,ks,1) = muoo
		 elseif(nameturb=='SST')then
			mb_qke(nb)%a4d(is,js,ks,1) = koo
			mb_qke(nb)%a4d(is,js,ks,2) = omgoo
		 elseif(nameturb=='KG')then
			mb_qke(nb)%a4d(is,js,ks,1) = koo
			mb_qke(nb)%a4d(is,js,ks,2) = goo
		 endif
	   endif

	  else !*tgh. �������

	   if(nlamtur>0)then
		 if(nameturb=='SA')then
			mb_qke(nb)%a4d(is,js,ks,1) = mb_qke(nb)%a4d(it,jt,kt,1)
		 elseif(nameturb=='SST')then
			mb_qke(nb)%a4d(is,js,ks,1) = mb_qke(nb)%a4d(it,jt,kt,1)
			mb_qke(nb)%a4d(is,js,ks,2) = mb_qke(nb)%a4d(it,jt,kt,2)
		 elseif(nameturb=='KG')then
			mb_qke(nb)%a4d(is,js,ks,1) = mb_qke(nb)%a4d(it,jt,kt,1)
			mb_qke(nb)%a4d(is,js,ks,2) = mb_qke(nb)%a4d(it,jt,kt,2)
		 endif
	   endif

	  endif
	enddo

	return
end subroutine BC_farfield_turbulence_tgh


!=============================================================================!
!=============================================================================!
subroutine DUVWT_halfNODE_TGH(nmax,ni,n1,uvwt,duvwt)
!*tgh. ֱ�Ӵӽڵ�ֵ�����ڵ��ϵ�������
!*tgh. ��������ϵ�ֵ
	use global_variables,only : nijk2nd
    implicit none

	real,parameter :: A1=27.d0,B1=-1.d0
	real,parameter :: A2=-22.d0,B2=17.d0,C2=9.d0,D2=-5.d0,E2=1.d0  !*TGH. OLD
!	real,parameter :: A2=-23.D0,B2=21.D0,C2=3.D0,D2=-1.D0,E2=0.D0  !*TGH NEW WCNS_E__BORDER
	real,parameter :: A3=-71.D0,B3=141.D0,C3=-93.D0,D3=23.D0,E3=0.D0  !*TGH NEW WCNS_E_5_BORDER

	integer :: nmax,i,m,ni,n1,n2
	real    :: uvwt(n1,-1:nmax+1),duvwt(n1,0:nmax)

	if( ni <= nijk2nd )then !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
!
! deal with the symmetric boundary condiction
		do m=1,n1
			do i=0,ni
				duvwt(m,i) = uvwt(m,i+1) - uvwt(m,i)
			enddo
		enddo
	else
		do m=1,n1
			do i=2,ni-2
				duvwt(m,i) = (A1*(uvwt(m,i+1) - uvwt(m,i)) + &
										  B1*(uvwt(m,i+2) - uvwt(m,i-1)))/24.d0
			enddo

			duvwt(m,1   ) =  (A2*uvwt(m,1)  + B2*uvwt(m,2)    + C2*uvwt(m,3)    + D2*uvwt(m,4)    + E2*uvwt(m,5)   )/24.d0
			duvwt(m,ni-1) = -(A2*uvwt(m,ni) + B2*uvwt(m,ni-1) + C2*uvwt(m,ni-2) + D2*uvwt(m,ni-3) + E2*uvwt(m,ni-4))/24.d0

			duvwt(m,0   ) =  (A3*uvwt(m,1)  + B3*uvwt(m,2)    + C3*uvwt(m,3)   + D3*uvwt(m,4)   + E3*uvwt(m,5)   )/24.d0 !*tgh. �������
			duvwt(m,ni  ) = -(A3*uvwt(m,ni) + B3*uvwt(m,ni-1) + C3*uvwt(m,ni-2)+ D3*uvwt(m,ni-3)+ E3*uvwt(m,ni-4))/24.d0 !*tgh. �������

		enddo

	endif
!
  return
end subroutine DUVWT_halfNODE_TGH
!=============================================================================!
!=============================================================================!
subroutine vis_highorder_tgh
!*TGH. ���ø߽׸�ʽ��ɢճ����
!*TGH. ������������dq��
!
!*TGH. Modified by TU Guohua, 2009.2
!
!*TGH. ��ĳ�������������С�ڵ���6����رո÷������ɢ����ά��
!*TGH. ���򲻹��Ż���Ϊ�˽�Լ�ڴ棬���������ظ����㣬��������!
!--------------------------------------------------------------
   use define_precision_mod
   use global_variables
   use stress_module
   use geometry_module
   use duvwt_module

    implicit none

    integer :: i,j,k,m
    real :: ev_l(nl)
    real,pointer,dimension( :,: ) :: fv,fv1

    real :: vis
    real :: cp,kcp
    real :: vis_l,vis_t,co_rdi
    real :: dtdx,dtdy,dtdz
    real :: um,vm,wm,qx,qy,qz
    real :: cs_init(ns),xi(ns),rdi(ns),hs(ns)
    real :: dfidx(ns),dfidy(ns),dfidz(ns)
    real :: ev(nl)
    real :: temp
    allocate( fv(nl,nmax),fv1(nl,0:nmax) )

    re = 1.0 / reynolds
    cp = 1.0/((gama-1.0)*moo*moo)

!----------------------------------------------
! Modified by TU Guohua
! old
!    do k=2,nk-1
!       do j=2,nj-1
!
! Modified
    if(ni <= nijk2d)goto 11  !*TGH. ��ά
     do k=1,nk
     do j=1,nj
! end Modified
!----------------------------------------------
          !�����һ���������ϵ�ճ��ͨ��


          do i =1,ni
			 vis_l = visl(i,j,k)
			 vis_t = vist(i,j,k)
			 vis = vis_l + vis_t

             if ( vol(i,j,k) > 0.0 ) then
			    call getgeo(i,j,k)     !* TGH. ȡ��������
			    call getuvwtder_TGH_4th(i,j,k)  !*TGH. �������ڼ���ռ�ĵ���
			    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                               dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)  !*TGH. ת���������ռ�ĵ���
				call stress(vis) !*TGH. ��Ӧ����txx,txy,txz,tyy,tyz,tzz

				kcp = ( vis_l/prl + vis_t/prt ) * cp

				qx =  kcp * dtdx
				qy =  kcp * dtdy
				qz =  kcp * dtdz

    !����ط�Ҫע��,��Ҫ������
				um = u(i,j,k)
				vm = v(i,j,k)
				wm = w(i,j,k)

!--  ���ӳ������ճ��ͨ����δ������ŵ��
                call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                       kx,ky,kz,qx,qy,qz,um,vm,wm,ev,nl)

             else
                ev(:) = 0.0
             endif

             do m=2,nl
                fv(m,i) = ev(m)
             enddo

          enddo !* i =1,ni

	      call VALUE_HALF_NODE(4,nmax,ni,fv(2:5,1:nmax),fv1(2:5,0:nmax))
		  call flux_dxyz(4,nmax,ni,fv1(2:5,0:nmax),fv(2:5,1:nmax))

          do i=1,ni
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) - re * fv(m,i)
             enddo
          enddo

       enddo
    enddo
11  continue

!----------------------------------------------
! Modified by TU Guohua
! old
!    do k=2,nk-1
!       do i=2,ni-1
!
! Modified
    if(nj <= nijk2d)goto 12  !*TGH. ��ά
    do k=1,nk
       do i=1,ni
! end Modified
!----------------------------------------------
          !�����һ���������ϵ�ճ��ͨ��


          do j = 1,nj
			 vis_l = visl(i,j,k)
			 vis_t = vist(i,j,k)
			 vis = vis_l + vis_t

             if ( vol(i,j,k) > 0.0 ) then
			    call getgeo(i,j,k)     !* TGH. ȡ��������
			    call getuvwtder_TGH_4th(i,j,k)  !*TGH. �������ڼ���ռ�ĵ���
			    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                               dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)  !*TGH. ת���������ռ�ĵ���
				call stress(vis) !*TGH. ��Ӧ����txx,txy,txz,tyy,tyz,tzz

				kcp = ( vis_l/prl + vis_t/prt ) * cp

				qx =  kcp * dtdx
				qy =  kcp * dtdy
				qz =  kcp * dtdz

    !����ط�Ҫע��,��Ҫ������
				um = u(i,j,k)
				vm = v(i,j,k)
				wm = w(i,j,k)

!--  ���ӳ������ճ��ͨ����δ������ŵ��
                call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                       ex,ey,ez,qx,qy,qz,um,vm,wm,ev,nl)

             else
                ev(:) = 0.0
             endif

             do m=2,nl
                fv(m,j) = ev(m)
             enddo
          enddo ! j=1,nj

	      call VALUE_HALF_NODE(4,nmax,nj,fv(2:5,1:nmax),fv1(2:5,0:nmax))
		  call flux_dxyz(4,nmax,nj,fv1(2:5,0:nmax),fv(2:5,1:nmax))


          do j=1,nj
             do m=2,nl
                dq(m,i,j,k) = dq(m,i,j,k) -  re * fv(m,j)
             enddo
          enddo

       enddo
    enddo
12  continue

!----------------------------------------------
! Modified by TU Guohua
! old
!    do i=2,ni-1
!       do j=2,nj-1
!
! Modified
    if(nk <= nijk2d)goto 13  !*TGH. ��ά

    do i=1,ni
       do j=1,nj
! end Modified
!----------------------------------------------
          !�����һ���������ϵ�ճ��ͨ��

          do k = 1,nk
			 vis_l = visl(i,j,k)
			 vis_t = vist(i,j,k)
			 vis = vis_l + vis_t

             if ( vol(i,j,k) > 0.0 ) then

 			    call getgeo(i,j,k)     !* TGH. ȡ��������
			    call getuvwtder_TGH_4th(i,j,k)  !*TGH. �������ڼ���ռ�ĵ���
			    call getdfdxyz(dudx,dudy,dudz,dvdx,dvdy,dvdz, &
                               dwdx,dwdy,dwdz,dtdx,dtdy,dtdz)  !*TGH. ת���������ռ�ĵ���
				call stress(vis) !*TGH. ��Ӧ����txx,txy,txz,tyy,tyz,tzz

				kcp = ( vis_l/prl + vis_t/prt ) * cp

				qx =  kcp * dtdx
				qy =  kcp * dtdy
				qz =  kcp * dtdz

    !����ط�Ҫע��,��Ҫ������
				um = u(i,j,k)
				vm = v(i,j,k)
				wm = w(i,j,k)

!--  ���ӳ������ճ��ͨ����δ������ŵ��
                call fluxvcd_point(txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz, &
                       cx,cy,cz,qx,qy,qz,um,vm,wm,ev,nl)

             else
                ev(:) = 0.0
             endif

             do m=2,nl
                fv(m,k) = ev(m)
             enddo

          enddo ! k=1,nk

	      call VALUE_HALF_NODE(4,nmax,nk,fv(2:5,1:nmax),fv1(2:5,0:nmax))
		  call flux_dxyz(4,nmax,nk,fv1(2:5,0:nmax),fv(2:5,1:nmax))


          do k=1,nk
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) - re* fv(m,k)
             enddo
          enddo


      enddo
    enddo
13  continue

    deallocate( fv,fv1 )

    return
end subroutine vis_highorder_tgh
!=============================================================================!
!=============================================================================!
subroutine getuvwtder_tgh_4th(i,j,k)

    use define_precision_mod
    use global_variables
    use duvwt_module
    implicit none
    real,parameter :: ac01= 8._prec/12._prec, ac02= -1._prec/12._prec ,ac03=0._prec !*tgh. �ڵ�
    real,parameter :: ac11=-11._prec/6._prec,  ac12=18._prec/6._prec  ,ac13=-9._prec/6._prec &
	                &,ac14=  2._prec/6._prec,  ac15=0._prec     !*tgh. �߽�
    real,parameter :: ac21= -2._prec/6._prec,  ac22=-3._prec/6._prec  ,ac23= 6._prec/6._prec &
	                &,ac24=  -1._prec/6._prec,  ac25=0._prec     !*tgh. �α߽�

    integer :: i,j,k
    real :: aa,bb,cc
    real :: u1,v1,w1,t1,u2,v2,w2,t2,u3,v3,w3,t3
    real :: AC1,AC2,AC3,AC4,AC5,u4,v4,w4,t4,u5,v5,w5,t5

	if(ni <= nijk2d)then !*TGH. ��ά
       ukc =  0._prec
       vkc =  0._prec
       wkc =  0._prec
       tkc =  0._prec

    elseif ( i == 1 ) then
       call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
       call getuvwt( 2 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( 3 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       call getuvwt( 4 ,j ,k ,u4 ,v4 ,w4 ,t4 )
       ukc =  ac11 * u1 + ac12 * u2 + ac13 * u3 + ac14 * u4
       vkc =  ac11 * v1 + ac12 * v2 + ac13 * v3 + ac14 * v4
       wkc =  ac11 * w1 + ac12 * w2 + ac13 * w3 + ac14 * w4
       tkc =  ac11 * t1 + ac12 * t2 + ac13 * t3 + ac14 * t4


    elseif ( i == ni ) then
       call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
       call getuvwt( i-1 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( i-2 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       call getuvwt( i-3 ,j ,k ,u4 ,v4 ,w4 ,t4 )
       ukc =  -( ac11 * u1 + ac12 * u2 + ac13 * u3 + ac14 * u4 )
       vkc =  -( ac11 * v1 + ac12 * v2 + ac13 * v3 + ac14 * v4 )
       wkc =  -( ac11 * w1 + ac12 * w2 + ac13 * w3 + ac14 * w4 )
       tkc =  -( ac11 * t1 + ac12 * t2 + ac13 * t3 + ac14 * t4 )


    elseif( i == 2 ) THEN
       call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
       call getuvwt( 1 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( 3 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       call getuvwt( 4 ,j ,k ,u4 ,v4 ,w4 ,t4 )
       ukc =  ac21 * u2 + ac22 * u1 + ac23 * u3 + ac24 * u4
       vkc =  ac21 * v2 + ac22 * v1 + ac23 * v3 + ac24 * v4
       wkc =  ac21 * w2 + ac22 * w1 + ac23 * w3 + ac24 * w4
       tkc =  ac21 * t2 + ac22 * t1 + ac23 * t3 + ac24 * t4

    elseif( i == ni-1 ) THEN
       call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
       call getuvwt( i+1 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( i-1 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       call getuvwt( i-2 ,j ,k ,u4 ,v4 ,w4 ,t4 )
       ukc = -( ac21 * u2 + ac22 * u1 + ac23 * u3 + ac24 * u4 )
       vkc = -( ac21 * v2 + ac22 * v1 + ac23 * v3 + ac24 * v4 )
       wkc = -( ac21 * w2 + ac22 * w1 + ac23 * w3 + ac24 * w4 )
       tkc = -( ac21 * t2 + ac22 * t1 + ac23 * t3 + ac24 * t4 )
	else

       call getuvwt( i-2 ,j ,k ,u1 ,v1 ,w1 ,t1 )
       call getuvwt( i-1 ,j ,k ,u2 ,v2 ,w2 ,t2 )
       call getuvwt( i+1 ,j ,k ,u3 ,v3 ,w3 ,t3 )
       call getuvwt( i+2 ,j ,k ,u4 ,v4 ,w4 ,t4 )
       ukc = ac01 * ( u3 - u2) + ac02 * ( u4 - u1 )
       vkc = ac01 * ( v3 - v2) + ac02 * ( v4 - v1 )
       wkc = ac01 * ( w3 - w2) + ac02 * ( w4 - w1 )
       tkc = ac01 * ( t3 - t2) + ac02 * ( t4 - t1 )
    endif

!*tgh. ����i����
!*tgh. ����j����

    if(nj <= nijk2d)then !*TGH. ��ά
       uet =  0._prec
       vet =  0._prec
       wet =  0._prec
       tet =  0._prec

	elseif ( j == 1 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,2 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,3 ,k ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,4 ,k ,u4 ,v4 ,w4 ,t4 )
       uet =  ac11 * u1 + ac12 * u2 + ac13 * u3 + ac14 * u4
       vet =  ac11 * v1 + ac12 * v2 + ac13 * v3 + ac14 * v4
       wet =  ac11 * w1 + ac12 * w2 + ac13 * w3 + ac14 * w4
       tet =  ac11 * t1 + ac12 * t2 + ac13 * t3 + ac14 * t4

    elseif ( j == nj ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,j-1 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j-2 ,k ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j-3 ,k ,u4 ,v4 ,w4 ,t4 )
       uet =  -( ac11 * u1 + ac12 * u2 + ac13 * u3 + ac14 * u4 )
       vet =  -( ac11 * v1 + ac12 * v2 + ac13 * v3 + ac14 * v4 )
       wet =  -( ac11 * w1 + ac12 * w2 + ac13 * w3 + ac14 * w4 )
       tet =  -( ac11 * t1 + ac12 * t2 + ac13 * t3 + ac14 * t4 )

    elseif ( j == 2 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,1 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,3 ,k ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,4 ,k ,u4 ,v4 ,w4 ,t4 )
       uet =  ac21 * u2 + ac22 * u1 + ac23 * u3 + ac24 * u4
       vet =  ac21 * v2 + ac22 * v1 + ac23 * v3 + ac24 * v4
       wet =  ac21 * w2 + ac22 * w1 + ac23 * w3 + ac24 * w4
       tet =  ac21 * t2 + ac22 * t1 + ac23 * t3 + ac24 * t4

    elseif ( j == nj-1 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,j+1 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j-1 ,k ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j-2 ,k ,u4 ,v4 ,w4 ,t4 )
       uet = -( ac21 * u2 + ac22 * u1 + ac23 * u3 + ac24 * u4 )
       vet = -( ac21 * v2 + ac22 * v1 + ac23 * v3 + ac24 * v4 )
       wet = -( ac21 * w2 + ac22 * w1 + ac23 * w3 + ac24 * w4 )
       tet = -( ac21 * t2 + ac22 * t1 + ac23 * t3 + ac24 * t4 )


    else
        call getuvwt( i ,j-2 ,k ,u1 ,v1 ,w1 ,t1 )
        call getuvwt( i ,j-1 ,k ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j+1 ,k ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j+2 ,k ,u4 ,v4 ,w4 ,t4 )
       uet = ac01 * ( u3 - u2) + ac02 * ( u4 - u1 )
       vet = ac01 * ( v3 - v2) + ac02 * ( v4 - v1 )
       wet = ac01 * ( w3 - w2) + ac02 * ( w4 - w1 )
       tet = ac01 * ( t3 - t2) + ac02 * ( t4 - t1 )
    endif

!*tgh. ����j����
!*tgh. ����k����
	if(nk <= nijk2d)then !*TGH. ��ά
       uct =  0._prec
       vct =  0._prec
       wct =  0._prec
       tct =  0._prec
    elseif ( k == 1 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,j ,2 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,3 ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j ,4 ,u4 ,v4 ,w4 ,t4 )
       uct =  ac11 * u1 + ac12 * u2 + ac13 * u3 + ac14 * u4
       vct =  ac11 * v1 + ac12 * v2 + ac13 * v3 + ac14 * v4
       wct =  ac11 * w1 + ac12 * w2 + ac13 * w3 + ac14 * w4
       tct =  ac11 * t1 + ac12 * t2 + ac13 * t3 + ac14 * t4
    elseif ( k == nk ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,j ,k-1 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,k-2 ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j ,k-3 ,u4 ,v4 ,w4 ,t4 )
       uct =  -( ac11 * u1 + ac12 * u2 + ac13 * u3 + ac14 * u4 )
       vct =  -( ac11 * v1 + ac12 * v2 + ac13 * v3 + ac14 * v4 )
       wct =  -( ac11 * w1 + ac12 * w2 + ac13 * w3 + ac14 * w4 )
       tct =  -( ac11 * t1 + ac12 * t2 + ac13 * t3 + ac14 * t4 )
    elseif ( k == 2 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,j ,1 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,3 ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j ,4 ,u4 ,v4 ,w4 ,t4 )
       uct =  ac21 * u2 + ac22 * u1 + ac23 * u3 + ac24 * u4
       vct =  ac21 * v2 + ac22 * v1 + ac23 * v3 + ac24 * v4
       wct =  ac21 * w2 + ac22 * w1 + ac23 * w3 + ac24 * w4
       tct =  ac21 * t2 + ac22 * t1 + ac23 * t3 + ac24 * t4
    elseif ( k == nk-1 ) then
        call getuvwt( i ,j ,k ,u1 ,v1 ,w1 ,t1)
        call getuvwt( i ,j ,k+1 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,k-1 ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j ,k-2 ,u4 ,v4 ,w4 ,t4 )
       uct = -( ac21 * u2 + ac22 * u1 + ac23 * u3 + ac24 * u4 )
       vct = -( ac21 * v2 + ac22 * v1 + ac23 * v3 + ac24 * v4 )
       wct = -( ac21 * w2 + ac22 * w1 + ac23 * w3 + ac24 * w4 )
       tct = -( ac21 * t2 + ac22 * t1 + ac23 * t3 + ac24 * t4 )

    else
        call getuvwt( i ,j ,k-2 ,u1 ,v1 ,w1 ,t1 )
        call getuvwt( i ,j ,k-1 ,u2 ,v2 ,w2 ,t2 )
        call getuvwt( i ,j ,k+1 ,u3 ,v3 ,w3 ,t3 )
        call getuvwt( i ,j ,k+2 ,u4 ,v4 ,w4 ,t4 )
       uct = ac01 * ( u3 - u2) + ac02 * ( u4 - u1 )
       vct = ac01 * ( v3 - v2) + ac02 * ( v4 - v1 )
       wct = ac01 * ( w3 - w2) + ac02 * ( w4 - w1 )
       tct = ac01 * ( t3 - t2) + ac02 * ( t4 - t1 )
    endif

    return
end subroutine getuvwtder_tgh_4th
!=============================================================================!
!=============================================================================!
subroutine timestep_tgh
   use define_precision_mod
   use global_variables,only : ntmst,nblocks,dtdt,vol,sml_vol &
                             &,mb_dtdt,mb_vol,mb_dim,nvis,PHYDTIME
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
   implicit none
   real(prec) :: dtmin,dtime_tem
   integer :: nb,i,j,k,ni,nj,nk,pnb

	dtmin=1000000000000.0
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
        call recast_GRID(nb)
        call recast_field(nb)

		call spectrum_tgh

       	call localdt0

	enddo

	if(ntmst == 1)then  !*TGH. ȫ��ʱ�䲽��ʱȡ��Сʱ�䲽
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
		 dtdt => mb_dtdt(nb)%a3d
		 vol  => mb_vol(nb)%a3d
		 ni = mb_dim(nb,1)
		 nj = mb_dim(nb,2)
		 nk = mb_dim(nb,3)

          do k=1,nk
          do j=1,nj
          do i=1,ni
             dtdt(i,j,k) = dtmin / max(vol(i,j,k), sml_vol)
          enddo
          enddo
         enddo
	  enddo
	elseif(.true.)then !tgh. Ѱ������ʱ�䲽���е���Сֵ
#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif
       dtdt => mb_dtdt(nb)%a3d
       vol  => mb_vol(nb)%a3d
       ni = mb_dim(nb,1)
       nj = mb_dim(nb,2)
       nk = mb_dim(nb,3)
       
       DTMIN=1000.0
       do k=2,nk-1
       do j=2,nj-1
       do i=2,ni-1
          dtmin = min(dtmin, dtdt(i,j,k)*vol(i,j,k))
       end do
       end do
       end do
    end do

	else !tgh. ȡ���һ����������ĵ��ʱ�䲽��Ϊ����ֵ
	  dtmin=dtdt(NI/2,NJ/2,NK/2)*vol(NI/2,NJ/2,NK/2)
	endif

	phydtime=dtmin !tgh. ����ʱ�䲽������Сʱ�䲽��

   return
end subroutine timestep_tgh

!=============================================================================!
!=============================================================================!
subroutine globaldt_tgh(ttmin)
!*tgh. ����CFLȷ����Сʱ�䲽��

    use define_precision_mod
    use global_variables,only : method,ni,nj,nk,sra,srb,src,srva,srvb,srvc,vol,cfl,sml_vol
    implicit none
    integer :: i,j,k
    real(prec) :: ra,rb,rc,rv,rabc,ttmin

	ttmin  = 100.0
    do k=2,nk-1
       do j=2,nj-1
          do i=2,ni-1
             ra = sra(i,j,k)
             rb = srb(i,j,k)
             rc = src(i,j,k)
             rv = srva(i,j,k) + srvb(i,j,k) + srvc(i,j,k)
!             rabc =( ra + rb + rc )/(vol(i,j,k)+sml_vol)
             rabc =( ra + rb + rc + rv )/max(vol(i,j,k),sml_vol) !  
             ttmin = min ( cfl/rabc , ttmin )
          enddo
       enddo
    enddo

end subroutine globaldt_tgh
!=============================================================================!
!=============================================================================!

subroutine post_for_2d
!* �ڶ�ά�����չ�����ƽ��
!* չ���ٶ�W=0���ұ���ΪW

    use define_precision_mod
    use global_variables,only : nblocks,r,u,v,w,p,t,ni,nj,nk,nijk2d
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
	real(prec) :: ra,ua,va,wa,pa
	integer :: nb,i,j,k,pnb

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
       call recast_field(nb)

	   if( ni <= nijk2d )then !*TGH. ��ά
		 do k=1,nk
	     do j=1,nj
			 ra = 0._prec
			 ua = 0._prec
			 va = 0._prec
			 wa = 0._prec
			 pa = 0._prec
			 do i=1,ni
				ra = ra + r(i,j,k)
				ua = ua + u(i,j,k)
				va = va + v(i,j,k)
!				wa = wa + w(i,j,k)
				pa = pa + p(i,j,k)
			 enddo
			 ra = ra/real(ni)
			 ua = ua/real(ni)
			 va = va/real(ni)
			 wa = 0._prec
			 pa = pa/real(ni)

			 do i=1,ni
			     r(i,j,k) = ra
			     u(i,j,k) = ua
			     v(i,j,k) = va
			     w(i,j,k) = wa
			     p(i,j,k) = pa
			 enddo

	     enddo
	     enddo

	   elseif(nj <= nijk2d)then !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��

		 do k=1,nk
	     do i=1,ni
			 ra = 0._prec
			 ua = 0._prec
			 va = 0._prec
			 wa = 0._prec
			 pa = 0._prec
			 do j=1,nj
				ra = ra + r(i,j,k)
				ua = ua + u(i,j,k)
				va = va + v(i,j,k)
!				wa = wa + w(i,j,k)
				pa = pa + p(i,j,k)
			 enddo
			 ra = ra/real(nj)
			 ua = ua/real(nj)
			 va = va/real(nj)
			 wa = 0._prec
			 pa = pa/real(nj)

			 do j=1,nj
			     r(i,j,k) = ra
			     u(i,j,k) = ua
			     v(i,j,k) = va
			     w(i,j,k) = wa
			     p(i,j,k) = pa
			 enddo

	     enddo
	     enddo


	   elseif(nk <= nijk2d)then !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��

		 do j=1,nj
	     do i=1,ni
			 ra = 0._prec
			 ua = 0._prec
			 va = 0._prec
			 wa = 0._prec
			 pa = 0._prec
			 do k=1,nk
				ra = ra + r(i,j,k)
				ua = ua + u(i,j,k)
				va = va + v(i,j,k)
!				wa = wa + w(i,j,k)
				pa = pa + p(i,j,k)
			 enddo
			 ra = ra/real(nk)
			 ua = ua/real(nk)
			 va = va/real(nk)
			 wa = 0._prec
			 pa = pa/real(nk)

			 do k=1,nk
			     r(i,j,k) = ra
			     u(i,j,k) = ua
			     v(i,j,k) = va
			     w(i,j,k) = wa
			     p(i,j,k) = pa
			 enddo

	     enddo
	     enddo

	   endif


	enddo


end subroutine post_for_2d

!=============================================================================!
!=============================================================================!

subroutine WCNSE5_VIS_virtual
  use define_precision_mod
  use global_variables,only : u,v,w,t,r,reynolds,visl,vist,dq,nvis     &
                            , kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,vol &
							, nl,nmax,ni,nj,nk,gama,moo,prl,prt
  use duvwt_all_field
  implicit none
!-----------------------------------------------------------------------------!
!	��ԭ WCNSE5_VIS Ҫ�õ���һ���ֵ                                          !
!	                                                                          !
!                                                                             !
!		��ƣ�Ϳ����                                                          !
!		���ԣ�Ϳ���� 2009.03                                                  !
!-----------------------------------------------------------------------------!
	integer :: i,j,k,m
	real(prec) :: duvwtdxyz (12,0:nmax),nxyz_line(3,0:nmax)
	real(prec) :: duvwt_line(12,nmax) !,uvwt_line(4,nmax)
	real(prec) :: duvwt_half(12,0:nmax),uvwt_half(4,0:nmax)
	real(prec) :: kxyz_line ( 9,nmax),kxyz_half(9,0:nmax)
	real(prec) :: vol_line(nmax),fv(nl,0:nmax),dfv(nl,nmax)
	real(prec) :: vol_half(0:nmax) !,txyz_half(9,0:nmax)
	real(prec) :: vslt1_half(0:nmax),vslt2_half(0:nmax),cp,cp_prl,cp_prt
	real(prec) :: vslt1_line(nmax),vslt2_line(nmax),re
	integer :: mvist
    real(prec) :: uvwt_lvir(4,-2:nmax+3) !Ϊ�˱���ʹ������ϵ�ֵ

	allocate(duvwt(12,ni,nj,nk))   ! overcome stack over problem
	allocate(duvwt_mid(12,0:ni,0:nj,0:nk))   !___ ֱ�Ӽ������1�׵���  2009.2.1

	re = 1.0/reynolds

	cp = 1.0/((gama-1.0)*moo*moo)
    cp_prl = cp/prl
    cp_prt = cp/prt

! to calculate the deriative of u,v,w,t in computing	coordinate

	call UVWT_DER_4th_virtual      !*TGH. ֱ����ڵ�1�׵��������߽��3����ܶ�Ϊ����Ҫ�����ֵ�����򣬲������
    call UVWT_DER_4th_half_virtual !*TGH. ֱ�������1�׵��������߽��3����ܶ�Ϊ����Ҫ�����ֵ�����򣬲������

!**************** I direction ****************
	do k=1,nk
	do j=1,nj

		do i=1,ni

!			uvwt_line(1,i) = u(i,j,k)
!			uvwt_line(2,i) = v(i,j,k)
!			uvwt_line(3,i) = w(i,j,k)
!			uvwt_line(4,i) = t(i,j,k)

			vslt1_line (i) = visl(i,j,k)
			vslt2_line (i) = visl(i,j,k)*cp_prl
		    do mvist=1,nvis !*����ʱ
			  vslt1_line (i) = vslt1_line(i) + vist(i,j,k)
			  vslt2_line (i) = vslt2_line(i) + vist(i,j,k)*cp_prt
			enddo

			kxyz_line(1,i) = kcx(i,j,k)
			kxyz_line(2,i) = etx(i,j,k)
			kxyz_line(3,i) = ctx(i,j,k)

			kxyz_line(4,i) = kcy(i,j,k)
			kxyz_line(5,i) = ety(i,j,k)
			kxyz_line(6,i) = cty(i,j,k)

			kxyz_line(7,i) = kcz(i,j,k)
			kxyz_line(8,i) = etz(i,j,k)
			kxyz_line(9,i) = ctz(i,j,k)

			vol_line(   i) = vol(i,j,k)

!			do m=1,12
			do m=5,12
				duvwt_line(m,i) = duvwt(m,i,j,k)
			enddo

		enddo

		do i=-2,ni+3
			uvwt_lvir(1,i) = u(i,j,k)
			uvwt_lvir(2,i) = v(i,j,k)
			uvwt_lvir(3,i) = w(i,j,k)
			uvwt_lvir(4,i) = t(i,j,k)
		enddo

        if( r(-2,j,k)   < 0.0 ) uvwt_lvir(4,-2)   = -1.0
        if( r(ni+3,j,k) < 0.0 ) uvwt_lvir(4,ni+3) = -1.0


		call VALUE_line_half_virtual(4,nmax,ni,uvwt_lvir,uvwt_half ) !��ԭʼ�����ڰ����ϵ�ֵ
		call VALUE_HALF_NODE_IJK(0,1,1,12,nmax,ni ,duvwt_line,duvwt_half) !�ѽڵ��ϵ�һ�׵�����ֵ������
		call VALUE_HALF_NODE(9 ,nmax,ni,kxyz_line ,kxyz_half )   !����������ֵ������
		call VALUE_HALF_NODE(1 ,nmax,ni,vol_line  ,vol_half  )
		call VALUE_HALF_NODE(1 ,nmax,ni,vslt1_line,vslt1_half  ) !*���꣺1:ni --> 0:ni
		call VALUE_HALF_NODE(1 ,nmax,ni,vslt2_line,vslt2_half  )

		do i=0,ni
		do m=1,4
			duvwt_half(m,i) = duvwt_mid(m,i,j,k)    !___ͬ����
		enddo
		enddo

		call DUVWT_DXYZ(nmax,ni,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)
!*TGH. ����ռ�����ϵ�һ�׵�������������ռ�ĵ���, ������ֵ��Χ0:nmax

		do i=0,ni
			nxyz_line(1,i) = kxyz_half(1,i)   !�������ȶ��������ռ��е�һ�׵���
			nxyz_line(2,i) = kxyz_half(4,i)
			nxyz_line(3,i) = kxyz_half(7,i)
		enddo


		call FLUX_VIS_LINE_new(nl,nmax,ni,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    )
		  !*TGH. ���߼���ճ��ͨ���ڼ���ռ�����ϵ�ͨ��

		call FLUX_DXYZ(nl,nmax,ni,fv,dfv)

		do i = 1,ni
		do m=1,nl
			dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,i)
		enddo
		enddo

	enddo
	enddo
!
!**************** J direction ****************
	do k=1,nk
	do i=1,ni
		do j=1,nj

!			uvwt_line(1,j) = u(i,j,k)
!			uvwt_line(2,j) = v(i,j,k)
!			uvwt_line(3,j) = w(i,j,k)
!			uvwt_line(4,j) = t(i,j,k)

			vslt1_line (j) = visl(i,j,k)
			vslt2_line (j) = visl(i,j,k)*cp_prl
		    do mvist=1,nvis !*����ʱ
			  vslt1_line (j) = vslt1_line(j) + vist(i,j,k)
			  vslt2_line (j) = vslt2_line(j) + vist(i,j,k)*cp_prt
			enddo

			kxyz_line(1,j) = kcx(i,j,k)
			kxyz_line(2,j) = etx(i,j,k)
			kxyz_line(3,j) = ctx(i,j,k)

			kxyz_line(4,j) = kcy(i,j,k)
			kxyz_line(5,j) = ety(i,j,k)
			kxyz_line(6,j) = cty(i,j,k)

			kxyz_line(7,j) = kcz(i,j,k)
			kxyz_line(8,j) = etz(i,j,k)
			kxyz_line(9,j) = ctz(i,j,k)

			vol_line(   j) = vol(i,j,k)

			do m=1,4
				duvwt_line(m,j) = duvwt(m,i,j,k)
			enddo
			do m=9,12
				duvwt_line(m,j) = duvwt(m,i,j,k)
			enddo

		enddo

		do j=-2,nj+3
			uvwt_lvir(1,j) = u(i,j,k)
			uvwt_lvir(2,j) = v(i,j,k)
			uvwt_lvir(3,j) = w(i,j,k)
			uvwt_lvir(4,j) = t(i,j,k)
		enddo

        if( r(i,-2,k)   < 0.0 ) uvwt_lvir(4,-2)   = -1.0
        if( r(i,nj+3,k) < 0.0 ) uvwt_lvir(4,nj+3) = -1.0

		call VALUE_line_half_virtual(4,nmax,nj,uvwt_lvir,uvwt_half )
		call VALUE_HALF_NODE_IJK(1,0,1,12,nmax,nj,duvwt_line,duvwt_half)
		call VALUE_HALF_NODE(9 ,nmax,nj,kxyz_line ,kxyz_half )
		call VALUE_HALF_NODE(1 ,nmax,nj,vol_line  ,vol_half  )

		call VALUE_HALF_NODE(1 ,nmax,nj,vslt1_line,vslt1_half  ) !*���꣺1:nj --> 0:nj
		call VALUE_HALF_NODE(1 ,nmax,nj,vslt2_line,vslt2_half  )

		do j=0,nj
		do m=5,8
			duvwt_half(m,j) = duvwt_mid(m,i,j,k)    !___ͬ����
		enddo
		enddo

		call DUVWT_DXYZ(nmax,nj,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)

		do j=0,nj
			nxyz_line(1,j) = kxyz_half(2,j)
			nxyz_line(2,j) = kxyz_half(5,j)
			nxyz_line(3,j) = kxyz_half(8,j)
		enddo

		call  FLUX_VIS_LINE_new(nl,nmax,nj,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    )

		call FLUX_DXYZ(nl,nmax,nj,fv,dfv)

		do j = 1,nj
		do m=1,nl
			dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,j)
		enddo
		enddo

	enddo
	enddo
!
!**************** K direction ****************
	do j=1,nj
	do i=1,ni
		do k=1,nk

!			uvwt_line(1,k) = u(i,j,k)
!			uvwt_line(2,k) = v(i,j,k)
!			uvwt_line(3,k) = w(i,j,k)
!			uvwt_line(4,k) = t(i,j,k)

			vslt1_line (k) = visl(i,j,k)
			vslt2_line (k) = visl(i,j,k)*cp_prl
		    do mvist=1,nvis !*����ʱ
			  vslt1_line (k) = vslt1_line(k) + vist(i,j,k)
			  vslt2_line (k) = vslt2_line(k) + vist(i,j,k)*cp_prt
			enddo

			kxyz_line(1,k) = kcx(i,j,k)
			kxyz_line(2,k) = etx(i,j,k)
			kxyz_line(3,k) = ctx(i,j,k)

			kxyz_line(4,k) = kcy(i,j,k)
			kxyz_line(5,k) = ety(i,j,k)
			kxyz_line(6,k) = cty(i,j,k)

			kxyz_line(7,k) = kcz(i,j,k)
			kxyz_line(8,k) = etz(i,j,k)
			kxyz_line(9,k) = ctz(i,j,k)

			vol_line(   k) = vol(i,j,k)

!			do m=1,12
			do m=1,8
				duvwt_line(m,k) = duvwt(m,i,j,k)
			enddo
		enddo

		do k=-2,nk+3
			uvwt_lvir(1,k) = u(i,j,k)
			uvwt_lvir(2,k) = v(i,j,k)
			uvwt_lvir(3,k) = w(i,j,k)
			uvwt_lvir(4,k) = t(i,j,k)
		enddo

        if( r(i,j,-2)   < 0.0 ) uvwt_lvir(4,-2)   = -1.0
        if( r(i,j,nk+3) < 0.0 ) uvwt_lvir(4,nk+3) = -1.0

		call VALUE_line_half_virtual(4,nmax,nk,uvwt_lvir,uvwt_half )
		call VALUE_HALF_NODE_IJK(1,1,0,12,nmax,nk,duvwt_line,duvwt_half)
		call VALUE_HALF_NODE(9 ,nmax,nk,kxyz_line ,kxyz_half )
		call VALUE_HALF_NODE(1 ,nmax,nk,vol_line  ,vol_half  )
		call VALUE_HALF_NODE(1 ,nmax,nk,vslt1_line,vslt1_half  ) !*���꣺1:nk --> 0:nk
		call VALUE_HALF_NODE(1 ,nmax,nk,vslt2_line,vslt2_half  )

		do k=0,nk
		do m=9,12
			duvwt_half(m,k) = duvwt_mid(m,i,j,k)    !___ͬ����
		enddo
		enddo

		call DUVWT_DXYZ(nmax,nk,12,duvwt_half,9,kxyz_half,vol_half,duvwtdxyz)

		do k=0,nk
			nxyz_line(1,k) = kxyz_half(3,k)
			nxyz_line(2,k) = kxyz_half(6,k)
			nxyz_line(3,k) = kxyz_half(9,k)
		enddo

		call FLUX_VIS_LINE_new(nl,nmax,nk,4,uvwt_half,3,nxyz_line, &
			                    12,duvwtdxyz,vslt1_half,vslt2_half,fv    )

		call FLUX_DXYZ(nl,nmax,nk,fv,dfv)

		do k = 1,nk
		do m=1,nl
			dq(m,i,j,k) = dq(m,i,j,k) - re*dfv(m,k)
		enddo
		enddo

	enddo
	enddo
!
	deallocate(duvwt)
	deallocate(duvwt_mid)
  return
end subroutine WCNSE5_VIS_virtual
!=============================================================================!
!=============================================================================!

subroutine UVWT_DER_4th
  use global_variables,only:u,v,w,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt
  implicit none
	integer :: i,j,k,m
	real    :: uvwt(4,-1:nmax+1),dtem(4,nmax)
!	real    :: duvwtdkc(nv,ni,nj,nk)

! I direction

	do k=1,nk
		do j=1,nj
			do i=0,ni+1

				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)

			enddo

			call DUVWT_NODE_LINE(nmax,ni,4,uvwt,dtem) !*TGH. ��ʱû�������ֵ
!			call DUVWT_NODE_LINE_virtual(nmax,ni,4,uvwt,dtem) !*TGH. Ҫ�����ֵ

			do i=1,ni
				do m=1,4

					duvwt(m,i,j,k)=dtem(m,i)

				enddo
			enddo

		enddo
	enddo
!
! J direction

	do k=1,nk
		do i=1,ni
			do j=0,nj+1

				uvwt(1,j) = u(i,j,k)
				uvwt(2,j) = v(i,j,k)
				uvwt(3,j) = w(i,j,k)
				uvwt(4,j) = t(i,j,k)

			enddo

			call DUVWT_NODE_LINE(nmax,nj,4,uvwt,dtem)
!			call DUVWT_NODE_LINE_virtual(nmax,nj,4,uvwt,dtem)

			do j=1,nj
				do m=1,4

					duvwt(m+4,i,j,k)=dtem(m,j)

				enddo
			enddo

		enddo
	enddo

! K direction

	do j=1,nj
		do i=1,ni
			do k=0,nk+1

				uvwt(1,k) = u(i,j,k)
				uvwt(2,k) = v(i,j,k)
				uvwt(3,k) = w(i,j,k)
				uvwt(4,k) = t(i,j,k)

			enddo

			call DUVWT_NODE_LINE(nmax,nk,4,uvwt,dtem)
!			call DUVWT_NODE_LINE_virtual(nmax,nk,4,uvwt,dtem)

			do k=1,nk
				do m=1,4

					duvwt(m+8,i,j,k)=dtem(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_4th

!=============================================================================!
!=============================================================================!

subroutine UVWT_DER_4th_virtual
!! ֱ����ڵ�1�׵��������߽��3����ܶ�Ϊ����Ҫ�����ֵ�����򣬲������
  use global_variables,only:u,v,w,r,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt
  implicit none
	integer :: i,j,k,m,st,ed
	real    :: uvwt(4,-2:nmax+3),dtem(4,nmax)
!	real    :: duvwtdkc(nv,ni,nj,nk)

! I direction

	do k=1,nk
		do j=1,nj

            if( r(-2,j,k)<0.0 )then
                uvwt(4,-2) = -1.0
                st = 0
            else
                st = -2
            endif

            if( r(ni+3,j,k)<0.0 )then
                uvwt(4,ni+3) = -1.0
                ed = ni+1
            else
                ed = ni+3
            endif

			do i=st,ed
				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)
			enddo


			call DUVWT_NODE_LINE_virtual(nmax,ni,4,uvwt,dtem) !*TGH. Ҫ�����ֵ

			do i=1,ni
				do m=1,4
					duvwt(m,i,j,k)=dtem(m,i)
				enddo
			enddo

		enddo
	enddo
!
! J direction

	do k=1,nk
		do i=1,ni
            if( r(i,-2,k)<0.0 )then
                uvwt(4,-2) = -1.0
                st = 0
            else
                st = -2
            endif

            if( r(i,nj+3,k)<0.0 )then
                uvwt(4,nj+3) = -1.0
                ed = nj+1
            else
                ed = nj+3
            endif

			do j=st,ed
				uvwt(1,j) = u(i,j,k)
				uvwt(2,j) = v(i,j,k)
				uvwt(3,j) = w(i,j,k)
				uvwt(4,j) = t(i,j,k)
			enddo

			call DUVWT_NODE_LINE_virtual(nmax,nj,4,uvwt,dtem)

			do j=1,nj
				do m=1,4
					duvwt(m+4,i,j,k)=dtem(m,j)
				enddo
			enddo

		enddo
	enddo

! K direction

	do j=1,nj
		do i=1,ni
            if( r(i,j,-2)<0.0 )then
                uvwt(4,-2) = -1.0
                st = 0
            else
                st = -2
            endif

            if( r(i,j,nk+3)<0.0 )then
                uvwt(4,nk+3) = -1.0
                ed = nk+1
            else
                ed = nk+3
            endif

			do k=st,ed
				uvwt(1,k) = u(i,j,k)
				uvwt(2,k) = v(i,j,k)
				uvwt(3,k) = w(i,j,k)
				uvwt(4,k) = t(i,j,k)
			enddo

			call DUVWT_NODE_LINE_virtual(nmax,nk,4,uvwt,dtem)

			do k=1,nk
				do m=1,4
					duvwt(m+8,i,j,k)=dtem(m,k)
				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_4th_virtual
!=============================================================================!
!=============================================================================!
subroutine UVWT_DER_4th_half_virtual    !____ֱ�������1�׵��� 2009.2.1
  use global_variables,only:u,v,w,r,t,ni,nj,nk,nmax
  use duvwt_all_field,only: duvwt_mid
  implicit none
	integer :: i,j,k,m,st,ed
	real    :: uvwt(4,-2:nmax+3),duvwt(4,0:nmax)
!	real    :: duvwtdkc(nv,0:ni,0:nj,0:nk)

! I direction

	do k=1,nk
		do j=1,nj
            if( r(-2,j,k)<0.0 )then
                uvwt(4,-2) = -1.0
                st = 0
            else
                st = -2
            endif

            if( r(ni+3,j,k)<0.0 )then
                uvwt(4,ni+3) = -1.0
                ed = ni+1
            else
                ed = ni+3
            endif

			do i=st,ed
				uvwt(1,i) = u(i,j,k)
				uvwt(2,i) = v(i,j,k)
				uvwt(3,i) = w(i,j,k)
				uvwt(4,i) = t(i,j,k)
			enddo

			call DUVWT_half_line_virtual(nmax,ni,4,uvwt,duvwt)     !*TGH. Ҫ������ϵ�ֵ

			do i=0,ni
				do m=1,4

					duvwt_mid(m,i,j,k)=duvwt(m,i)

				enddo
			enddo

		enddo
	enddo
!
! J direction

	do k=1,nk
		do i=1,ni
            if( r(i,-2,k)<0.0 )then
                uvwt(4,-2) = -1.0
                st = 0
            else
                st = -2
            endif

            if( r(i,nj+3,k)<0.0 )then
                uvwt(4,nj+3) = -1.0
                ed = nj+1
            else
                ed = nj+3
            endif

			do j=st,ed
				uvwt(1,j) = u(i,j,k)
				uvwt(2,j) = v(i,j,k)
				uvwt(3,j) = w(i,j,k)
				uvwt(4,j) = t(i,j,k)
			enddo

			call DUVWT_half_line_virtual(nmax,nj,4,uvwt,duvwt)     !*TGH. Ҫ������ϵ�ֵ

			do j=0,nj
				do m=1,4

					duvwt_mid(m+4,i,j,k)=duvwt(m,j)

				enddo
			enddo

		enddo
	enddo

! K direction

	do j=1,nj
		do i=1,ni
            if( r(i,j,-2)<0.0 )then
                uvwt(4,-2) = -1.0
                st = 0
            else
                st = -2
            endif

            if( r(i,j,nk+3)<0.0 )then
                uvwt(4,nk+3) = -1.0
                ed = nk+1
            else
                ed = nk+3
            endif

			do k=st,ed
				uvwt(1,k) = u(i,j,k)
				uvwt(2,k) = v(i,j,k)
				uvwt(3,k) = w(i,j,k)
				uvwt(4,k) = t(i,j,k)
			enddo

			call DUVWT_half_line_virtual(nmax,nk,4,uvwt,duvwt)     !*TGH. Ҫ������ϵ�ֵ

			do k=0,nk
				do m=1,4

					duvwt_mid(m+8,i,j,k)=duvwt(m,k)

				enddo
			enddo

		enddo
	enddo
!
  return
end subroutine UVWT_DER_4th_half_virtual

!=============================================================================!

!=============================================================================!

subroutine DUVWT_NODE_LINE(nmax,ni,n1,uvwt,duvwt)  !___DERINODE
!---------------------------------------------------------------!
!* ����һ�����Ͻ���1�׵������������
!---------------------------------------------------------------!
	use global_variables,only : nijk2nd
    use define_precision_mod
    implicit none

	real(prec) :: A1,B1,dd12,dd6
	real(prec) :: A2,B2,C2,D2,E2
!	real(prec) :: AC2,BC2,CC2,DC2,EC2

	integer :: nmax,i,m,ni,n1,n2
	real(prec) :: uvwt(n1,-1:nmax+1),duvwt(n1,nmax)

    A1=  8.0_prec ; B1= -1._prec ; dd12=12._prec ; dd6 =6._prec
    A2= -3.0_prec ; B2=-10._prec ; C2=  18._prec ; D2= -6._prec ; E2=  1._prec

	if( ni <= nijk2nd ) then   !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
		do i=2,ni-1
		do m=1,n1
			duvwt(m,i) = 0.5*(uvwt(m,i+1) - uvwt(m,i-1))
		enddo
		enddo

		do m=1,n1
			duvwt(m,1) = uvwt(m,2 ) - uvwt(m,1)
			duvwt(m,ni)= uvwt(m,ni) - uvwt(m,ni-1)
		enddo

	else
		do i=3,ni-2
		do m=1,n1
				duvwt(m,i) = (A1*(uvwt(m,i+1) - uvwt(m,i-1)) + &
										  B1*(uvwt(m,i+2) - uvwt(m,i-2)))/12.d0
		enddo
		enddo

		do m=1,n1
			duvwt(m,2   ) =  (A2*uvwt(m,1) + B2*uvwt(m,2) + C2*uvwt(m,3) + &
			                  D2*uvwt(m,4) + E2*uvwt(m,5))/dd12
			duvwt(m,ni-1) = -(A2*uvwt(m,ni  ) + B2*uvwt(m,ni-1) + C2*uvwt(m,ni-2) + &
			                  D2*uvwt(m,ni-3) + E2*uvwt(m,ni-4))/dd12


			duvwt(m,1   ) =  (-11.d0*uvwt(m,1) +18.d0*uvwt(m,2)   -9*uvwt(m,3)   +2*uvwt(m,4))/dd6
			duvwt(m,ni  ) = -(-11.d0*uvwt(m,ni)+18.d0*uvwt(m,ni-1)-9*uvwt(m,ni-2)+2*uvwt(m,ni-3))/dd6
		enddo

	endif
!
  return
end subroutine DUVWT_NODE_LINE

!=============================================================================!
!=============================================================================!
subroutine DUVWT_NODE_LINE_virtual(nmax,ni,n1,uvwt,duvwt)  !___DERINODE
!---------------------------------------------------------------!
!* ����һ�����Ͻ���1�׵��������ݵ�3���ܶ��¶��Ƿ�Ҫ�����
!---------------------------------------------------------------!
	use global_variables,only : nijk2nd
    use define_precision_mod
    implicit none

	real(prec) :: A1,B1,dd12,dd6
	real(prec) :: A2,B2,C2,D2,E2
!	real(prec) :: AC2,BC2,CC2,DC2,EC2

	integer :: nmax,i,m,ni,n1,n2,st,ed
	real(prec) :: uvwt(n1,-2:nmax+3),duvwt(n1,nmax)

    A1=  8.0_prec ; B1= -1._prec ; dd12=12._prec ; dd6 =6._prec
    A2= -3.0_prec ; B2=-10._prec ; C2=  18._prec ; D2= -6._prec ; E2=  1._prec

	if( ni <= nijk2nd ) then   !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
		do i=1,ni
		do m=1,n1
			duvwt(m,i) = 0.5*(uvwt(m,i+1) - uvwt(m,i-1))
		enddo
		enddo
	else
        st = 1
        ed = ni
        if( uvwt(4,-2)  < 0.0 )then
          st = 3
		  do m=1,n1
			duvwt(m,1   ) =  (-11.*uvwt(m,1) +18.*uvwt(m,2)   - 9.*uvwt(m,3) + 2.*uvwt(m,4))/dd6
			duvwt(m,2   ) =  (A2*uvwt(m,1) + B2*uvwt(m,2) + C2*uvwt(m,3) + &
			                  D2*uvwt(m,4) + E2*uvwt(m,5) )/dd12
		  enddo
        endif

        if( uvwt(4,ni+3)< 0.0 )then
          ed = ni-2
          do m=1,n1
			duvwt(m,ni  ) = -(-11.*uvwt(m,ni)+18.*uvwt(m,ni-1)-9.0*uvwt(m,ni-2)+2.*uvwt(m,ni-3))/dd6
			duvwt(m,ni-1) = -(A2*uvwt(m,ni  ) + B2*uvwt(m,ni-1) + C2*uvwt(m,ni-2) + &
			                  D2*uvwt(m,ni-3) + E2*uvwt(m,ni-4) )/dd12
		  enddo
        endif

		do i=st,ed
		do m=1,n1
				duvwt(m,i) = ( A1*(uvwt(m,i+1) - uvwt(m,i-1)) + &
						       B1*(uvwt(m,i+2) - uvwt(m,i-2)) )/12.d0
		enddo
		enddo

	endif
!
  return
end subroutine DUVWT_NODE_LINE_virtual

!=============================================================================!
!=============================================================================!

subroutine DUVWT_half_line_virtual(nmax,ni,n1,uvwt,duvwt)
!---------------------------------------------------------------!
!* ����һ�����ϰ����1�׵��������ݵ�3���ܶȷ����ж��Ƿ�Ҫ�����
!---------------------------------------------------------------!
	use global_variables,only : nijk2nd
    use define_precision_mod
    implicit none

	integer :: nmax,i,m,ni,n1,n2,st,ed
	real(prec) :: A1,B1,dd24
	real(prec) :: A2,B2,C2,D2,E2,a3,b3,c3,d3
	real(prec) :: uvwt(n1,-2:nmax+3),duvwt(n1,0:nmax)

	A1= 27._prec ; B1=-1._prec ; dd24 = 24._prec
    A2=-22._prec ; B2=17._prec ; C2   = 9._prec ; D2=-5._prec ; E2=1._prec !*TGH. OLD
    A3=-71._prec ; B3=141._prec; C3   =-93._prec; D3=23._prec 

	if( ni <= nijk2nd )then  !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��
		do m=1,n1
		do i=0,ni
				duvwt(m,i) = uvwt(m,i+1) - uvwt(m,i)
		enddo
		enddo
	else

        st = 0
        ed = ni
        if( uvwt(4,-2)  < 0.0 )then
          st = 2
		  do m=1,n1
			duvwt(m,0   ) =  (A3*uvwt(m,1)  + B3*uvwt(m,2)    + C3*uvwt(m,3)    + D3*uvwt(m,4) )/dd24 
			duvwt(m,1   ) =  (A2*uvwt(m,1)  + B2*uvwt(m,2)    + C2*uvwt(m,3)    + D2*uvwt(m,4)    + E2*uvwt(m,5)   )/dd24
		  enddo
        endif

        if( uvwt(4,ni+3)< 0.0 )then
          ed = ni-2
		do m=1,n1
			duvwt(m,ni  ) = -(A3*uvwt(m,ni) + B3*uvwt(m,ni-1) + C3*uvwt(m,ni-2) + D3*uvwt(m,ni-3) )/dd24
			duvwt(m,ni-1) = -(A2*uvwt(m,ni) + B2*uvwt(m,ni-1) + C2*uvwt(m,ni-2) + D2*uvwt(m,ni-3) + E2*uvwt(m,ni-4))/dd24
		enddo
        endif

		do i=st,ed
		do m=1,n1
				duvwt(m,i) = ( A1*(uvwt(m,i+1) - uvwt(m,i  )) + &
						       B1*(uvwt(m,i+2) - uvwt(m,i-1))  )/dd24
		enddo
		enddo

	endif
!
  return
end subroutine DUVWT_half_line_virtual

!=============================================================================!
!=============================================================================!

subroutine VALUE_line_half_virtual(n,nmax,ni,q,q_half)
!---------------------------------------------------------------!
!* ��ֵһ�����ϰ����ֵ�����ݵ�3���ܶ��ж��Ƿ�Ҫ�����
!---------------------------------------------------------------!
	use global_variables,only : nijk2nd
    use define_precision_mod
	implicit none

	real(prec) A1,B1,dd16
	real(prec) A2,B2,C2,D2,a3,b3,c3,d3
	integer :: nmax,n,ni,m,i,st,ed
	real    :: q(n,-2:nmax+3),q_half(n,0:nmax)

	a1=9._prec ; B1=-1._prec ; dd16= 16._prec
    A2=5._prec ; B2=15._prec ; C2=   -5._prec ; D2=1._prec
    A3=35._prec; B3=-35._prec; C3=   21._prec ; D3=-5._prec

	if(ni <= nijk2nd)then  !*TGH. ע�⣬�����Ƕ�ά����Ϊ2��

		do i=0,ni
		do m=1,n
				q_half(m,i) = 0.5*(q(m,i) + q(m,i+1))
		enddo
		enddo

	else
       st = 0
       ed = ni
       if( q(4,-2) < 0.0 )then
          st = 2
		  do m=1,n
			q_half(m,0   ) = (A3*q(m,1   )+B3*q(m,2   )+C3*q(m,3   )+D3*q(m,4   ))/dd16 
			q_half(m,1   ) = (A2*q(m,1   )+B2*q(m,2   )+C2*q(m,3   )+D2*q(m,4   ))/dd16
		  enddo
        endif
        if( q(4,ni+3) < 0.0 )then
          ed = ni-2
		  do m=1,n
			q_half(m,ni  ) = (A3*q(m,ni  )+B3*q(m,ni-1)+C3*q(m,ni-2)+D3*q(m,ni-3))/dd16
			q_half(m,ni-1) = (A2*q(m,ni  )+B2*q(m,ni-1)+C2*q(m,ni-2)+D2*q(m,ni-3))/dd16
		  enddo
        endif


		do i=st,ed
		do m=1,n
			q_half(m,i) = (A1*(q(m,i) + q(m,i+1)) + B1*(q(m,i+2) + q(m,i-1)))/dd16
		enddo
		enddo
	endif

	return
end subroutine VALUE_line_half_virtual

!=============================================================================!
!=============================================================================!

!=============================================================================!
subroutine dif_cic_new(nb,nr) !�ж�����ֵʱ���ǵ���������ʱΪ��Сֵ�����
!*TGH. �����Խ��޸�rhs����dq��������ϵ�ֵ��¼�Խ����϶�Ӧ��Ŀ����rhs
!*TGH. ��Ҫ�õ��Խ����Լ�����ϵ�rhs1�����������õ��Խ����ϵ�dq��
!*tgh. �����Խ���Ȼ�ü�ƽ��
!*TGH. ʹ��ǰ����Recast������������
!* ����ϵ������Ϊeps���ڡ�subroutine Aps_dq_new�� ��
!* ��ƣ�Ϳ����                             ----------------------------------!
!* ���ԣ�Ϳ����  2009.3                     ----------------------------------!

	use define_precision_mod
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

             is1 = i + s_lr3d(1)
             js1 = j + s_lr3d(2)  !*TGH. �Խ��棨Դ��������һ�������
             ks1 = k + s_lr3d(3)

             is2 = i - s_lr3d(1)
             js2 = j - s_lr3d(2)  !*TGH. �Խ��棨Դ����������һ�������
             ks2 = k - s_lr3d(3)

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
			 nvolt = mb_vol(nbt)%a3d(it,jt,kt)


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
!			    dqq(n1) = mb_dq(nb)%a4d(n1,i,j,k)    !*TGH. Դ���ϵĶԽ���
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
!			     mb_dq(nb)%a4d(n1,i,j,k) = (qsp(n1)+qsn(n1))    !*TGH. ??? ʹ���������RHS?��dq?
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


			 if(rm<1.e-20 .or. pm<1.e-20 .or. rm>1.e10 .or. pm>1.e10 )then
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
end subroutine dif_cic_new

!=============================================================================!
!=============================================================================!

subroutine Aps_dq_new(prim,nxyz,s_lr,dq,f,npn)
!*TGH. ��A~��A+������ʾ��Ȼ��ֱ�Ӷ�RHS����dq������ת��
!*TGH.  f = A^-.Dq ��npn=1ʱ�� ��f = A^+.Dq ��npn=-1ʱ��
!* ��ƣ�Ϳ����                             ----------------------------------!
!* ���ԣ�Ϳ����  2009.3                     ----------------------------------!
!* By TU Guohua

	use define_precision_mod
    use global_const,only:nl,ns,ms1,beta1,tref,gama,sml_sss
    implicit none
    integer :: npn,m,s_lr
    real(prec) :: hint !,gama
    real(prec) :: prim(nl),dq(nl),hs(ns),as(ns)
    real(prec) :: nx,ny,nz,nt,ct,cgm,cgm1,nxyz(4),f(nl)
    real(prec) :: l1,l4,l5,x1,x2
    real(prec) :: dh,dc,c2dc,ae,af
    real(prec) :: rm,um,vm,wm,pm,cm,c2,v2,tm,hm
	real(prec) :: eps,eps0,absl1,absl4,absl5,dl145

	eps0 = 0.0001_prec !����ֵΪ��Сֵʱ����

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)


	hint = gama/(gama-1.0)*pm/rm
    c2 = gama*pm/rm
    cm = sqrt(c2)
    v2 = um*um + vm*vm + wm*wm

    hm = hint + 0.5*v2

	  nx = nxyz(1)
	  ny = nxyz(2)
	  nz = nxyz(3)
	  nt = nxyz(4)

    ct = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
	eps = eps0*cgm

    l1 = ct
    l4 = ct + cm * cgm
    l5 = ct - cm * cgm

    l1 = l1*s_lr
    l4 = l4*s_lr
    l5 = l5*s_lr

!------------------
!    old
!    l1 = 0.5*( 1.0 + npn*sign(1.0,l1) )
!    l4 = 0.5*( 1.0 + npn*sign(1.0,l4) )
!    l5 = 0.5*( 1.0 + npn*sign(1.0,l5) )
!
!------------------
!    new
	 absl1 = abs(l1)
	 absl4 = abs(l4)
	 absl5 = abs(l5)

	 dl145 = (absl1+eps) * 2._prec
	 l1 = (absl1 + npn*l1 + eps) / dl145

	 dl145 = (absl4+eps) * 2._prec
	 l4 = (absl4 + npn*l4 + eps) / dl145

	 dl145 = (absl5+eps) * 2._prec
	 l5 = (absl5 + npn*l5 + eps) / dl145
!------------------


    x1 = ( 2.0*l1 - l4 - l5 )/( 2.0 * c2 )
    x2 = ( l4 - l5 )/( 2.0 * cm )

    cgm1 = 1.0/cgm
    nx = nx * cgm1
    ny = ny * cgm1
    nz = nz * cgm1
    ct = ( ct - nt ) * cgm1
    ae = gama - 1.0

       af = 0.5*ae*v2

    dc = ct * dq(1) - nx * dq(2) - ny * dq(3) - nz * dq(4)
    dh = af * dq(1) - ae * ( um * dq(2) + vm * dq(3) + wm * dq(4) - dq(5) )


    c2dc = c2 * dc

    f(1) = l1 * dq(1)             -    dh   * x1           -    dc   * x2
    f(2) = l1 * dq(2) + ( nx*c2dc - um*dh ) * x1 + ( nx*dh - um*dc ) * x2
    f(3) = l1 * dq(3) + ( ny*c2dc - vm*dh ) * x1 + ( ny*dh - vm*dc ) * x2
    f(4) = l1 * dq(4) + ( nz*c2dc - wm*dh ) * x1 + ( nz*dh - wm*dc ) * x2
    f(5) = l1 * dq(5) + ( ct*c2dc - hm*dh ) * x1 + ( ct*dh - hm*dc ) * x2

    return
end subroutine Aps_dq_new
!_____________________________________________________________________!

!=============================================================================!
!=============================================================================!


subroutine GCL_analysis_pre
!* �ó�����ʱʹ�ã��ȵ�����Ϻ󣬿ɼ���preprecoss��
	use define_precision_mod
	use global_variables,only : mb_dim,nl,nmax,nblocks
	use store_GCL
	implicit none
	integer :: m

	allocate( mb_gcl0(nblocks) )
	do m=1,nblocks
	    allocate( mb_gcl0(m)%a4d(4,mb_dim(m,1),mb_dim(m,2),mb_dim(m,3)) )
	enddo

end subroutine GCL_analysis_pre

!=====================================================================!
!=====================================================================!

subroutine GCL_analysis
	use define_precision_mod
	use store_GCL
	use global_variables,only : nblocks,nscheme
	implicit none
	integer :: nb,fileid
	character(len=120) :: filegcl,title_gcl
	call GCL_analysis_pre

	fileid=20
	filegcl='output/tecgcl.dat'
    title_gcl ='variables = "x" "y" "z" "Ix" "Iy" "Iz" "It"'

    open(fileid,file=trim(filegcl),status='unknown')


    do nb = 1,nblocks

        call recast_grid(nb)

		gcl0 => mb_gcl0(nb)%a4d

		if(nscheme ==4 )then
		    call gcl_wcnse5(nb)
		elseif(nscheme == 1 )then
			call gcl_m2nd(nb)
		else
			write(*,*)'û�ж�Ӧ��������ʽ�ļ����غ��ɷ���ģ��'
			write(*,*)'����nscheme = 4��1'
			STOP
		endif

	    call tecplot_gcl(nb,fileid,title_gcl)

	enddo

	close(fileid)

	call GCL_analysis_post

	stop 'Finished GCL Analisis'

end subroutine GCL_analysis

!=====================================================================!
!=====================================================================!
subroutine GCL_analysis_post
!* �ó�����ʱʹ�ã��ȵ�����Ϻ󣬿ɼ���preprecoss��
	use define_precision_mod
	use global_variables,only : nblocks,ni,nj,nk
	use store_GCL
	implicit none
	integer :: nb,i,j,k
	real(prec) :: ixm,iym,izm,itm,ixyzt,ixaver,iyaver,izaver,itaver,ixall,iyall,izall,itall,ixyztall
	real(prec) :: ixabs,iyabs,izabs,itabs
	integer :: im,jm,km,ima,jma,kma,nbm
	real(prec) :: points,pointsall


	open(2,file='output/gcl_result.dat')

	ixall = 0.0
	iyall = 0.0
	izall = 0.0
	itall = 0.0
	ixyztall = 0.0 !*���ֵ
	pointsall = 0.0

	do nb=1,nblocks

        call recast_grid(nb)

		gcl0 => mb_gcl0(nb)%a4d

		ixm = 0.0
		iym = 0.0
		izm = 0.0
		itm = 0.0
		ixyzt = 0.0

		ixaver= 0.0
		iyaver= 0.0
		izaver= 0.0
		itaver= 0.0

		do k=1,nk
		do j=1,nj
		do i=1,ni

			ixabs = abs(gcl0(1,i,j,k))
			iyabs = abs(gcl0(2,i,j,k))
			izabs = abs(gcl0(3,i,j,k))
			itabs = abs(gcl0(4,i,j,k))

			ixaver = ixaver + ixabs
			iyaver = iyaver + iyabs
			izaver = izaver + izabs
			itaver = itaver + itabs

			if(ixabs > ixm )ixm = ixabs
			if(iyabs > ixm )iym = iyabs
			if(izabs > ixm )ixm = izabs
			if(itabs > ixm )itm = itabs
			if(ixabs+iyabs+izabs > ixyzt )then
				ixyzt = ixabs+iyabs+izabs+itabs !* �����غ��ɾ���ֵ֮�͵����ֵ
				im = i
				jm = j
				km = k
				if(ixyzt > ixyztall)then
					ixyztall = ixyzt        !* ���п鼸���غ��ɾ���ֵ֮�͵����ֵ
				    nbm = nb
					ima = im
					jma = jm
					kma = km
				endif
			endif

		enddo
		enddo
		enddo

		ixall = ixall + ixaver
		iyall = iyall + iyaver
		izall = izall + izaver
		itall = itall + itaver

		points = real(ni*nj*nk)
		pointsall = pointsall + points
		ixaver = ixaver/points
		iyaver = iyaver/points
		izaver = izaver/points
		itaver = itaver/points

		write(2,9 )nb,im,jm,km
		write(2,10)nb,ixm,iym,izm,itm
		write(2,11)nb,ixaver,iyaver,izaver,itaver
		write(2,8 )ixyzt,ixaver+iyaver+izaver+itaver
		write(2,*)
		write(2,*)

9		format(1x,'Block No.=',i3,', Max_err_point at (i,j,k):',3i4)
10		format(1x,'Block No.=',i3,', Ix,Iy,Iz,It (max)  =',4e13.5)
11		format(1x,'Block No.=',i3,', Ix,Iy,Iz,It (aver) =',4e13.5)
8       format(1x,'Ix+Iy+Iz+It (nb, Max) =',e13.5', Ix+Iy+Yz+It (nb, aver) =',e13.5)

	enddo

	ixall = ixall/pointsall
	iyall = iyall/pointsall
	izall = izall/pointsall
	itall = itall/pointsall

	write(2,*)'========All Blocks======='
	write(2,12)nbm,ima,jma,kma
	write(2,13)ixall,iyall,izall,itall
	write(2,14)ixyztall,ixall+iyall+izall+itall

12	format(1x,'Max_err_position (nb,i,j,k) = ',4i4)
13	format(1x,'Ix,Iy,Iz,It (all,Aver) =',4e13.5)
14	format(1x,'Ix+Iy+Iz+It (all,Max)  =',e13.5,', Ix+Iy+Iz (all,aver)=',e13.5)
	close(2)

!********************************!
!* �ͷ��ڴ�
	do nb=1,nblocks
	    deallocate( mb_gcl0(nb)%a4d )
	enddo
	deallocate( mb_gcl0 )

end subroutine GCL_analysis_post

!=====================================================================!
!=====================================================================!

subroutine gcl_wcnse5(nb)
!* ��������WCNS-E-5����õ�GCL
	use define_precision_mod
	use store_GCL
	use global_variables,only : kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,gcl0 &
	                          , ni,nj,nk,nmax
	implicit none
	integer :: nb,m,i,j,k
	real(prec) :: nxyz(3,nmax),nxyz_half(3,0:nmax),dnxyz(3,nmax)

	do k=1,nk
	do j=1,nj
	do i=1,ni
	do m=1,4
		gcl0(m,i,j,k) = 0._prec
	enddo
	enddo
	enddo
	enddo

!* ����Ϊ I ����Լ����غ��ɵĹ���
	do k=1,nk
	do j=1,nj
		do i=1,ni
		    nxyz(1,i) = kcx(i,j,k)
		    nxyz(2,i) = kcy(i,j,k)
		    nxyz(3,i) = kcz(i,j,k)
		enddo

		call VALUE_HALF_NODE(3,nmax,ni,nxyz,nxyz_half)
		call FLUX_DXYZ(3,nmax,ni,nxyz_half,dnxyz)

		do i=1,ni
		do m=1,3
			gcl0(m,i,j,k) = gcl0(m,i,j,k) + dnxyz(m,i)
		enddo
		enddo
	enddo
	enddo

!* ���� I ����Լ����غ��ɵĹ���
!* ����Ϊ J ����Լ����غ��ɵĹ���
	do k=1,nk
	do i=1,ni
		do j=1,nj
		    nxyz(1,j) = etx(i,j,k)
		    nxyz(2,j) = ety(i,j,k)
		    nxyz(3,j) = etz(i,j,k)
		enddo

		call VALUE_HALF_NODE(3,nmax,nj,nxyz,nxyz_half)
		call FLUX_DXYZ(3,nmax,nj,nxyz_half,dnxyz)

		do j=1,nj
		do m=1,3
			gcl0(m,i,j,k) = gcl0(m,i,j,k) + dnxyz(m,j)
		enddo
		enddo
	enddo
	enddo
!* ���� J ����Լ����غ��ɵĹ���
!* ����Ϊ K ����Լ����غ��ɵĹ���
	do j=1,nj
	do i=1,ni
		do k=1,nk
		    nxyz(1,k) = ctx(i,j,k)
		    nxyz(2,k) = cty(i,j,k)
		    nxyz(3,k) = ctz(i,j,k)
		enddo

		call VALUE_HALF_NODE(3,nmax,nk,nxyz,nxyz_half)
		call FLUX_DXYZ(3,nmax,nk,nxyz_half,dnxyz)

		do k=1,nk
		do m=1,3
			gcl0(m,i,j,k) = gcl0(m,i,j,k) + dnxyz(m,k)
		enddo
		enddo
	enddo
	enddo
!* ���� K ����Լ����غ��ɵĹ���

end subroutine gcl_wcnse5

!=====================================================================!
!=====================================================================!

subroutine gcl_M2nd(nb)
!* �������ö��׸�ʽ����õ�GCL
	use define_precision_mod
	use store_GCL
	use global_variables,only : kcx,kcy,kcz,etx,ety,etz,ctx,cty,ctz,gcl0 &
	                          , ni,nj,nk,nmax
	implicit none
	integer :: nb,m,i,j,k
	real(prec) :: nxyz(3,nmax),nxyz_half(3,0:nmax),dnxyz(3,nmax)

	do k=1,nk
	do j=1,nj
	do i=1,ni
	do m=1,4
		gcl0(m,i,j,k) = 0._prec
	enddo
	enddo
	enddo
	enddo

!* ����Ϊ I ����Լ����غ��ɵĹ���
	do k=1,nk
	do j=1,nj
		do i=1,ni
		    nxyz(1,i) = kcx(i,j,k)
		    nxyz(2,i) = kcy(i,j,k)
		    nxyz(3,i) = kcz(i,j,k)
		enddo
		call VALUE_HALF_2nd(3,nmax,ni,nxyz,nxyz_half)
		call FLUX_DXYZ_2nd(3,nmax,ni,nxyz_half,dnxyz)

		do i=1,ni
		do m=1,3
			gcl0(m,i,j,k) = gcl0(m,i,j,k) + dnxyz(m,i)
		enddo
		enddo
	enddo
	enddo

!* ���� I ����Լ����غ��ɵĹ���
!* ����Ϊ J ����Լ����غ��ɵĹ���
	do k=1,nk
	do i=1,ni
		do j=1,nj
		    nxyz(1,j) = etx(i,j,k)
		    nxyz(2,j) = ety(i,j,k)
		    nxyz(3,j) = etz(i,j,k)
		enddo

		call VALUE_HALF_2nd(3,nmax,nj,nxyz,nxyz_half)
		call FLUX_DXYZ_2nd(3,nmax,nj,nxyz_half,dnxyz)

		do j=1,nj
		do m=1,3
			gcl0(m,i,j,k) = gcl0(m,i,j,k) + dnxyz(m,j)
		enddo
		enddo
	enddo
	enddo
!* ���� J ����Լ����غ��ɵĹ���
!* ����Ϊ K ����Լ����غ��ɵĹ���
	do j=1,nj
	do i=1,ni
		do k=1,nk
		    nxyz(1,k) = ctx(i,j,k)
		    nxyz(2,k) = cty(i,j,k)
		    nxyz(3,k) = ctz(i,j,k)
		enddo

		call VALUE_HALF_2nd(3,nmax,nk,nxyz,nxyz_half)
		call FLUX_DXYZ_2nd(3,nmax,nk,nxyz_half,dnxyz)

		do k=1,nk
		do m=1,3
			gcl0(m,i,j,k) = gcl0(m,i,j,k) + dnxyz(m,k)
		enddo
		enddo
	enddo
	enddo
!* ���� K ����Լ����غ��ɵĹ���

end subroutine gcl_M2nd

!=====================================================================!
!=====================================================================!


subroutine tecplot_gcl(nb,fileid,title_gcl)
	use store_GCL
    use global_variables,only : mb_dim,ni,nj,nk,nblocks,mb_bc,x,y,z
    implicit none
	integer :: nb,fileid
	character(len=120) :: title_gcl
	integer :: nrmax,bctype,i,j,k,m,n,nr
	integer :: nim,njm,nkm
	integer :: s_st(3),s_ed(3)
    logical :: ldraw

	nim = (1+ni)/2    !* I�������λ��
	njm = (1+nj)/2    !* J�������λ��
	nkm = (1+nk)/2    !* K�������λ��

	if( (nim/=1) .and. (nim /=ni) )then
	   i=nim
         write(fileid,*)trim(title_gcl)

         write(fileid,*)'zone i=',1,', j=',nj,', k=',nk,', f=point'

         write(fileid,*)((( x(i,j,k),y(i,j,k),z(i,j,k),gcl0(1,i,j,k),gcl0(2,i,j,k),gcl0(3,i,j,k),gcl0(4,i,j,k) &
			                 & ,i=nim,nim),j=1,nj),k=1,nk)

	endif

	if( (njm/=1) .and. (njm /=nj) )then
		j=njm
          write(fileid,*)trim(title_gcl)

         write(fileid,*)'zone i=',ni,', j=',1,', k=',nk,', f=point'

         write(fileid,*)((( x(i,j,k),y(i,j,k),z(i,j,k),gcl0(1,i,j,k),gcl0(2,i,j,k),gcl0(3,i,j,k),gcl0(4,i,j,k) &
			                 & ,i=1,ni),j=njm,njm),k=1,nk)
	endif

	if( (nkm/=1) .and. (nkm /=nk) )then

	   k=nkm
         write(fileid,*)trim(title_gcl)

         write(fileid,*)'zone i=',ni,', j=',nj,', k=',1,', f=point'

         write(fileid,*)((( x(i,j,k),y(i,j,k),z(i,j,k),gcl0(1,i,j,k),gcl0(2,i,j,k),gcl0(3,i,j,k),gcl0(4,i,j,k) &
			                 & ,i=1,ni),j=1,nj),k=nkm,nkm)
	endif


    nrmax = mb_bc(nb)%nregions
    do nr=1,nrmax
        bctype = mb_bc(nb)%bc(nr)%bctype
		ldraw  = ( bctype    > 0  )

		if ( ldraw ) then
             do m=1,3
                s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
                s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
                if ( s_st(m) > s_ed(m) ) then
					n=s_st(m)
					s_st(m) = s_ed(m)
					s_ed(m) = n
                endif

             enddo

             write(fileid,*)trim(title_gcl)

             write(fileid,*)'zone i=',s_ed(1)-s_st(1)+1,', j=',s_ed(2)-s_st(2)+1,', k=',s_ed(3)-s_st(3)+1,', f=point'

             write(fileid,*)((( x(i,j,k),y(i,j,k),z(i,j,k),gcl0(1,i,j,k),gcl0(2,i,j,k),gcl0(3,i,j,k),gcl0(4,i,j,k) &
			                 & ,i=s_st(1),s_ed(1)),j=s_st(2),s_ed(2)),k=s_st(3),s_ed(3))
		endif

	enddo

end subroutine tecplot_gcl

!=====================================================================!
!=====================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timestep_match
!----------------------------------------------------------------
!
!  ƥ��Խӱ߽������ʱ�䲽�����Ա�֤�Խӱ߽總����ʱ�䲽��Э��
!  ���������ʱ�䲽����׼
!
!---------------------------------------------------------------

    use define_precision_mod
    use global_variables, only : nblocks,mb_bc,mb_dtdt,mb_vol,mb_dim,ntmst,dtdt,vol,timedt_rate
    implicit none
    integer :: nb,nr,nrmax,bctype
	integer :: ni,nj,nk
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,it,jt,kt,i,j,k,nt,m
    integer :: is1,js1,ks1,it1,jt1,kt1
    integer :: s_st(3),s_ed(3)
	integer :: s_lr3d(3),t_lr3d(3)
	integer :: nbmin,imin,jmin,kmin
	real(prec) :: dtmin


	if ( (ntmst /= 0) .or. (timedt_rate < 1.) ) then ! �ǵ���ʱ�䲽�������иöγ������
		return
	endif

	dtmin=1.0e6

    do nb = 1,nblocks
      ni = mb_dim(nb,1)
      nj = mb_dim(nb,2)
      nk = mb_dim(nb,3)
	  dtdt => mb_dtdt(nb)%a3d
	  vol  => mb_vol(nb)%a3d
		
	  do k=1,nk
	  do j=1,nj
	  do i=1,ni
		if(dtmin > dtdt(i,j,k)*vol(i,j,k))then

		  dtmin = dtdt(i,j,k)*vol(i,j,k) 

		  nbmin=nb
		  imin=i
		  jmin=j
		  kmin=k
		endif
	  enddo
	  enddo
	  enddo

	enddo
	
!	write(*,*)'dtmin',dtmin/mb_vol(nbmin)%a3d(imin,jmin,kmin)


    do nb = 1,nblocks
      ni = mb_dim(nb,1)
      nj = mb_dim(nb,2)
      nk = mb_dim(nb,3)
	  dtdt => mb_dtdt(nb)%a3d
	  vol  => mb_vol(nb)%a3d
		
	  do k=1,nk
	  do j=1,nj
	  do i=1,ni
		if(dtdt(i,j,k)*vol(i,j,k) > timedt_rate*dtmin)then
		  dtdt(i,j,k) = timedt_rate*dtmin/vol(i,j,k)
		endif

	  enddo
	  enddo
	  enddo

	enddo
    do nb = 1,nblocks
      ni = mb_dim(nb,1)
      nj = mb_dim(nb,2)
      nk = mb_dim(nb,3)

      nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
    !���崦���߽�����
      do nr = 1,nrmax
         bctype = mb_bc(nb)%bc(nr)%bctype    !* �Խ����� *!
         if( bctype < 0 ) then               !* �Խӱ߽�
			
           do m=1,3
             s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
             s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
			 s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)
			 t_lr3d(m) = mb_bc(nb)%bc(nr)%t_lr3d(m)
           enddo

		   nbs   = mb_bc(nb)%bc(nr)%nbs             !��� 
		   nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
!		   s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!		   t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�

           do i = s_st(1),s_ed(1)
           do j = s_st(2),s_ed(2)
           do k = s_st(3),s_ed(3)
			  is = i
			  js = j
			  ks = k
              it = mb_bc(nb)%bc(nr)%image(i,j,k )
              jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )
              kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )
			 
			  is1 = i - s_lr3d(1)
			  js1 = j - s_lr3d(2)
			  ks1 = k - s_lr3d(3)

			  it1 = i - t_lr3d(1)
			  jt1 = j - t_lr3d(2)
			  kt1 = k - t_lr3d(3)

			  mb_dtdt(nbs)%a3d(is,js,ks) = min( mb_dtdt(nbs)%a3d(is,js,ks), &
			                                    mb_dtdt(nbt)%a3d(it,jt,kt), &
			                                    mb_dtdt(nbt)%a3d(is1,js1,ks1), &
			                                    mb_dtdt(nbt)%a3d(it1,jt1,kt1) )
			  mb_dtdt(nbs)%a3d(is1,js1,ks1) = mb_dtdt(nbs)%a3d(is,js,ks)
           enddo
           enddo
           enddo

         endif
      enddo

	enddo

		
    return
end subroutine timestep_match
!_____________________________________________________________________!
subroutine read_bc_periodic(fileid,nb,nr,s_nd,s_st,s_ed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! �������ڱ߽�Ķ�Ӧ���ϵ
!
!  �߽��ļ���ʽ
!     +iss      +ise      +jss      -jse      -kss       -kse        8   !Դ��
!     +ies      +ite      +jts      -jte      -kts       -kte       nbt   !���ڱ߽��Ŀ���
!	����ÿ�����򶼴ӵ����궼start��eend
!   �����š�������ʾ������򣬶��Ǳ�ʾԴ����Ŀ���������߶�Ӧ��ϵ
!    �����߶�Ӧ��ϵ�� (+,+) <----> (+,+); (+,-) <----> (+,-); (-,-) <----> (-,-)
!    !����(+,+) ����ʾ���ڱ߽�ķ���
!
!
! By, TU, Guohua   20090926
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use global_variables
    implicit none
    integer :: fileid
    integer :: ndif,nr,nrmax,ntemp,js1,js2,ks1,ks2,ls1,ls2
    integer :: m,n,i,j,k,imax,jmax,kmax,nbt,nb
    integer :: co,op(2),t_coor(2)
    integer :: s_nd,s_lr,s_fix,idelt,jdelt,kdelt
    integer :: t_nd,t_lr,t_fix,s_t_dirction(3,3)
    integer :: s_index,t_index,s_sign(3),t_sign(3),st_sign(3)
    integer,dimension(1:3) :: s_st,s_ed,t_st,t_ed
    
    read(1,*)t_st(1),t_ed(1),t_st(2),t_ed(2),t_st(3),t_ed(3),nbt
    do m=1,3
       mb_bc(nb)%bc(nr)%t_lr3d(m) = 0
       if( t_st(m) > 0 .and. t_ed(m) >0 ) then
           t_nd = m   !���ڱ߽�ķ���
           if( min (t_st(m), t_ed(m)) == (abs(t_st(m)-t_ed(m))+1) )then
               t_lr = -1
               mb_bc(nb)%bc(nr)%t_lr3d(m) = -1
           else
               t_lr =  1
               mb_bc(nb)%bc(nr)%t_lr3d(m) = 1
           endif

		   if(T_lr == -1)THEN
                T_fix = min( T_st(m), T_ed(m) )
		   ELSE
                T_fix = MAX( T_st(m), T_ed(m) )
		   ENDIF
       endif
    enddo

    mb_bc(nb)%bc(nr)%nbt   = nbt   !�Խӿ�
    mb_bc(nb)%bc(nr)%t_lr  = t_lr  !�Խ�λ�ã�-1---����߶Խӣ�1---���ұ߶Խ�
    mb_bc(nb)%bc(nr)%t_nd  = t_nd  !�Խӷ���1��2��3----�ֱ����I��J��K����Խ�
    mb_bc(nb)%bc(nr)%t_fix = t_fix !�Խ����λ�ã��������J����Խӣ��Խ���ΪJ=t_fix����(�ο���)

    !ȷ���Խӷ���
    do m=1,3
       do n=1,3
          s_t_dirction(m,n) = 0
       enddo
    enddo

	s_t_dirction(t_nd,s_nd) = 1

    do m=1,3
       if ( s_st(m)*s_ed(m) < 0  )then
          do n=1,3
             if ( t_st(n)*t_ed(n) < 0 ) then
                s_t_dirction(n        , m        ) = 1
                s_t_dirction(6-n-t_nd , 6-m-s_nd ) = 1
                goto 10
             endif
          enddo
       endif
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
end subroutine read_bc_periodic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!_____________________________________________________________________!
subroutine boundary_8_tgh(nb,nr,bctype) !���ڱ߽�
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ʹ�÷�Χ�����޲�ַ�, û�л�ѧ��Ӧ
!  ����ע�⣺����ģ�ͱ߽绹������
!
!  �߽��ļ���ʽ
!     +iss      +ise      -jss      +jse      -kss       -kse        8   !Դ��
!     -ies      +ite      -jts      -jte      +kts       +kte       nbt    !���ڱ߽��Ŀ���
!	����ÿ�����򶼴ӵ����궼start��eend
!   �����š�������ʾ������򣬶��Ǳ�ʾԴ����Ŀ���������߶�Ӧ��ϵ
!    �����߶�Ӧ��ϵ�� (+,+) <----> (+,+); (+,-) <----> (+,-); (-,-) <----> (-,-)
!    !����(+,+) ����ʾ���ڱ߽�ķ���
!
!
!
! By, TU, Guohua   20090926
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use define_precision_mod
    use global_variables, only : mb_bc,nl,method,mb_r,mb_u,mb_v,mb_w,mb_p,mb_t,phytime , &
	                           & R,U,V,W,P,T,X,Y,NI,NJ
    implicit none
    integer :: nwhole,nd,half_dim,whole_dim
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: nb,nr,bctype
    integer :: is,js,ks,it,jt,kt,i,j,k,m,n
	integer :: ist,jst,kst,its,jts,kts
    integer :: s_st(3),s_ed(3),s_lr3d(3),t_lr3d(3)
	logical :: periodic_1point
	REAL(PREC) :: EKC,RB,TAOB,TE,PI
	real(prec) :: u00,v00,time, xcore,ycore,XP,YP

!__________________________________________________________!
!�߽����ȷֵ
  if(.false.)THEN
	return
  ENDIF
!����.�߽����ȷֵ
!__________________________________________________________!



	if( method /= 1)then
		write(*,*)'�������޲�ַ�����subroutine boundary_8_tgh�����ܴ����������ڱ߽�'
		STOP
	endif

    t_lr = 0
    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)
       s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)
       t_lr3d(m) = mb_bc(nb)%bc(nr)%t_lr3d(m)
	enddo

	if(abs(s_lr3d(1))+abs(s_lr3d(2))+abs(s_lr3d(3)) >1)then
		write(*,*)'���Գ���'
		write(*,*)'�������ڱ߽����'
		pause 'subroutine boundary_8_tgh'
	endif

    s_nd  = mb_bc(nb)%bc(nr)%s_nd     !�߽��淽��:1,2,3��Ӧ��i,j,k
    s_lr  = mb_bc(nb)%bc(nr)%s_lr     !���ұ߽�-1,1��Ӧ�����ұ߽�
    s_fix = mb_bc(nb)%bc(nr)%s_fix    !�̶�����(fixed_coor)
    nbs   = mb_bc(nb)%bc(nr)%nbs      !��� 
    nbt   = mb_bc(nb)%bc(nr)%nbt      !���
	if( (s_st(s_nd)-s_ed(s_nd) == 0) )then
		 periodic_1point = .true.
	else
		 periodic_1point = .false.
	endif
	
	if(periodic_1point)then  !ֻ��1������Ϊ���ڱ߽磬ȡ���ڶ�Ӧ���ƽ��ֵ
       do is = s_st(1),s_ed(1)
       do js = s_st(2),s_ed(2)
       do ks = s_st(3),s_ed(3)
		 it = mb_bc(nb)%bc(nr)%image(is,js,ks)
		 jt = mb_bc(nb)%bc(nr)%jmage(is,js,ks)
		 kt = mb_bc(nb)%bc(nr)%kmage(is,js,ks)

		 ist = is + s_lr3d(1)
		 jst = js + s_lr3d(2)
		 kst = ks + s_lr3d(3)

		 its = it - s_lr3d(1)
		 jts = jt - s_lr3d(2)
		 kts = kt - s_lr3d(3)


		 mb_r(nbt)%a3d(it,jt,kt) =0.5* (mb_r(nbs)%a3d(is,js,ks) + mb_r(nbt)%a3d(it,jt,kt))
		 mb_u(nbt)%a3d(it,jt,kt) =0.5* (mb_u(nbs)%a3d(is,js,ks) + mb_u(nbt)%a3d(it,jt,kt))
		 mb_v(nbt)%a3d(it,jt,kt) =0.5* (mb_v(nbs)%a3d(is,js,ks) + mb_v(nbt)%a3d(it,jt,kt))
		 mb_w(nbt)%a3d(it,jt,kt) =0.5* (mb_w(nbs)%a3d(is,js,ks) + mb_w(nbt)%a3d(it,jt,kt))
		 mb_p(nbt)%a3d(it,jt,kt) =0.5* (mb_p(nbs)%a3d(is,js,ks) + mb_p(nbt)%a3d(it,jt,kt))
         mb_r(nbs)%a3d(is,js,ks) =  mb_r(nbt)%a3d(it,jt,kt)
         mb_u(nbs)%a3d(is,js,ks) =  mb_u(nbt)%a3d(it,jt,kt) 
         mb_v(nbs)%a3d(is,js,ks) =  mb_v(nbt)%a3d(it,jt,kt) 
         mb_w(nbs)%a3d(is,js,ks) =  mb_w(nbt)%a3d(it,jt,kt) 
         mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbt)%a3d(it,jt,kt)

         mb_r(nbs)%a3d(ist,jst,kst) =  mb_r(nbt)%a3d(its,jts,kts)  !������ϵ�ֵ
         mb_u(nbs)%a3d(ist,jst,kst) =  mb_u(nbt)%a3d(its,jts,kts) 
         mb_v(nbs)%a3d(ist,jst,kst) =  mb_v(nbt)%a3d(its,jts,kts) 
         mb_w(nbs)%a3d(ist,jst,kst) =  mb_w(nbt)%a3d(its,jts,kts) 
         mb_p(nbs)%a3d(ist,jst,kst) =  mb_p(nbt)%a3d(its,jts,kts)
		 
	   enddo
	   enddo
	   enddo

	else !����1������Ϊ���ڱ߽�

       do i = s_st(1), s_ed(1)
       do j = s_st(2), s_ed(2)
       do k = s_st(3), s_ed(3)
		 is = i + s_lr3d(1)
		 js = j + s_lr3d(2)
		 ks = k + s_lr3d(3)

		 it = mb_bc(nb)%bc(nr)%image(i,j,k) - t_lr3d(1)
		 jt = mb_bc(nb)%bc(nr)%jmage(i,j,k) - t_lr3d(2)
		 kt = mb_bc(nb)%bc(nr)%kmage(i,j,k) - t_lr3d(3)

         mb_r(nbs)%a3d(is,js,ks) =  mb_r(nbt)%a3d(it,jt,kt)
         mb_u(nbs)%a3d(is,js,ks) =  mb_u(nbt)%a3d(it,jt,kt) 
         mb_v(nbs)%a3d(is,js,ks) =  mb_v(nbt)%a3d(it,jt,kt) 
         mb_w(nbs)%a3d(is,js,ks) =  mb_w(nbt)%a3d(it,jt,kt) 
         mb_p(nbs)%a3d(is,js,ks) =  mb_p(nbt)%a3d(it,jt,kt)


!				call BC_peroidic_turbulence(nbs,is,js,ks,it,jt,kt,i,j,k) ! turbulence model


       enddo
       enddo
       enddo

	endif

    return    
end subroutine boundary_8_tgh


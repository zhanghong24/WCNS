!_____________________________________________________________________!
subroutine read_bc_connect(fileid,nb,nr,s_nd,s_st,s_ed)
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
end subroutine read_bc_connect
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine boundary_n1(nb,nr,bctype) !�Խӱ߽�����
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
            mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)
            mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
            mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
            mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
            mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)
										
			    	call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC

       enddo
       enddo
       enddo

       call dif_average(nb,nr,bctype,s_st,s_ed) !* TGH. �߽��ϵ�ֵ���������ƽ����ͬʱ���д�������ģʽ�߽� *!

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
            mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)
            mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
            mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
            mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
            mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)

   			    call BC_connect_turbulence(nbs,nbt,is,js,ks,it,jt,kt,i,j,k) ! turbulence model BC
                call BC_face_dif_turbulence(nb,nb,is,js,ks,ist,jst,kst,i,j,k)

       enddo
       enddo
       enddo

	endif
     
    return
end subroutine boundary_n1
!_____________________________________________________________________!
subroutine boundary_n1_vir3(nb,nr,bctype) !�Խӱ߽�����
!*TGH. ���ʱ���Ȱ�ԭʼ���������ȡ��Ӧ�ص�����ֵ����������mb_qkeҲȡ��Ӧ�ص����ϵ�ֵ��
!*TGH. Ȼ��call dif_average ��������call BC_face_dif_turbulence�����߽��ϵ�ֵ���ڱ߽������ֵ�ļ�ƽ��
!*TGH. �Ѿ��޸ĵ��ʺ������߽����������ǻ����ʺϻ�ѧ��Ӧ�����

    use global_variables
    implicit none
    integer :: nb,nr,bctype,m,n
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd
!	integer :: t_fix,t_lr
    integer :: is,js,ks,i,j,k,nt
    integer :: it,jt,kt,it0,jt0,kt0,ist,jst,kst
    integer :: s_st(3),s_ed(3),s_lr3d(3),t_lr3d(3)
    integer :: nv_s,nv_t,nfile_par_s,nfile_par_t,num_var_s,num_var_t,nchem_s,nchem_t
    real    :: prim_s1(nl),q_1(nl),zf1,ttt1
	logical :: logi_turb
	
	if(nvis==1 .and. nlamtur>=0)then
		logi_turb = .true.
	else
		logi_turb = .false.
	endif

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
	   s_lr3d(m) = mb_bc(nb)%bc(nr)%s_lr3d(m)
	   t_lr3d(m) = mb_bc(nb)%bc(nr)%t_lr3d(m)
    enddo
    s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
    s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!    s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
    nbs   = nb !mb_bc(nb)%bc(nr)%nbs             !��� 

    nbt   = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
    t_nd  = mb_bc(nb)%bc(nr)%t_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k 
!    t_lr  = mb_bc(nb)%bc(nr)%t_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
!    t_fix = mb_bc(nb)%bc(nr)%t_fix           !�̶�����(fixed_coor)


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

	   if(logi_turb)then
	      mb_vist(nbs)%a3d(is,js,ks) = mb_vist(nbt)%a3d(it,jt,kt)

		  mb_vist(nbs)%a3d(i,j,k)    =(mb_vist(nbs)%a3d(is,js,ks) + &
		                               mb_vist(nbs)%a3d(ist,jst,kst) ) * 0.5_prec
!		  mb_vist(nbs)%a3d(i,j,k) = mb_vist(nbs)%a3d(is,js,ks) + mb_vist(nbs)%a3d(ist,jst,kst) &
!		         - 0.5*( mb_vist(nbs)%a3d(ist-s_lr3d(1),jst-s_lr3d(2),kst-s_lr3d(3)) + &
!						 mb_vist(nbt)%a3d(it -t_lr3d(1),jt -t_lr3d(2),kt -t_lr3d(3)) )  ! 2�׾���
	     do m=1,nlamtur
	        mb_qke(nbs)%a4d(is,js,ks,m)  =  mb_qke(nbt)%a4d(it,jt,kt,m)
		    mb_qke(nbs)%a4d(i,j,k,m)     = (mb_qke(nbs)%a4d(is,js,ks,m) + &
		                                    mb_qke(nbs)%a4d(ist,jst,kst,m) ) * 0.5_prec

		 enddo
	  endif

       mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)  !��һ���ص��㸳ֵ
       mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
       mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
       mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
       mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)


       is = is + s_lr3d(1)  !��2���ص����λ��
       js = js + s_lr3d(2)
       ks = ks + s_lr3d(3)

       it = it - t_lr3d(1)
       jt = jt - t_lr3d(2)
       kt = kt - t_lr3d(3)

       mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)   !��2���ص��㸳ֵ
       mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
       mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
       mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
       mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)

       is = is + s_lr3d(1)  !��3���ص����λ��
       js = js + s_lr3d(2)
       ks = ks + s_lr3d(3)

       it = it - t_lr3d(1)
       jt = jt - t_lr3d(2)
       kt = kt - t_lr3d(3)

       mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt) * dble(1 - 2*cic1)  !��3���ص��㸳ֵ
       mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
       mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)
       mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
       mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)

    enddo
    enddo
    enddo
     
    return
end subroutine boundary_n1_vir3
!_____________________________________________________________________!
subroutine recast_grid(nb)
    use global_variables,only : ni,nj,nk,x,y,z,vol,  &
                                kcx   ,kcy   ,kcz   ,kct   ,&
                                etx   ,ety   ,etz   ,ett   ,&
                                ctx   ,cty   ,ctz   ,ctt   ,&
                                mb_x  ,mb_y  ,mb_z  ,       &
                                mb_kcx,mb_kcy,mb_kcz,mb_kct,&
                                mb_etx,mb_ety,mb_etz,mb_ett,&
                                mb_ctx,mb_cty,mb_ctz,mb_ctt,&
                                mb_vol,mb_dim
    use global_variables,only : nbself
    implicit none
    integer :: nb
    ni = mb_dim(nb,1)
    nj = mb_dim(nb,2)
    nk = mb_dim(nb,3)

    x => mb_x(nb)%a3d
    y => mb_y(nb)%a3d
    z => mb_z(nb)%a3d

    vol => mb_vol(nb)%a3d

    kcx => mb_kcx(nb)%a3d
    kcy => mb_kcy(nb)%a3d
    kcz => mb_kcz(nb)%a3d
    kct => mb_kct(nb)%a3d

    etx => mb_etx(nb)%a3d
    ety => mb_ety(nb)%a3d
    etz => mb_etz(nb)%a3d
    ett => mb_ett(nb)%a3d

    ctx => mb_ctx(nb)%a3d
    cty => mb_cty(nb)%a3d
    ctz => mb_ctz(nb)%a3d
    ctt => mb_ctt(nb)%a3d
    
    nbself = nb      !!lhy
    
    return
end subroutine recast_grid
!_____________________________________________________________________!
subroutine recast_field(nb)
    use global_variables
    integer :: nb
    ni = mb_dim(nb,1)
    nj = mb_dim(nb,2)
    nk = mb_dim(nb,3)

    r => mb_r(nb)%a3d
    u => mb_u(nb)%a3d
    v => mb_v(nb)%a3d
    w => mb_w(nb)%a3d
    p => mb_p(nb)%a3d
    dtdt => mb_dtdt(nb)%a3d

    visl => mb_visl(nb)%a3d
    vist => mb_vist(nb)%a3d

    q  => mb_q(nb)%a4d
    dq => mb_dq(nb)%a4d

   t => mb_t(nb)%a3d
   c => mb_c(nb)%a3d

   sra => mb_sra(nb)%a3d
   srb => mb_srb(nb)%a3d
   src => mb_src(nb)%a3d
   srva => mb_srva(nb)%a3d
   srvb => mb_srvb(nb)%a3d
   srvc => mb_srvc(nb)%a3d

    if ( nlamtur /= 0 ) then   !��������
       qke => mb_qke(nb)%a4d
	     ds_turbulence => mb_distance(nb)%a3D    !  Normal Distence from Solid Surface
    endif                                      !  Modified by liangzai
    
    nbself = nb      !!lhy

    return
end subroutine recast_field
!_____________________________________________________________________!
subroutine update(nb)
    use global_variables,only:nlamtur,nwallfun
    implicit none
    integer :: nb

    !!call update_ns(nb)
    call update_limited(nb)

    return
end subroutine update
!_____________________________________________________________________!
!_____________________________________________________________________!
subroutine update_ns(nb)
    use global_const,only : pmin_limit,pmax_limit,rmax_limit,rmin_limit
    use global_variables, &
    only : method,ni,nj,nk,nm,gama,q,dq,r,u,v,w,p,roo,uoo,voo,woo,poo,moo
    implicit none
    integer :: nb,m,m1,i,j,k,n_count,i0,j0,k0
    integer :: ii,jj,kk,np
	integer :: iend,jend,kend !* Added by TU Guohua
    real :: gamam1,p_min,r_min,p_max,r_max
    real :: rm,um,vm,wm,pm,em,rm1,rnp
    
    gamam1 = gama - 1.0
    
    p_min = pmin_limit
    p_max = pmax_limit
    
    r_min = rmin_limit
    r_max = rmax_limit

  	n_count = 0

!    m1 = 1 + method
!*TGH. ע����ʼ����
   m1=1
   iend=ni
   jend=nj
   kend=nk
    do k=m1,kend 
      do j=m1,jend    
         do i=m1,iend  
             do m=1,nm
               q(m,i,j,k) = q(m,i,j,k) + dq(m,i,j,k)
             enddo
          enddo
       enddo
    enddo

    do k=m1,kend 
      do j=m1,jend    
         do i=m1,iend  
             rm = q(1,i,j,k)
             rm1 = 1.0/rm
             um = q(2,i,j,k) * rm1
             vm = q(3,i,j,k) * rm1
             wm = q(4,i,j,k) * rm1
             em = q(5,i,j,k)
             pm = gamam1 * ( em - 0.5*rm*( um*um + vm*vm + wm*wm ) )
             if ( pm <= p_min .or. rm <= r_min  .or. &
		  	    & pm >= p_max .or. rm >= r_max   ) then  !________

                if ( pm <= 0.0 .or. rm <= 0.0 ) then
                   !if(n_count < 3 ) write(*,901)'p<0 at:',nb,i,j,k,pm,rm
901                                     format(1x,a7,i3,3i4,2e14.5)
				   
				   n_count = n_count + 1

                   rm = 0.0
                   um = 0.0
                   vm = 0.0
                   wm = 0.0
                   pm = 0.0
                   np = 0
				       do ii=-1,1
                       do jj =-1,1
                       do kk =-1,1
                            if ( ( abs(ii) +abs(jj) +abs(kk) ) == 1 ) then
                               np = np + 1
							   i0=i+ii !min(ni,max(1,i+ii))
							   j0=j+jj !min(nj,max(1,j+jj))
							   k0=k+kk !min(nk,max(1,k+kk))
                               rm = rm + r(i0,j0,k0)
                               um = um + u(i0,j0,k0)
                               vm = vm + v(i0,j0,k0)
                               wm = wm + w(i0,j0,k0)
                               pm = pm + p(i0,j0,k0)
                            endif
                       enddo
                       enddo
                       enddo

				     rnp=real(np)
                   rm = rm/rnp
                   um = um/rnp
                   vm = vm/rnp
                   wm = wm/rnp
                   pm = pm/rnp
!				   rm = min(r_max,  max(r_min,rm))
!				   um = min(2.,     max(-2.  ,um))
!				   vm = min(2.,     max(-2.  ,vm))
!				   wm = min(2.,     max(-2.  ,wm))
!				   pm = min(p_max,  max(p_min,pm))

!				   rm = min(5.*roo,  max(0.2*roo,rm))
!				   um = min(2.,     max(-2.    ,um))
!				   vm = min(2.,     max(-2.    ,vm))
!				   wm = min(2.,     max(-2.    ,wm))
!				   pm = min(5.*poo,  max(0.2*poo,pm))

                else

                   !if(n_count < 3 )write(*,902)'p<p_core at:',nb,i,j,k,pm,rm
902                        format(1x,a13,i3,3i4,2e14.5)
				           n_count = n_count + 1
                   rm = 0.0
                   um = 0.0
                   vm = 0.0
                   wm = 0.0
                   pm = 0.0
                   np = 0
				       do ii=-1,1
                       do jj =-1,1
                       do kk =-1,1
                            if ( ( abs(ii) +abs(jj) +abs(kk) ) == 1 ) then
                               np = np + 1
							   i0=i+ii !min(ni,max(1,i+ii))
							   j0=j+jj !min(nj,max(1,j+jj))
							   k0=k+kk !min(nk,max(1,k+kk))
                               rm = rm + r(i0,j0,k0)
                               um = um + u(i0,j0,k0)
                               vm = vm + v(i0,j0,k0)
                               wm = wm + w(i0,j0,k0)
                               pm = pm + p(i0,j0,k0)
                            endif
                       enddo
                       enddo
                       enddo
					  rnp=real(np)
                   rm = rm/rnp
                   um = um/rnp
                   vm = vm/rnp
                   wm = wm/rnp
                   pm = pm/rnp

!				   rm = min(r_max,  max(r_min,rm))
!				   um = min(2.,     max(-2.  ,um))
!				   vm = min(2.,     max(-2.  ,vm))
!				   wm = min(2.,     max(-2.  ,wm))
!				   pm = min(p_max,  max(p_min,pm))

!				   rm = min(5.*roo,  max(0.2*roo,rm))
!				   um = min(2.,     max(-2.    ,um))
!				   vm = min(2.,     max(-2.    ,vm))
!				   wm = min(2.,     max(-2.    ,wm))
!				   pm = min(5.*poo,  max(0.2*poo,pm))

                endif

                em = pm/(gama-1.0) + 0.5*rm*( um*um + vm*vm + wm*wm )
                q(1,i,j,k) = rm
                q(2,i,j,k) = rm * um
                q(3,i,j,k) = rm * vm
                q(4,i,j,k) = rm * wm
                q(5,i,j,k) = em
             endif
             
             r(i,j,k) = rm
             u(i,j,k) = um
             v(i,j,k) = vm
             w(i,j,k) = wm
             p(i,j,k) = pm
          enddo
       enddo
    enddo
    !if(n_count > 0 )write(*,*)'��',nb,'������쳣����Ϊ',n_count
    
    if (n_count > 0) call beep_warning(500,100)
	
	if(n_count > ni*nj*nk/64 .or. n_count >= 500)then
		write(*,*)nb, '����쳣���������涨ֵ'
		write(*,*)'================����ǿ�ƽ���==============='
		stop
	endif
	

    contains
    
    subroutine beep_warning(frequency, duration)
       !use ifport
       integer(4),intent(in) :: frequency, duration
       !call beepqq(frequency, duration)
    end subroutine beep_warning

end subroutine update_ns
!_____________________________________________________________________!
subroutine update_limited(nb) !��������������
    use global_const,only : pmin_limit,pmax_limit,rmax_limit,rmin_limit
    use global_variables, &
    only : method,nl,ni,nj,nk,nm,gama,q,dq,r,u,v,w,p,roo,uoo,voo,woo,poo,moo
    implicit none
    integer :: nb,m,m1,i,j,k,n_count,i0,j0,k0
    integer :: ii,jj,kk,np
	integer :: iend,jend,kend !* Added by TU Guohua
    real :: gamam1,p_min,r_min,p_max,r_max
    real :: rm,um,vm,wm,pm,em,rm1,rnp
	real :: qtem(nm),Falpha

    gamam1 = gama - 1.0
    
    p_min = pmin_limit
    p_max = pmax_limit
    
    r_min = rmin_limit
    r_max = rmax_limit

  	n_count = 0

!    m1 = 1 + method
!*TGH. ע����ʼ����
   m1=1
   iend=ni
   jend=nj
   kend=nk
    do k=m1,kend 
    do j=m1,jend    
    do i=m1,iend
								 
        do m=1,nm
             qtem(m) = q(m,i,j,k) + dq(m,i,j,k)
        enddo
        rm = qtem(1)
        rm1 = 1.0/rm
        um = qtem(2) * rm1
        vm = qtem(3) * rm1
        wm = qtem(4) * rm1
        em = qtem(5)
        pm = gamam1 * ( em - 0.5*rm*( um*um + vm*vm + wm*wm ) )
		
		Falpha = 1./( 1.+ 2.*max(0.0,-0.3+abs((pm-p(i,j,k))/p(i,j,k))) )

        do m=1,nm
             qtem(m) = q(m,i,j,k) + dq(m,i,j,k)*Falpha

			 q(m,i,j,k) = qtem(m)
        enddo
        rm = qtem(1)
        rm1 = 1.0/rm
        um = qtem(2) * rm1
        vm = qtem(3) * rm1
        wm = qtem(4) * rm1
        em = qtem(5)
        pm = gamam1 * ( em - 0.5*rm*( um*um + vm*vm + wm*wm ) )

        if ( pm <= p_min .or. rm <= r_min  .or. &
		  	    & pm >= p_max .or. rm >= r_max   ) then  !________

                if ( pm <= 0.0 .or. rm <= 0.0 ) then
                   !if(n_count < 3 ) write(*,901)'p<0 at:',nb,i,j,k,pm,rm
901                                     format(1x,a7,i3,3i4,2e14.5)
				   
				   n_count = n_count + 1

                   rm = 0.0
                   um = 0.0
                   vm = 0.0
                   wm = 0.0
                   pm = 0.0
                   np = 0
				       do ii=-1,1
                       do jj =-1,1
                       do kk =-1,1
                            if ( ( abs(ii) +abs(jj) +abs(kk) ) == 1 ) then
                               np = np + 1
							   i0=i+ii !min(ni,max(1,i+ii))
							   j0=j+jj !min(nj,max(1,j+jj))
							   k0=k+kk !min(nk,max(1,k+kk))
                               rm = rm + r(i0,j0,k0)
                               um = um + u(i0,j0,k0)
                               vm = vm + v(i0,j0,k0)
                               wm = wm + w(i0,j0,k0)
                               pm = pm + p(i0,j0,k0)
                            endif
                       enddo
                       enddo
                       enddo

				     rnp=real(np)
                   rm = rm/rnp
                   um = um/rnp
                   vm = vm/rnp
                   wm = wm/rnp
                   pm = pm/rnp
!				   rm = min(r_max,  max(r_min,rm))
!				   um = min(2.,     max(-2.  ,um))
!				   vm = min(2.,     max(-2.  ,vm))
!				   wm = min(2.,     max(-2.  ,wm))
!				   pm = min(p_max,  max(p_min,pm))

!				   rm = min(5.*roo,  max(0.2*roo,rm))
!				   um = min(2.,     max(-2.    ,um))
!				   vm = min(2.,     max(-2.    ,vm))
!				   wm = min(2.,     max(-2.    ,wm))
!				   pm = min(5.*poo,  max(0.2*poo,pm))

                else

                   !if(n_count < 3 )write(*,902)'p<p_core at:',nb,i,j,k,pm,rm
902                        format(1x,a13,i3,3i4,2e14.5)
				           n_count = n_count + 1
                   rm = 0.0
                   um = 0.0
                   vm = 0.0
                   wm = 0.0
                   pm = 0.0
                   np = 0
				       do ii=-1,1
                       do jj =-1,1
                       do kk =-1,1
                            if ( ( abs(ii) +abs(jj) +abs(kk) ) == 1 ) then
                               np = np + 1
							   i0=i+ii !min(ni,max(1,i+ii))
							   j0=j+jj !min(nj,max(1,j+jj))
							   k0=k+kk !min(nk,max(1,k+kk))
                               rm = rm + r(i0,j0,k0)
                               um = um + u(i0,j0,k0)
                               vm = vm + v(i0,j0,k0)
                               wm = wm + w(i0,j0,k0)
                               pm = pm + p(i0,j0,k0)
                            endif
                       enddo
                       enddo
                       enddo
					  rnp=real(np)
                   rm = rm/rnp
                   um = um/rnp
                   vm = vm/rnp
                   wm = wm/rnp
                   pm = pm/rnp

!				   rm = min(r_max,  max(r_min,rm))
!				   um = min(2.,     max(-2.  ,um))
!				   vm = min(2.,     max(-2.  ,vm))
!				   wm = min(2.,     max(-2.  ,wm))
!				   pm = min(p_max,  max(p_min,pm))

!				   rm = min(5.*roo,  max(0.2*roo,rm))
!				   um = min(2.,     max(-2.    ,um))
!				   vm = min(2.,     max(-2.    ,vm))
!				   wm = min(2.,     max(-2.    ,wm))
!				   pm = min(5.*poo,  max(0.2*poo,pm))

                endif

                em = pm/(gama-1.0) + 0.5*rm*( um*um + vm*vm + wm*wm )
                q(1,i,j,k) = rm
                q(2,i,j,k) = rm * um
                q(3,i,j,k) = rm * vm
                q(4,i,j,k) = rm * wm
                q(5,i,j,k) = em
             endif

             r(i,j,k) = rm
             u(i,j,k) = um
             v(i,j,k) = vm
             w(i,j,k) = wm
             p(i,j,k) = pm
          enddo
       enddo
    enddo
    !if(n_count > 0 )write(*,*)'��',nb,'������쳣����Ϊ',n_count
    
    if (n_count > 0) call beep_warning(500,100)
	
	if(n_count > ni*nj*nk/64 .or. n_count >= 500)then
		write(*,*)nb, '����쳣���������涨ֵ'
		write(*,*)'================����ǿ�ƽ���==============='
		stop
	endif


    contains
    
    subroutine beep_warning(frequency, duration)
       !use ifport
       integer(4),intent(in) :: frequency, duration
       !call beepqq(frequency, duration)
    end subroutine beep_warning

end subroutine update_limited
!_____________________________________________________________________!


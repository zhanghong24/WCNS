!-----------------------------------------------------------------------------!
subroutine set_nonconnect_boundary_dq_0(nb)   !___ֻ�����Խӱ߽�
    use global_variables
    implicit none
    integer :: nb,nr,nrmax,bctype,m,k

    nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
    !���崦���߽�����
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if( bctype > 0 ) then               !�ǶԽӱ߽�
          call boundary_dq_to_0(nb,nr,bctype)
       endif
    enddo
		
    return
end subroutine set_nonconnect_boundary_dq_0

!_____________________________________________________________________!
subroutine store_time(nb)
    use global_variables
    implicit none
    integer :: i,j,k,m,cmethod,nb

    do k=1,nk
       do j=1,nj
          do i=1,ni
             mb_dtdt(nb)%a3d(i,j,k) = dtdt(i,j,k)
          enddo
       enddo
    enddo

    return
end subroutine store_time
!_____________________________________________________________________!
subroutine bc_average  !���޲�ֶԽӱ߽�ֵ
    use global_variables
    implicit none
    integer :: s_st(3),s_ed(3),s_lr3d(3),nb,nbt,nr,bctype,nrmax
    integer :: i,j,k,is,js,ks,it,jt,kt,m
    real :: prim_s1(nl),q_1(nl)
    real :: zf1

    do nb = 1, nblocks

    nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
    !���崦���߽�����
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if( bctype < 0 ) then               !�Խӱ߽�

		   nbt = mb_bc(nb)%bc(nr)%nbt             !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������

        do m=1,3
           s_st(m)    = mb_bc(nb)%bc(nr)%s_st(m)
           s_ed(m)    = mb_bc(nb)%bc(nr)%s_ed(m)
           s_lr3d(m)  = mb_bc(nb)%bc(nr)%s_lr3d(m)
        enddo
        do i = s_st(1),s_ed(1)
           do j = s_st(2),s_ed(2)
              do k = s_st(3),s_ed(3)
                 is = i
                 js = j
                 ks = k
                 it = mb_bc(nb)%bc(nr)%image(i,j,k )
                 jt = mb_bc(nb)%bc(nr)%jmage(i,j,k )
                 kt = mb_bc(nb)%bc(nr)%kmage(i,j,k )

                 mb_r(nb)%a3d(i,j,k) = 0.5*( mb_r(nb)%a3d(is,js,ks) + mb_r(nbt)%a3d(it,jt,kt) )
                 mb_u(nb)%a3d(i,j,k) = 0.5*( mb_u(nb)%a3d(is,js,ks) + mb_u(nbt)%a3d(it,jt,kt) )
                 mb_v(nb)%a3d(i,j,k) = 0.5*( mb_v(nb)%a3d(is,js,ks) + mb_v(nbt)%a3d(it,jt,kt) )
                 mb_w(nb)%a3d(i,j,k) = 0.5*( mb_w(nb)%a3d(is,js,ks) + mb_w(nbt)%a3d(it,jt,kt) )
                 mb_p(nb)%a3d(i,j,k) = 0.5*( mb_p(nb)%a3d(is,js,ks) + mb_p(nbt)%a3d(it,jt,kt) )
                 mb_t(nb)%a3d(i,j,k) = 0.5*( mb_t(nb)%a3d(is,js,ks) + mb_t(nbt)%a3d(it,jt,kt) )

                 mb_r(nbt)%a3d(it,jt,kt) =  mb_r(nb)%a3d(is,js,ks)
                 mb_u(nbt)%a3d(it,jt,kt) =  mb_u(nb)%a3d(is,js,ks)
                 mb_v(nbt)%a3d(it,jt,kt) =  mb_v(nb)%a3d(is,js,ks)
                 mb_w(nbt)%a3d(it,jt,kt) =  mb_w(nb)%a3d(is,js,ks)
                 mb_p(nbt)%a3d(it,jt,kt) =  mb_p(nb)%a3d(is,js,ks)
                 mb_t(nbt)%a3d(it,jt,kt) =  mb_t(nb)%a3d(is,js,ks)

              enddo
           enddo
        enddo

		endif
	enddo
	enddo

    return
end subroutine bc_average
!_____________________________________________________________________!
subroutine Aps_dq(prim,nxyz,s_lr,dq,f,npn)
!*TGH. ��A~��A+������ʾ��Ȼ��ֱ�Ӷ�RHS����dq������ת��
!*TGH.  f = A^-.Dq ��npn=1ʱ�� ��f = A^+.Dq ��npn=-1ʱ��

    use global_const,only:nl,ns,ms1,beta1,tref,gama,sml_sss
    implicit none
    integer :: npn,m,s_lr
    real :: hint
    real :: prim(nl),dq(nl),hs(ns),as(ns)
    real :: nx,ny,nz,nt,ct,cgm,cgm1,nxyz(4),f(nl)
    real :: l1,l4,l5,x1,x2
    real :: dh,dc,c2dc,ae,af
    real :: rm,um,vm,wm,pm,cm,c2,v2,tm,hm

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

    l1 = ct
    l4 = ct + cm * cgm
    l5 = ct - cm * cgm

    l1 = l1*s_lr
    l4 = l4*s_lr
    l5 = l5*s_lr

    l1 = 0.5*( 1.0 + npn*sign(1.0,l1) )
    l4 = 0.5*( 1.0 + npn*sign(1.0,l4) )
    l5 = 0.5*( 1.0 + npn*sign(1.0,l5) )

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
end subroutine Aps_dq
!_____________________________________________________________________!
subroutine pri_to_con(nb)
    use global_variables
    implicit none
    integer :: nb,m,m1,i,j,k,n_count
    real :: rm,um,vm,wm,pm,em,gamam1

    gamam1 = gama - 1.0
    do k=-1,nk+1  
      do j=-1,nj+1   
         do i=-1,ni+1  
             rm = mb_r(nb)%a3d(i,j,k)
             um = mb_u(nb)%a3d(i,j,k)
             vm = mb_v(nb)%a3d(i,j,k)
             wm = mb_w(nb)%a3d(i,j,k)
             pm = mb_p(nb)%a3d(i,j,k)
             mb_q(nb)%a4d(1,i,j,k) = rm
             mb_q(nb)%a4d(2,i,j,k) = rm*um
             mb_q(nb)%a4d(3,i,j,k) = rm*vm
             mb_q(nb)%a4d(4,i,j,k) = rm*wm
             mb_q(nb)%a4d(5,i,j,k) = pm/gamam1 + 0.5*rm*( um*um + vm*vm + wm*wm )
          enddo
       enddo
    enddo

    return
end subroutine pri_to_con
!-----------------------------------------------------------------------------!
subroutine update_cic(nb)   !___ֻ�����Խӱ߽�
!*TGH. �Ѿ�������˶Խӱ߽��ϵ�dq���غ������ʱ��������������RHS�������±߽��ϵ�ԭʼ����
   use global_variables
    implicit none
    integer :: nb,nr,nrmax,bctype,m,k

    nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
    !���崦���߽�����
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if( bctype < 0 ) then               !�Խӱ߽�
          call update_ns_cic(nb,nr)   !*tgh. ���¶Խӱ߽��ϵ�ֵ
!*TGH. �Ѿ�������˶Խӱ߽��ϵ�dq���غ������ʱ��������������RHS�������±߽��ϵ�ԭʼ����
       endif
    enddo
		
    return
end subroutine update_cic
!_____________________________________________________________________!
subroutine update_ns_cic(nb,nr)   !___���¶Խӱ߽��q
!*TGH. �Ѿ�������˶Խӱ߽��ϵ�dq���غ������ʱ��������������RHS�������±߽��ϵ�ԭʼ����
    use global_variables
    integer :: nb,nr,bctype,m,n,n1
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,i,j,k,nt,ntarg
    integer :: it0,jt0,kt0,it,jt,kt
    integer :: s_st(3),s_ed(3)
    integer :: nv_s,nv_t,nfile_par_s,nfile_par_t,num_var_s,num_var_t,nchem_s,nchem_t
    real    :: prim_s1(nl),q_1(nl),zf1,ttt1,f(nl)

    integer :: m1,n_count
    integer :: ii,jj,kk,np,i0,j0,k0
    real :: gamam1,p_core,d_core
    real :: rm,um,vm,wm,pm,em,rm1,rnp
    gamam1 = gama - 1.0
    p_core = 0.0000005
    d_core = 0.0000010
  	n_count = 0

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
    enddo
    s_nd  = mb_bc(nb)%bc(nr)%s_nd            !�߽��淽��:1,2,3��Ӧ��i,j,k
    s_lr  = mb_bc(nb)%bc(nr)%s_lr            !���ұ߽�-1,1��Ӧ�����ұ߽�
    s_fix = mb_bc(nb)%bc(nr)%s_fix           !�̶�����(fixed_coor)
    nbs   = mb_bc(nb)%bc(nr)%nbs             !��� 

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             do m=1,nm
               q(m,i,j,k) = q(m,i,j,k) + dq(m,i,j,k)
             enddo
          enddo
       enddo
    enddo

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)

             rm = q(1,i,j,k)
             rm1 = 1.0/rm
             um = q(2,i,j,k) * rm1
             vm = q(3,i,j,k) * rm1
             wm = q(4,i,j,k) * rm1
             em = q(5,i,j,k)
             pm = gamam1 * ( em - 0.5*rm*( um*um + vm*vm + wm*wm ) )
             if ( pm <= p_core .or. rm <= d_core) then  !________
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
                            if ( ( abs(ii) +abs(jj) +abs(kk) ) /= 0 ) then
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
                   rm = ABS(rm)/rnp
                   um = um/rnp
                   vm = vm/rnp
                   wm = wm/rnp
                   pm = ABS(pm)/rnp
				   rm = min(5.*roo, max(0.2*roo,rm))
				   um = min(2.,     max(-2.    ,um))
				   vm = min(2.,     max(-2.    ,vm))
				   wm = min(2.,     max(-2.    ,wm))
				   pm = min(5.*poo, max(0.2*poo,pm))


!                   rm = 1.
!                   um = 1.
!                   vm = 0.
!                   wm = 0.
!                   pm = poo

                else  !*tgh. end ѹ���ܶȳ�������

                   !if(n_count < 3 )!write(*,902)'p<p_core at:',nb,i,j,k,pm,rm
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
                            if ( ( abs(ii) +abs(jj) +abs(kk) ) /= 0 ) then
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
                   rm = ABS(rm)/rnp
                   um = um/rnp
                   vm = vm/rnp
                   wm = wm/rnp
                   pm = ABS(pm)/rnp
				   rm = min(5.*roo, max(0.2*roo,rm))
				   um = min(2.,     max(-2.    ,um))
				   vm = min(2.,     max(-2.    ,vm))
				   wm = min(2.,     max(-2.    ,wm))
				   pm = min(5.*poo, max(0.2*poo,pm))


!                   rm = 1.
!                   um = 1.
!                   vm = 0.
!                   wm = 0.
!                   pm = poo

                  endif !*tgh. end ѹ���ܶȼ�С����

                em = pm/(gama-1.0) + 0.5*rm*( um*um + vm*vm + wm*wm )
                q(1,i,j,k) = rm
                q(2,i,j,k) = rm * um
                q(3,i,j,k) = rm * vm
                q(4,i,j,k) = rm * wm
                q(5,i,j,k) = em
             endif !*tgh. end ѹ���ܶ�����

             r(i,j,k) = rm
             u(i,j,k) = um
             v(i,j,k) = vm
             w(i,j,k) = wm
             p(i,j,k) = pm
          enddo
       enddo
    enddo
    !if(n_count > 0 )write(*,*)'��',nb,'������쳣����Ϊ',n_count

    return
end subroutine update_ns_cic
!-----------------------------------------------------------------------------!
subroutine boundary_connect1(nb)   !___ֻ�����Խӱ߽�
    use global_variables
    implicit none
    integer :: nb,nr,nrmax,bctype,m,k

    nrmax = mb_bc(nb)%nregions             !���鹲��nrmax���߽���Ҫ����
    !���崦���߽�����
    do nr = 1,nrmax
       bctype = mb_bc(nb)%bc(nr)%bctype
       if( bctype < 0 ) then               !�Խӱ߽�
          call boundary_n1_nocic(nb,nr,bctype)
          !*TGH. �Խӱ߽�ʱ������ϵ�ֵ
          !*TGH. ����ԭʼ���� �� ��Ӧ������������һ����ֵ
          !*TGH. ����������ֱ��ȡԴ��Խ����ϵ�������
       endif
    enddo
		
    return
end subroutine boundary_connect1
!_____________________________________________________________________!
subroutine boundary_n1_nocic(nb,nr,bctype) !�Խӱ߽����� ���������Խӱ߽� ֻ�������޲��
!*TGH. �Խӱ߽�ʱ������ϵ�ֵ
!*TGH. ����ԭʼ���� �� ��Ӧ������������һ����ֵ
!*TGH. ����������ֱ��ȡԴ��Խ����ϵ�������

    use global_variables
    implicit none
    integer :: nb,nr,bctype,m,n,n1
    integer :: nbs,s_nd,s_fix,s_lr
    integer :: nbt,t_nd,t_fix,t_lr
    integer :: is,js,ks,i,j,k,nt,ntarg
    integer :: it0,jt0,kt0,it,jt,kt
    integer :: s_st(3),s_ed(3)
    integer :: nv_s,nv_t,nfile_par_s,nfile_par_t,num_var_s,num_var_t,nchem_s,nchem_t
    real    :: prim_s1(nl),q_1(nl),zf1,ttt1,f(nl)

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

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
             do n=1,1
                is = i + mb_bc(nb)%bc(nr)%s_lr3d(1)*n
                js = j + mb_bc(nb)%bc(nr)%s_lr3d(2)*n  !*TGH. ��һ��
                ks = k + mb_bc(nb)%bc(nr)%s_lr3d(3)*n
                nt = n - 1 + method
                it0= mb_bc(nb)%bc(nr)%image(i,j,k )
                jt0= mb_bc(nb)%bc(nr)%jmage(i,j,k )    !*TGH. Ŀ����ϵ�����
                kt0= mb_bc(nb)%bc(nr)%kmage(i,j,k )
                it = it0 - mb_bc(nb)%bc(nr)%t_lr3d(1)*nt
                jt = jt0 - mb_bc(nb)%bc(nr)%t_lr3d(2)*nt  !*TGH. Ŀ�����������һ��
                kt = kt0 - mb_bc(nb)%bc(nr)%t_lr3d(3)*nt

                mb_r(nbs)%a3d(is,js,ks) = mb_r(nbt)%a3d(it,jt,kt)
                mb_u(nbs)%a3d(is,js,ks) = mb_u(nbt)%a3d(it,jt,kt)
                mb_v(nbs)%a3d(is,js,ks) = mb_v(nbt)%a3d(it,jt,kt)  !*TGH. ����ԭʼ���� �� ��Ӧ������������һ����ֵ
                mb_w(nbs)%a3d(is,js,ks) = mb_w(nbt)%a3d(it,jt,kt)
                mb_p(nbs)%a3d(is,js,ks) = mb_p(nbt)%a3d(it,jt,kt)

                mb_kct(nbs)%a3d(is,js,ks) = mb_kct(nbs)%a3d(i,j,k)  
                mb_kcx(nbs)%a3d(is,js,ks) = mb_kcx(nbs)%a3d(i,j,k) !*TGH. ����������ֱ��ȡԴ��Խ����ϵ�������
                mb_kcy(nbs)%a3d(is,js,ks) = mb_kcy(nbs)%a3d(i,j,k)  
                mb_kcz(nbs)%a3d(is,js,ks) = mb_kcz(nbs)%a3d(i,j,k)  

                mb_ett(nbs)%a3d(is,js,ks) = mb_ett(nbs)%a3d(i,j,k)  
                mb_etx(nbs)%a3d(is,js,ks) = mb_etx(nbs)%a3d(i,j,k)  
                mb_ety(nbs)%a3d(is,js,ks) = mb_ety(nbs)%a3d(i,j,k)  
                mb_etz(nbs)%a3d(is,js,ks) = mb_etz(nbs)%a3d(i,j,k)  

                mb_ctt(nbs)%a3d(is,js,ks) = mb_ctt(nbs)%a3d(i,j,k)  
                mb_ctx(nbs)%a3d(is,js,ks) = mb_ctx(nbs)%a3d(i,j,k)  
                mb_cty(nbs)%a3d(is,js,ks) = mb_cty(nbs)%a3d(i,j,k)  
                mb_ctz(nbs)%a3d(is,js,ks) = mb_ctz(nbs)%a3d(i,j,k)
								  
                mb_vol(nbs)%a3d(is,js,ks) = mb_vol(nbs)%a3d(i,j,k)  !*TGH. ����Jacobianֱ��ȡԴ��Խ����ϵ�Jacobian

             enddo
          enddo
       enddo
    enddo

     
    return
end subroutine boundary_n1_nocic
!_____________________________________________________________________!
subroutine boundary_dq_to_0(nb,nr,bctype) !�Խӱ߽����� ���������Խӱ߽� ֻ�������޲��
    use global_variables
    implicit none
    integer :: nb,nr,bctype,m,n,n1
    integer :: i,j,k
    integer :: s_st(3),s_ed(3)
    real    :: prim_s1(nl),q_1(nl),zf1,ttt1,f(nl)

    do m=1,3
       s_st(m) = mb_bc(nb)%bc(nr)%s_st(m)    !��ʼ������(�ɶ�����)
       s_ed(m) = mb_bc(nb)%bc(nr)%s_ed(m)    !��ֹ������(�ɶ�����)
    enddo

    do i = s_st(1),s_ed(1)
       do j = s_st(2),s_ed(2)
          do k = s_st(3),s_ed(3)
                do n1=1,5
								dq(n1,i,j,k) = 0.d0
				enddo
          enddo
       enddo
    enddo
     
    return
end subroutine boundary_dq_to_0

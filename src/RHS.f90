
subroutine r_h_s_vis
    use define_precision_mod
    use global_variables,only : nblocks,nscheme,nvis,nlamtur,csrv,nchem,nchem_source
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
    REAL :: TST,TED
    integer :: nb,pnb


#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb=1,nblocks
#endif
        call recast_grid(nb)
        call recast_field(nb)

	    if ( nlamtur >= 0 ) call turbulence(nb) !��������

        !*TGH. ���濪ʼճ������ɢ
         call WCNSE5_VIS_VIRTUAL
          !*TGH. end ճ������ɢ

	enddo !nb = 1,nblocks

end subroutine r_h_s_vis

!=============================================================================!
!=============================================================================!

subroutine r_h_s_invis
    use define_precision_mod
    use global_variables,only : nblocks,nscheme,nlimiter,nvis
#ifdef PARALLEL
    use mod_parallels,only : pnblocks,pnbindexs
#endif
    implicit none
	integer :: nb,pnb
    real,external :: f_zero,minmod,vanleer,min_3u,min_3q,min_van

#ifdef PARALLEL
    do pnb=1,pnblocks
       nb = pnbindexs(pnb)
#else
    do nb = 1,nblocks
#endif

        call recast_grid(nb)
        call recast_field(nb) !*TGH. ���¶Ⱥ�����Ҳ������recast

		if ( nlimiter == 0 ) then
		   call invcode(f_zero)
		elseif ( nlimiter == 1 ) then
           call invcode(minmod)
		elseif ( nlimiter == 2 ) then
           call invcode(vanleer)
		elseif ( nlimiter == 3 ) then
           call invcode(min_3u)
		elseif ( nlimiter == 4 ) then
           call invcode(min_van)
		endif

	enddo !* nb = 1,nblocks

    return


end subroutine r_h_s_invis

!=============================================================================!
subroutine get_q_only
    use define_precision_mod
    use global_variables,only : r,u,v,w,p,q,dq,gama, &
	                          & nl,ni,nj,nk,method
    implicit none
    integer :: i,j,k,m,ijk
    integer :: nn(2)
    real :: prim(nl),q_q(nl)
    real :: a2
    real :: rm,um,vm,wm,pm,em

    do k=1,nk
       do j=1,nj
          do i=1,ni
             prim(1) = r(i,j,k)
             prim(2) = u(i,j,k)
             prim(3) = v(i,j,k)
             prim(4) = w(i,j,k)
             prim(5) = p(i,j,k)
             call prim_to_q(prim,q_q,gama)

			 do m=1,nl
                q(m,i,j,k) = q_q(m)
		     enddo

          enddo
       enddo
    enddo

    nn(1) = 0
    nn(2) = ni+1

    do k=1,nk
       do j=1,nj
          do i=1,2
		     ijk=nn(i)
             prim(1) = r(ijk,j,k)
             prim(2) = u(ijk,j,k)
             prim(3) = v(ijk,j,k)
             prim(4) = w(ijk,j,k)
             prim(5) = p(ijk,j,k)
             call prim_to_q_bc(prim,q_q,gama)   !prim_to_q(prim,q_q,gama)

			 do m=1,nl
                q(m,nn(i),j,k) = q_q(m)
			 enddo

          enddo
       enddo
    enddo
    nn(2) = nj+1
    do k=1,nk
       do i=1,ni
          do j=1,2
		     ijk=nn(j)
             prim(1) = r(i,ijk,k)
             prim(2) = u(i,ijk,k)
             prim(3) = v(i,ijk,k)
             prim(4) = w(i,ijk,k)
             prim(5) = p(i,ijk,k)
             call prim_to_q_bc(prim,q_q,gama)   !prim_to_q(prim,q_q,gama)

			 do m=1,nl
                q(m,i,nn(j),k) = q_q(m)
			 enddo

          enddo
       enddo
    enddo

    nn(2) = nk+1
    do j=1,nj
       do i=1,ni
          do k=1,2
		     ijk=nn(k)
             prim(1) = r(i,j,ijk)
             prim(2) = u(i,j,ijk)
             prim(3) = v(i,j,ijk)
             prim(4) = w(i,j,ijk)
             prim(5) = p(i,j,ijk)
             call prim_to_q_bc(prim,q_q,gama)   !prim_to_q(prim,q_q,gama)

			 do m=1,nl
                q(m,i,j,nn(k)) = q_q(m)
			 enddo

          enddo
       enddo
    enddo

    return
end subroutine get_q_only

!=============================================================================!
!=============================================================================!

subroutine get_c_only
!*tgh. ֻ�����٣������¶�
    use global_variables,only : r,p,c, &
	                          & nl,ni,nj,nk,gama,method
    implicit none
    integer :: i,j,k
    integer :: nn(2)
    real :: a2
    real :: rm,pm

    do k=1,nk
       do j=1,nj
          do i=1,ni
             rm = r(i,j,k)
             pm = p(i,j,k)

             a2 = gama * pm / rm
             c(i,j,k) = sqrt( a2 )

          enddo
       enddo
    enddo

    nn(1) = 0
    nn(2) = ni+1

    do k=1,nk
       do j=1,nj
          do i=1,2
             rm = r(nn(i),j,k)
             pm = p(nn(i),j,k)
             a2 = gama * pm /  rm
             c(nn(i),j,k) = sqrt( a2 )
          enddo
       enddo
    enddo

    nn(2) = nj+1
    do k=1,nk
       do i=1,ni
          do j=1,2
             rm = r(i,nn(j),k)
             pm = p(i,nn(j),k)
             a2 = gama * pm /  rm
             if ( a2 <= 0.0 ) then
                write(*,*)i,j,k,' p=',p(i,nn(j),k),' r=',rm
             endif
             c(i,nn(j),k) = sqrt( a2 )

          enddo
       enddo
    enddo

    nn(2) = nk+1
    do j=1,nj
       do i=1,ni
          do k=1,2
             rm = r(i,j,nn(k))
             pm = p(i,j,nn(k))
             a2 = gama * pm /  rm
             c(i,j,nn(k)) = sqrt( a2 )
          enddo
       enddo
    enddo

    return
end subroutine get_c_only

!=============================================================================!
!=============================================================================!
subroutine get_t_only
    use global_variables,only : r,p,t,moo,&
	                          & nl,ni,nj,nk,gama,method
    implicit none
    integer :: i,j,k
    integer :: nn(2)
    real :: a2,moo2
    real :: rm,pm

	moo2=moo*moo
    do k=1,nk
       do j=1,nj
          do i=1,ni
             rm = r(i,j,k)
             pm = p(i,j,k)

             a2 = gama * pm / rm
             t(i,j,k) = moo2 * a2

          enddo
       enddo
    enddo

    nn(1) = 0
    nn(2) = ni+1

    do k=1,nk
       do j=1,nj
          do i=1,2
             rm = r(nn(i),j,k)
             pm = p(nn(i),j,k)
             a2 = gama * pm /  rm
             t(nn(i),j,k) = moo2 * a2
          enddo
       enddo
    enddo

    nn(2) = nj+1
    do k=1,nk
       do i=1,ni
          do j=1,2
             rm = r(i,nn(j),k)
             pm = p(i,nn(j),k)
             a2 = gama * pm /  rm
             t(i,nn(j),k) = moo2 * a2

          enddo
       enddo
    enddo

    nn(2) = nk+1
    do j=1,nj
       do i=1,ni
          do k=1,2
             rm = r(i,j,nn(k))
             pm = p(i,j,nn(k))
             a2 = gama * pm /  rm
             t(i,j,nn(k)) = moo2 * a2
          enddo
       enddo
    enddo

    return
end subroutine get_t_only
!=============================================================================!
!=============================================================================!
subroutine get_c_t_only
    use global_variables,only : r,p,c,t,moo,&
	                          & nl,ni,nj,nk,gama,method
    implicit none
    integer :: i,j,k
    integer :: n1(2),n2(2),n3(2)
    real :: a2,moo2
    real :: rm,pm

	moo2=moo*moo
    do k=1,nk
       do j=1,nj
          do i=1,ni
             rm = r(i,j,k)
             pm = p(i,j,k)

             a2 = gama * pm / rm
			 c(i,j,k) = sqrt(a2)
             t(i,j,k) = moo2 * a2

          enddo
       enddo
    enddo

    n3(1) = -2
    n2(1) = -1
    n1(1) = 0
    n1(2) = ni+1
    n2(2) = ni+2
    n3(2) = ni+3
    do k=1,nk
       do j=1,nj
          do i=1,2
             rm = r(n1(i),j,k)
             pm = p(n1(i),j,k)
             a2 = gama * pm /  rm
			 c(n1(i),j,k) = sqrt(a2)
             t(n1(i),j,k) = moo2 * a2

             rm = r(n2(i),j,k) 
             pm = p(n2(i),j,k)
             a2 = gama * pm /  rm
             t(n2(i),j,k) = moo2 * a2
                          
             rm = abs( r(n3(i),j,k) )
             pm =      p(n3(i),j,k)
             t(n3(i),j,k) = moo2 * a2
          enddo
          enddo
       enddo

    n1(2) = nj+1
    n2(2) = nj+2
    n3(2) = nj+3
    do k=1,nk
       do i=1,ni
          do j=1,2
             rm = r(i,n1(j),k)
             pm = p(i,n1(j),k)
             a2 = gama * pm /  rm
			 c(i,n1(j),k) = sqrt(a2)
             t(i,n1(j),k) = moo2 * a2

             rm = r(i,n2(j),k)
             pm = p(i,n2(j),k)
             a2 = gama * pm /  rm
             t(i,n2(j),k) = moo2 * a2

             rm = abs( r(i,n3(j),k) )
             pm =      p(i,n3(j),k)
             a2 = gama * pm /  rm
             t(i,n3(j),k) = moo2 * a2
          enddo
       enddo
    enddo

    n1(2) = nk+1
    n2(2) = nk+2
    n3(2) = nk+3
    do j=1,nj
       do i=1,ni
          do k=1,2
             rm = r(i,j,n1(k))
             pm = p(i,j,n1(k))
             a2 = gama * pm /  rm
			 c(i,j,n1(k)) = sqrt(a2)
             t(i,j,n1(k)) = moo2 * a2

             rm = r(i,j,n2(k))
             pm = p(i,j,n2(k))
             a2 = gama * pm /  rm
             t(i,j,n2(k)) = moo2 * a2
             
             rm = abs( r(i,j,n3(k)) )
             pm =      p(i,j,n3(k))
             a2 = gama * pm /  rm
             t(i,j,n3(k)) = moo2 * a2
          enddo
       enddo
    enddo

    return
end subroutine get_c_t_only
!=============================================================================!
!=============================================================================!
subroutine get_t_ini
    use global_variables,only : r,p,c,t,moo,&
	                          & nl,ni,nj,nk,gama,method
    implicit none
    integer :: i,j,k
    integer :: nn(4)
    real :: a2,moo2
    real :: rm,pm

	moo2=moo*moo
    do k=1,nk
       do j=1,nj
          do i=1,ni
             rm = r(i,j,k)
             pm = p(i,j,k)

             a2 = gama * pm / rm
             t(i,j,k) = moo2 * a2

          enddo
       enddo
    enddo
end subroutine get_t_ini
!=============================================================================!
!=============================================================================!

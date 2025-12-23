!_____________________________________________________________________!
subroutine invcode( limiter )
    use global_variables,only : nscheme,nflux
    implicit none
    real,external :: limiter
    external  nnd_pv,WCNS_E_5,nnd_sq,zero_order,WCNS_E_5_J,WCNS_E_5_41,WCNS_E_5_42,WCNS_E_5_43, &
              WCNS_E_5_45,WCNS_E_5_47
    external  flux_SW,flux_Roe,flux_vanleer,flux_SW_mod,flux_RoeMx
!	external  flux_AUSMPW,flux_AUSM_plus,flux_AUSMPW_plus
    
	if ( nscheme == 1 ) then
       if ( nflux == 1 ) call inviscd3d( limiter,nnd_pv,flux_SW)
       !!if ( nflux == 1 ) call inviscd3d( limiter,nnd_pv,flux_SW_mod)
       if ( nflux == 2 ) call inviscd3d( limiter,nnd_pv,flux_vanleer)
       if ( nflux == 3 ) call inviscd3d( limiter,nnd_pv,flux_Roe)

    elseif ( nscheme == 4 ) then
			if(nflux == 1) call inviscd3d( limiter,WCNS_E_5,flux_SW)
			if(nflux == 2) call inviscd3d( limiter,WCNS_E_5,flux_vanleer)
			if(nflux == 3) call inviscd3d( limiter,WCNS_E_5,flux_Roe)
!    elseif ( nscheme == 40 ) then
!			if(nflux == 1) call inviscd3d( limiter,WCNS_E_5_J,flux_SW)
!			if(nflux == 2) call inviscd3d( limiter,WCNS_E_5_J,flux_vanleer)
!			if(nflux == 3) call inviscd3d( limiter,WCNS_E_5_J,flux_Roe)
    elseif ( nscheme == 41 ) then
			if(nflux == 1) call inviscd3d( limiter,WCNS_E_5_41,flux_SW)
			if(nflux == 2) call inviscd3d( limiter,WCNS_E_5_41,flux_vanleer)
			if(nflux == 3) call inviscd3d( limiter,WCNS_E_5_41,flux_Roe)
    elseif ( nscheme == 42 ) then
			if(nflux == 1) call inviscd3d( limiter,WCNS_E_5_42,flux_SW)
			!!if(nflux == 1) call inviscd3d( limiter,WCNS_E_5_42,flux_SW_mod)
			if(nflux == 2) call inviscd3d( limiter,WCNS_E_5_42,flux_vanleer)
			if(nflux == 3) call inviscd3d( limiter,WCNS_E_5_42,flux_Roe)
			if(nflux == 4) call inviscd3d( limiter,WCNS_E_5_42,flux_RoeMx)
!    elseif ( nscheme == 43 ) then
!			if(nflux == 1) call inviscd3d( limiter,WCNS_E_5_43,flux_SW)
!			if(nflux == 2) call inviscd3d( limiter,WCNS_E_5_43,flux_vanleer)
!			if(nflux == 3) call inviscd3d( limiter,WCNS_E_5_43,flux_Roe)
    elseif ( nscheme == 45 ) then
			if(nflux == 1) call inviscd3d( limiter,WCNS_E_5_45,flux_SW)
			if(nflux == 2) call inviscd3d( limiter,WCNS_E_5_45,flux_vanleer)
			if(nflux == 3) call inviscd3d( limiter,WCNS_E_5_45,flux_Roe)
    elseif ( nscheme == 47 ) then
			if(nflux == 1) call inviscd3d( limiter,WCNS_E_5_47,flux_SW)
			if(nflux == 2) call inviscd3d( limiter,WCNS_E_5_47,flux_vanleer)
			if(nflux == 3) call inviscd3d( limiter,WCNS_E_5_47,flux_Roe)
			if(nflux == 4) call inviscd3d( limiter,WCNS_E_5_47,flux_RoeMx)
    else
		write(*,*)'在选择计算格式时出错'
		write(*,*)'建立选1, 4，41,42，或47'
		stop

    endif
   
   return
end subroutine invcode
!_____________________________________________________________________!
subroutine inviscd3d( limiter,flux_line,flux_type )
!*TGH. 若某个方向的网格总数小于等于6，则不计算该方向（二维）
!*TGH. Modified by TU Guohua, 2009.2
!*TGH. 注意 观察修改后的正确性

    use global_variables
#ifdef PARALLEL
    use mod_parallels,only : mb_pids
#endif
    implicit none
    integer :: i,j,k,m,cmethd
    real,external :: limiter
    external  flux_line,flux_type
    real :: fc(1:nl,1:nmax),q_line(1:nl+1,-2:nmax+3,2),trxyz(5,nmax) !*TGH. NL=NM=5 *!
    real :: trxyz1(5,-3:nmax+4)   !!lhy
    real :: pr,pl   !!lhy
    real,pointer :: rpres(:,:,:,:) !!lhy
    integer :: nb,nr,idir,inrout,nbt
    integer :: i1,it0,jt0,kt0,it,jt,kt
    
    if (nflux == 4) then
       allocate(rpres(-1:ni+1,-1:nj+1,-1:nk+1,3))
       allocate(rpmin(0:nmax))
       
       do k=0,nk+1
       do j=0,nj+1
          do i=-1,ni+1
             pl = p(i,j,k)
             pr = p(i+1,j,k)
             rpres(i,j,k,1) = min(pl/pr,pr/pl)
          end do
 	      if (r(-2,j,k) < small) then
 	         rpres(0,j,k,1) = rpres(1,j,k,1)
 	         rpres(-1,j,k,1) = rpres(1,j,k,1)
	      end if

 	      if (r(ni+3,j,k) < small) then
 	         rpres(ni,j,k,1) = rpres(ni-1,j,k,1)
 	         rpres(ni+1,j,k,1) = rpres(ni-1,j,k,1)
	      end if
       end do
       end do       
       
       do k=0,nk+1
       do i=0,ni+1
          do j=-1,nj+1
             pl = p(i,j,k)
             pr = p(i,j+1,k)
             rpres(i,j,k,2) = min(pl/pr,pr/pl)
          end do
 	      if (r(i,-2,k) < small) then
 	         rpres(i,0,k,2) = rpres(i,1,k,2)
 	         rpres(i,-1,k,2) = rpres(i,1,k,2)
	      end if

 	      if (r(i,nj+3,k) < small) then
 	         rpres(i,nj,k,2) = rpres(i,nj-1,k,2)
 	         rpres(i,nj+1,k,2) = rpres(i,nj-1,k,2)
	      end if
       end do
       end do       
       
       do j=0,nj+1
       do i=0,ni+1
          do k=-1,nk+1
             pl = p(i,j,k)
             pr = p(i,j,k+1)
             rpres(i,j,k,3) = min(pl/pr,pr/pl)
          end do
 	      if (r(i,j,-2) < small) then
 	         rpres(i,j,0,3) = rpres(i,j,1,3)
 	         rpres(i,j,-1,3) = rpres(i,j,1,3)
	      end if

 	      if (r(i,j,nk+3) < small) then
 	         rpres(i,j,nk,3) = rpres(i,j,nk-1,3)
 	         rpres(i,j,nk+1,3) = rpres(i,j,nk-1,3)
	      end if
       end do
       end do         
    end if
    
    if ( nchem == 0 .and. nchem_source == 0 ) then
       do i=-2,nmax+3
          do m=1,2
             q_line(nl+1,i,m) = gama
          enddo
       enddo
    endif

    cmethd = 1 + method

     if(ni <= nijk2d)goto 11 !*tgh. 不计算该方向，二维
     do k= 1,nk         !___cic
       do j= 1,nj      !___cic

          do i= -2,ni+3    !1,ni+1
             q_line(1,i,1) = r(i,j,k)
             q_line(2,i,1) = u(i,j,k)
             q_line(3,i,1) = v(i,j,k)
             q_line(4,i,1) = w(i,j,k)
             q_line(5,i,1) = p(i,j,k)
             do m=6,nl
                q_line(m,i,1) = fs(m-5,i,j,k) !* TGH. 化学反应组份 ??? *!
             enddo
             do m=1,nl
                q_line(m,i,2) = q(m,i,j,k)
             enddo
             q_line(nl+1,i,2) = t(i,j,k)
          enddo
          
          if (nflux == 4) then
             do i=0,ni
                rpmin(i) = rpres(i,j,k,1)
                rpmin(i) = min(rpres(i,j,k,2),rpres(i  ,j+1,k,2), &
                               rpres(i-1,j,k,2),rpres(i-1,j+1,k,2),rpmin(i))
                rpmin(i) = min(rpres(i  ,j,k,3),rpres(i  ,j,k+1,3), &
                               rpres(i-1,j,k,3),rpres(i-1,j,k+1,3),rpmin(i))
             end do
          end if
          
          select case(nscheme)
          case(45,47)
             do i=1,ni
                trxyz1(1,i) = kcx(i,j,k)
                trxyz1(2,i) = kcy(i,j,k)
                trxyz1(3,i) = kcz(i,j,k)
                trxyz1(4,i) = kct(i,j,k)
                trxyz1(5,i) = vol(i,j,k)
             enddo
             
             nb = nbself
             
             nr = mb_flg(nb,1)%a3d(1,j,k)
             if (mb_bc(nb)%bc(nr)%bctype < 0) then
                nbt    = mb_bc(nb)%bc(nr)%nbt
                idir   = mb_bc(nb)%bc(nr)%t_nd
                inrout = mb_bc(nb)%bc(nr)%t_lr
                
    			it0 = mb_bc(nb)%bc(nr)%image(1,j,k)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(1,j,k)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(1,j,k)
    			
                m = 3*(idir-1) + 1
                do i=-3,0
                   i1 = 1 - i
                   it = it0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(1)
                   jt = jt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(2)
                   kt = kt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(3)
#ifdef PARALLEL                   
                   if (mb_pids(nb) /= mb_pids(nbt)) then              
                      trxyz1(1,i) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m) 
                      trxyz1(2,i) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+1) 
                      trxyz1(3,i) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+2) 
                      trxyz1(4,i) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,idir+9) 
                   else
#endif             
                      if (idir == 1) then
                         trxyz1(1,i) = inrout*mb_kcx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,i) = inrout*mb_kcy(nbt)%a3d(it,jt,kt)
                         trxyz1(3,i) = inrout*mb_kcz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,i) = inrout*mb_kct(nbt)%a3d(it,jt,kt)
                      else if (idir == 2) then
                         trxyz1(1,i) = inrout*mb_etx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,i) = inrout*mb_ety(nbt)%a3d(it,jt,kt)
                         trxyz1(3,i) = inrout*mb_etz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,i) = inrout*mb_ett(nbt)%a3d(it,jt,kt)
                      else
                         trxyz1(1,i) = inrout*mb_ctx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,i) = inrout*mb_cty(nbt)%a3d(it,jt,kt)
                         trxyz1(3,i) = inrout*mb_ctz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,i) = inrout*mb_ctt(nbt)%a3d(it,jt,kt)
                      end if
#ifdef PARALLEL                   
                   end if
#endif                   
                end do
             else
                do i=-3,0
                   trxyz1(1:5,i) = trxyz1(1:5,1)
                end do
             end if
             
             nr = mb_flg(nb,2)%a3d(ni,j,k)
             if (mb_bc(nb)%bc(nr)%bctype < 0) then
                nbt    = mb_bc(nb)%bc(nr)%nbt
                idir   = mb_bc(nb)%bc(nr)%t_nd
                inrout = mb_bc(nb)%bc(nr)%t_lr
                
    			it0 = mb_bc(nb)%bc(nr)%image(ni,j,k)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(ni,j,k)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(ni,j,k)
    			
                m = 3*(idir-1) + 1
                do i=ni+1,ni+4
                   i1 = i - ni
                   it = it0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(1)
                   jt = jt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(2)
                   kt = kt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(3)
#ifdef PARALLEL                   
                   if (mb_pids(nb) /= mb_pids(nbt)) then              
                      trxyz1(1,i) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m) 
                      trxyz1(2,i) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+1) 
                      trxyz1(3,i) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+2) 
                      trxyz1(4,i) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,idir+9) 
                   else
#endif            
                      if (idir == 1) then
                         trxyz1(1,i) = -inrout*mb_kcx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,i) = -inrout*mb_kcy(nbt)%a3d(it,jt,kt)
                         trxyz1(3,i) = -inrout*mb_kcz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,i) = -inrout*mb_kct(nbt)%a3d(it,jt,kt)
                      else if (idir == 2) then
                         trxyz1(1,i) = -inrout*mb_etx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,i) = -inrout*mb_ety(nbt)%a3d(it,jt,kt)
                         trxyz1(3,i) = -inrout*mb_etz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,i) = -inrout*mb_ett(nbt)%a3d(it,jt,kt)
                      else
                         trxyz1(1,i) = -inrout*mb_ctx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,i) = -inrout*mb_cty(nbt)%a3d(it,jt,kt)
                         trxyz1(3,i) = -inrout*mb_ctz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,i) = -inrout*mb_ctt(nbt)%a3d(it,jt,kt)
                      end if
#ifdef PARALLEL                   
                   end if
#endif    
                end do
             else
                do i=ni+1,ni+4
                   trxyz1(1:5,i) = trxyz1(1:5,ni)
                end do
             end if
                          
             call flux_line(xk,xb,ni,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz1,q_line,fc) 
          case default
             do i=1,ni
                trxyz(1,i) = kcx(i,j,k)
                trxyz(2,i) = kcy(i,j,k)
                trxyz(3,i) = kcz(i,j,k)
                trxyz(4,i) = kct(i,j,k)
                trxyz(5,i) = vol(i,j,k)
             enddo
             call flux_line(xk,xb,ni,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz,q_line,fc) 
          end select

          do i=1,ni     !___cic
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) +  fc(m,i)
             enddo
          enddo

       enddo
    enddo

11  continue

    if(nj <= nijk2d)goto 12 !*tgh. 不计算该方向，二维
    do k=1,nk         !___cic
       do i=1,ni      !___cic


          do j=-2,nj+3  !1,nj+1
             q_line(1,j,1) = r(i,j,k)
             q_line(2,j,1) = u(i,j,k)
             q_line(3,j,1) = v(i,j,k)
             q_line(4,j,1) = w(i,j,k)
             q_line(5,j,1) = p(i,j,k)
             do m=6,nl
                q_line(m,j,1) = fs(m-5,i,j,k)
             enddo
             do m=1,nl
                q_line(m,j,2) = q(m,i,j,k)
             enddo

             q_line(nl+1,j,2) = t(i,j,k)
          enddo
          
          if (nflux == 4) then
             do j=0,nj
                rpmin(j) = rpres(i,j,k,2)
                rpmin(j) = min(rpres(i,j  ,k,1),rpres(i+1,j  ,k,1), &
                               rpres(i,j-1,k,1),rpres(i+1,j-1,k,1),rpmin(j))
                rpmin(j) = min(rpres(i,j  ,k,3),rpres(i,j  ,k+1,3), &
                               rpres(i,j-1,k,3),rpres(i,j-1,k+1,3),rpmin(j))
             end do
          end if

          select case(nscheme)
          case(45,47)
             do j=1,nj
                trxyz1(1,j) = etx(i,j,k)
                trxyz1(2,j) = ety(i,j,k)
                trxyz1(3,j) = etz(i,j,k)
                trxyz1(4,j) = ett(i,j,k)
                trxyz1(5,j) = vol(i,j,k)
             end do
             
             nb = nbself
             
             nr = mb_flg(nb,3)%a3d(i,1,k)
             if (mb_bc(nb)%bc(nr)%bctype < 0) then
                nbt    = mb_bc(nb)%bc(nr)%nbt
                idir   = mb_bc(nb)%bc(nr)%t_nd
                inrout = mb_bc(nb)%bc(nr)%t_lr
                
    			it0 = mb_bc(nb)%bc(nr)%image(i,1,k)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(i,1,k)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(i,1,k)
    			
                m = 3*(idir-1) + 1
                do j=-3,0
                   i1 = 1 - j
                   it = it0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(1)
                   jt = jt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(2)
                   kt = kt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(3)
#ifdef PARALLEL                   
                   if (mb_pids(nb) /= mb_pids(nbt)) then              
                      trxyz1(1,j) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m) 
                      trxyz1(2,j) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+1) 
                      trxyz1(3,j) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+2) 
                      trxyz1(4,j) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,idir+9)
                   else
#endif    
                      if (idir == 1) then
                         trxyz1(1,j) = inrout*mb_kcx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,j) = inrout*mb_kcy(nbt)%a3d(it,jt,kt)
                         trxyz1(3,j) = inrout*mb_kcz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,j) = inrout*mb_kct(nbt)%a3d(it,jt,kt)
                      else if (idir == 2) then
                         trxyz1(1,j) = inrout*mb_etx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,j) = inrout*mb_ety(nbt)%a3d(it,jt,kt)
                         trxyz1(3,j) = inrout*mb_etz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,j) = inrout*mb_ett(nbt)%a3d(it,jt,kt)
                      else
                         trxyz1(1,j) = inrout*mb_ctx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,j) = inrout*mb_cty(nbt)%a3d(it,jt,kt)
                         trxyz1(3,j) = inrout*mb_ctz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,j) = inrout*mb_ctt(nbt)%a3d(it,jt,kt)
                      end if
#ifdef PARALLEL                   
                   end if
#endif    
                end do
             else
                do j=-3,0
                   trxyz1(1:5,j) = trxyz1(1:5,1)
                end do
             end if
             
             nr = mb_flg(nb,4)%a3d(i,nj,k)
             if (mb_bc(nb)%bc(nr)%bctype < 0) then
                nbt    = mb_bc(nb)%bc(nr)%nbt
                idir   = mb_bc(nb)%bc(nr)%t_nd
                inrout = mb_bc(nb)%bc(nr)%t_lr
                
    			it0 = mb_bc(nb)%bc(nr)%image(i,nj,k)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(i,nj,k)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(i,nj,k)
    			
                m = 3*(idir-1) + 1
                do j=nj+1,nj+4
                   i1 = j - nj
                   it = it0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(1)
                   jt = jt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(2)
                   kt = kt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(3)
#ifdef PARALLEL                   
                   if (mb_pids(nb) /= mb_pids(nbt)) then              
                      trxyz1(1,j) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m) 
                      trxyz1(2,j) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+1) 
                      trxyz1(3,j) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+2) 
                      trxyz1(4,j) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,idir+9)
                   else
#endif             
                      if (idir == 1) then
                         trxyz1(1,j) = -inrout*mb_kcx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,j) = -inrout*mb_kcy(nbt)%a3d(it,jt,kt)
                         trxyz1(3,j) = -inrout*mb_kcz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,j) = -inrout*mb_kct(nbt)%a3d(it,jt,kt)
                      else if (idir == 2) then
                         trxyz1(1,j) = -inrout*mb_etx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,j) = -inrout*mb_ety(nbt)%a3d(it,jt,kt)
                         trxyz1(3,j) = -inrout*mb_etz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,j) = -inrout*mb_ett(nbt)%a3d(it,jt,kt)
                      else
                         trxyz1(1,j) = -inrout*mb_ctx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,j) = -inrout*mb_cty(nbt)%a3d(it,jt,kt)
                         trxyz1(3,j) = -inrout*mb_ctz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,j) = -inrout*mb_ctt(nbt)%a3d(it,jt,kt)
                      end if
#ifdef PARALLEL                   
                   end if
#endif    
                end do
             else
                do j=nj+1,nj+4
                   trxyz1(1:5,j) = trxyz1(1:5,nj)
                end do
             end if
             call flux_line(xk,xb,nj,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz1,q_line,fc)
          case default
             do j=1,nj
                trxyz(1,j) = etx(i,j,k)
                trxyz(2,j) = ety(i,j,k)
                trxyz(3,j) = etz(i,j,k)
                trxyz(4,j) = ett(i,j,k)
                trxyz(5,j) = vol(i,j,k)
             end do
             call flux_line(xk,xb,nj,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz,q_line,fc)
          end select

          do j=1,nj     !___cic
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) +  fc(m,j)
             enddo
          enddo


       enddo
    enddo
12  continue

    if(nk <= nijk2d)goto 13 !*tgh. 不计算该方向，二维
    do j=1,nj       !___cic
       do i=1,ni    !___cic


          do k= -2,nk+3    !1,nk+1
             q_line(1,k,1) = r(i,j,k)
             q_line(2,k,1) = u(i,j,k)
             q_line(3,k,1) = v(i,j,k)
             q_line(4,k,1) = w(i,j,k)
             q_line(5,k,1) = p(i,j,k)
             do m=6,nl
                q_line(m,k,1) = fs(m-5,i,j,k)
             enddo
             do m=1,nl
                q_line(m,k,2) = q(m,i,j,k)
             enddo
             q_line(nl+1,k,2) = t(i,j,k)
          enddo

          if (nflux == 4) then
             do k=0,nk
                rpmin(k) = rpres(i,j,k,3)
                rpmin(k) = min(rpres(i,j,k  ,1),rpres(i+1,j,k  ,1), &
                               rpres(i,j,k-1,1),rpres(i+1,j,k-1,1),rpmin(k))
                rpmin(k) = min(rpres(i,j,k  ,2),rpres(i,j+1,k  ,2), &
                               rpres(i,j,k-1,2),rpres(i,j+1,k-1,2),rpmin(k))
             end do
          end if

          
          select case(nscheme)
          case(45,47)
             do k=1,nk
                trxyz1(1,k) = ctx(i,j,k)
                trxyz1(2,k) = cty(i,j,k)
                trxyz1(3,k) = ctz(i,j,k)
                trxyz1(4,k) = ctt(i,j,k)
                trxyz1(5,k) = vol(i,j,k)
             enddo
             
             nb = nbself
             
             nr = mb_flg(nb,5)%a3d(i,j,1)
             if (mb_bc(nb)%bc(nr)%bctype < 0) then
                nbt    = mb_bc(nb)%bc(nr)%nbt
                idir   = mb_bc(nb)%bc(nr)%t_nd
                inrout = mb_bc(nb)%bc(nr)%t_lr
                
    			it0 = mb_bc(nb)%bc(nr)%image(i,j,1)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,1)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,1)
    			
                m = 3*(idir-1) + 1
                do k=-3,0
                   i1 = 1 - k
                   it = it0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(1)
                   jt = jt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(2)
                   kt = kt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(3)
#ifdef PARALLEL    
                   if (mb_pids(nb) /= mb_pids(nbt)) then              
                      trxyz1(1,k) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m) 
                      trxyz1(2,k) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+1) 
                      trxyz1(3,k) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+2) 
                      trxyz1(4,k) = inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,idir+9)
                   else
#endif                   
                      if (idir == 1) then
                         trxyz1(1,k) = inrout*mb_kcx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,k) = inrout*mb_kcy(nbt)%a3d(it,jt,kt)
                         trxyz1(3,k) = inrout*mb_kcz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,k) = inrout*mb_kct(nbt)%a3d(it,jt,kt)
                      else if (idir == 2) then
                         trxyz1(1,k) = inrout*mb_etx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,k) = inrout*mb_ety(nbt)%a3d(it,jt,kt)
                         trxyz1(3,k) = inrout*mb_etz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,k) = inrout*mb_ett(nbt)%a3d(it,jt,kt)
                      else
                         trxyz1(1,k) = inrout*mb_ctx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,k) = inrout*mb_cty(nbt)%a3d(it,jt,kt)
                         trxyz1(3,k) = inrout*mb_ctz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,k) = inrout*mb_ctt(nbt)%a3d(it,jt,kt)
                      end if
#ifdef PARALLEL    
                  end if
#endif                   
                end do
             else
                do k=-3,0
                   trxyz1(1:5,k) = trxyz1(1:5,1)
                end do
             end if
             
             nr = mb_flg(nb,6)%a3d(i,j,nk)
             if (mb_bc(nb)%bc(nr)%bctype < 0) then
                nbt    = mb_bc(nb)%bc(nr)%nbt
                idir   = mb_bc(nb)%bc(nr)%t_nd
                inrout = mb_bc(nb)%bc(nr)%t_lr
                
    			it0 = mb_bc(nb)%bc(nr)%image(i,j,nk)
    			jt0 = mb_bc(nb)%bc(nr)%jmage(i,j,nk)
    			kt0 = mb_bc(nb)%bc(nr)%kmage(i,j,nk)
    			
                m = 3*(idir-1) + 1
                do k=nk+1,nk+4
                   i1 = k - nk
                   it = it0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(1)
                   jt = jt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(2)
                   kt = kt0 - i1*mb_bc(nb)%bc(nr)%t_lr3d(3)
#ifdef PARALLEL                   
                   if (mb_pids(nb) /= mb_pids(nbt)) then              
                      trxyz1(1,k) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m) 
                      trxyz1(2,k) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+1) 
                      trxyz1(3,k) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,m+2) 
                      trxyz1(4,k) = -inrout*mb_bc(nb)%bc(nr)%sxyzpack(it,jt,kt,idir+9)
                   else
#endif                   
                      if (idir == 1) then
                         trxyz1(1,k) = -inrout*mb_kcx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,k) = -inrout*mb_kcy(nbt)%a3d(it,jt,kt)
                         trxyz1(3,k) = -inrout*mb_kcz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,k) = -inrout*mb_kct(nbt)%a3d(it,jt,kt)
                      else if (idir == 2) then
                         trxyz1(1,k) = -inrout*mb_etx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,k) = -inrout*mb_ety(nbt)%a3d(it,jt,kt)
                         trxyz1(3,k) = -inrout*mb_etz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,k) = -inrout*mb_ett(nbt)%a3d(it,jt,kt)
                      else
                         trxyz1(1,k) = -inrout*mb_ctx(nbt)%a3d(it,jt,kt) 
                         trxyz1(2,k) = -inrout*mb_cty(nbt)%a3d(it,jt,kt)
                         trxyz1(3,k) = -inrout*mb_ctz(nbt)%a3d(it,jt,kt)
                         trxyz1(4,k) = -inrout*mb_ctt(nbt)%a3d(it,jt,kt)
                      end if
#ifdef PARALLEL                   
                   end if
#endif                   
                end do
             else
                do k=nk+1,nk+4
                   trxyz1(1:5,k) = trxyz1(1:5,nk)
                end do
             end if
             call flux_line(xk,xb,nk,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz1,q_line,fc)
          case default          
             do k=1,nk
                trxyz(1,k) = ctx(i,j,k)
                trxyz(2,k) = cty(i,j,k)
                trxyz(3,k) = ctz(i,j,k)
                trxyz(4,k) = ctt(i,j,k)
                trxyz(5,k) = vol(i,j,k)
             enddo
             call flux_line(xk,xb,nk,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz,q_line,fc)
          end select

          do k=1,nk      !___cic
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) + fc(m,k)
             enddo
          enddo


       enddo
    enddo
13  continue

    if (nflux == 4) then
       deallocate(rpmin,rpres)
    end if

    return
end subroutine inviscd3d
!_____________________________________________________________________!
subroutine zero_order(xk,xb,ni,nmax,nl,method,cmethd,flux_type, &
                           limiter,efix,trxyz,q_line,dfc)
    implicit none
    integer :: i,m,cmethd,ni,nmax,nl,method
    real,external :: limiter
    external  flux_type
    real,dimension( 1:nl ) :: dql,dqr,f
    real :: nx,ny,nz,nt,xk,xb,efix,gamaeql,gamaeqr
    real :: c1,c2,dfc(1:nl,1:nmax)
    real :: fc(1:nl,1:nmax),q_line(1:nl+1,-1:nmax+1,2),trxyz(5,nmax)
    real :: d_q(1:nl+1,0:nmax+1),d_qcl(1:nl+1,0:nmax),d_qcr(1:nl+1,0:nmax)

!----------------------------------------------------------
!* Modified by TU Guohua
	if(ni.le.4)then
      do i=1,ni
       do m=1,nl
          dfc(m,i) = 0.0
       enddo
      enddo
      return
	endif
!* end. Modified by TU Guohua, 2009.2
!----------------------------------------------------------

    do i=cmethd,ni      !2,ni
       nx = trxyz(1,i)
       ny = trxyz(2,i)
       nz = trxyz(3,i)
       nt = trxyz(4,i)
       do m=1,method
          nx = 0.50 * ( trxyz(1,i-1) + nx )
          ny = 0.50 * ( trxyz(2,i-1) + ny )
          nz = 0.50 * ( trxyz(3,i-1) + nz )
          nt = 0.50 * ( trxyz(4,i-1) + nt )
       enddo
       do m=1,nl
          dql(m) = q_line(m,i-1,1)
          dqr(m) = q_line(m,i  ,1)
       enddo

       m = nl + 1
       gamaeql = q_line(m,i-1,1)
       gamaeqr = q_line(m,i  ,1) 

       call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)

       do m=1,nl
          fc(m,i) = f(m)
       enddo
    enddo

    do i=cmethd,ni-1     !2,ni-1
       do m=1,nl
          dfc(m,i) =  fc(m,i+1) - fc(m,i)
       enddo
    enddo
    return
end subroutine zero_order  
!_____________________________________________________________________!
subroutine center(xk,xb,ni,nmax,nl,method,cmeth,flux_type, &
                          limiter,efix,trxyz,q_line,dfc)
    implicit none
    integer :: i,j,k,m,cmeth,cmeth1,ni,nmax,nl,method,n
    real,external :: limiter,flux_type
    real :: dql(nl),dqr(nl),fl(nl),fr(nl)
    real :: fc(1:nl,1:nmax),q_line(1:nl+1,-1:nmax+1,2),trxyz(5,nmax)
    real :: xk,efix,xb,nx,ny,nz,nt
    real :: f12(nl),prim(nl),em,dfc(1:nl,1:nmax),gama
    
    cmeth1= cmeth-1

    do i=1,ni  
       do m=1,nl
          prim(m)=  q_line(m,i,1) 
        enddo
				gama = q_line(nl+1,i,1)
        nx = trxyz(1,i)
        ny = trxyz(2,i)
        nz = trxyz(3,i)
        nt = trxyz(4,i)
        call flux(prim,gama,nl,nx,ny,nz,nt,f12)  !有限体积为右边通量，差分为格点处通量
        do m=1,nl
           fc(m,i) = f12(m)
        enddo
        do n=1,1-method
           do m=1,nl
              prim(m)=  q_line(m,i-1,1) 
           enddo
					 gama = q_line(nl+1,i-1,1)
           call flux(prim,gama,nl,nx,ny,nz,nt,f12)  !有限体积为左边通量
           do m=1,nl
              fc(m,i) = fc(m,i) + f12(m)
           enddo
        enddo
     enddo

     do i=cmeth,ni-1                         !2,ni-1     !有限差分METHOD＝1，有限体积METHOD＝0
        do m=1,nl
           dfc(m,i) = 0.5*(fc(m,i+1) - fc(m,i-method))    !有限差分为单元左右格点处通量之差
        enddo
     enddo

    return
end subroutine center
!_____________________________________________________________________!
subroutine nnd_sq(xk,xb,ni,nmax,nl,method,cmeth,flux_type, &
                          limiter,efix,trxyz,q_line,dfc)
    implicit none
    integer :: i,j,k,m,cmeth,cmeth1,ni,nmax,nl,method
    real,external :: limiter
    external  flux_type
    real :: dql(nl),dqr(nl),fl(nl),fr(nl),prim(nl)
    real :: gykb(4),dfc(1:nl,1:nmax),flr(nl),beta,bet1,bet2
    real :: fc(1:nl,1:nmax),q_line(1:nl+1,-1:nmax+1,2),trxyz(5,nmax)
    real :: d_q(1:nl,0:nmax+1),d_qc(1:nl,0:nmax),xk,efix,xb

    call center(xk,xb,ni,nmax,nl,method,cmeth,flux_type, &
                         limiter,efix,trxyz,q_line,dfc)
    cmeth1= cmeth-1
    do i=method,ni+1 
       do m=1,nl
          d_q(m,i) = q_line(m,i,2) - q_line(m,i-1,2)
       enddo
    enddo
    
    do j=1,method
       do m=1,nl
          d_q(m,cmeth1) = d_q(m,cmeth)
          d_q(m,ni+1  ) = d_q(m,ni)
       enddo
    enddo
    do i=cmeth1,ni
       do m=1,nl
          d_qc(m,i) = limiter( d_q(m,i+1), d_q(m,i) )
       enddo
    enddo

    do i=cmeth,ni
       do j=1,4
          gykb(j) = trxyz(j,i)
          do m=1,method
             gykb(j) = 0.50 * ( trxyz(j,i-1) + gykb(j) )
          enddo
       enddo

       do m=1,nl
          prim(m)= 0.50 * ( q_line(m,i-1,1) + q_line(m,i,1) )
          dql(m) = d_q(m,i) - d_qc(m,i-1)
          dqr(m) = d_q(m,i) - d_qc(m,i  )
       enddo

       call mxdq1(prim,gykb,dql,dqr,fl,fr,efix)

       do m=1,nl
          fc(m,i) = - 0.5 * ( fl(m) - fr(m) )
       enddo
    enddo

    do i=cmeth,ni-1     !2,ni-1
       do m=1,nl
          dfc(m,i) = dfc(m,i) + fc(m,i+1) - fc(m,i)
       enddo
    enddo

    return
end subroutine nnd_sq   !
!_____________________________________________________________________!
subroutine nnd_pv(xk,xb,ni,nmax,nl,method,cmethd,flux_type, &
                           limiter,efix,trxyz,q_line,dfc)
!* TGH. xk=-1.0; xb=1.0  *!
!----------------------------------------------------------------------
!*
!* TGH. 把求得的导数从2~MAX-1 扩展到1~max
!* TGH. 在边界上用一阶单侧差分
!* TGH. 原结构在第2点和第max-1点要用到虚点上的值
!* TGH. 修改后不用虚点的值
!* TGH. Euler边界：边界点1阶单侧，边界第二点1阶迎风
!* TGH.            Euler时边界附近的压力等值线比较难看，为了提高精度可以修改第二点的插值格式
!* TGH. N-S边界：  边界点1阶单侧，边界第2~m点2阶中心或高阶偏心，需要根据情况调节m和精度
!----------------------------------------------------------------------
!* TGH. 还待完善
!----------------------------------------------------------------------
	use define_precision_mod
	use global_variables,only : nlimiter,nvis
    implicit none
    integer :: i,m,cmethd,ni,nmax,nl,method
    real,external :: limiter
    external  flux_type
    real,dimension( 1:nl ) :: dql,dqr,f
    real :: nx,ny,nz,nt,xk,xb,efix,gamaeql,gamaeqr
    real :: c1,c2,dfc(1:nl,1:nmax)
    real :: fc(1:nl,1:nmax+1),q_line(1:nl+1,-2:nmax+3,2),trxyz(5,nmax)
    real :: d_q(1:nl+1,0:nmax+2),d_qcl(1:nl+1,0:nmax+1),d_qcr(1:nl+1,0:nmax+1)
    logical :: st_log,ed_log,euler_log
    
    euler_log = .false.
    st_log    = .false.
    ed_log    = .false.

     if(nvis==0)           euler_log = .true.
     if( (nlimiter /= 0) .and. (q_line(1,-2,1)  <=0.) ) st_log = .true.
     if( (nlimiter /= 0) .and. (q_line(1,ni+3,1)<=0.) ) ed_log = .true.

    c1 = 0.25 * ( 1.0 + xk )   !*tgh. c1=0.0
    c2 = 0.25 * ( 1.0 - xk )   !*tgh. c2=0.5
    do i=0,ni+2  !* 1~max+1
       do m=1,nl
          d_q(m,i) = q_line(m,i,1) - q_line(m,i-1,1)
       enddo
    enddo
 
     do i=0,ni+1     !* 0~max+1 
       do m=1,nl+1
          d_qcl(m,i) = limiter( d_q(m,i+1), d_q(m,i  ) ) !* 用到的d_q：0~max+2
          d_qcr(m,i) = limiter( d_q(m,i  ), d_q(m,i+1) )
       enddo
    enddo

    do i=1,ni+1 !*TGH. 1~max+1， 用到的q_line为1~max+1
	   if(i==1)then
          nx = 1.5_prec * trxyz(1,i) - 0.5_prec*trxyz(1,2)
          ny = 1.5_prec * trxyz(2,i) - 0.5_prec*trxyz(2,2)
          nz = 1.5_prec * trxyz(3,i) - 0.5_prec*trxyz(3,2)
          nt = 1.5_prec * trxyz(4,i) - 0.5_prec*trxyz(4,2)
	   elseif(i==ni+1)then
          nx = 1.5_prec * trxyz(1,ni) - 0.5_prec*trxyz(1,ni-1)
          ny = 1.5_prec * trxyz(2,ni) - 0.5_prec*trxyz(2,ni-1)
          nz = 1.5_prec * trxyz(3,ni) - 0.5_prec*trxyz(3,ni-1)
          nt = 1.5_prec * trxyz(4,ni) - 0.5_prec*trxyz(4,ni-1)
	   else
          nx = 0.5_prec * (trxyz(1,i) + trxyz(1,i-1) )
          ny = 0.5_prec * (trxyz(2,i) + trxyz(2,i-1) )
          nz = 0.5_prec * (trxyz(3,i) + trxyz(3,i-1) )
          nt = 0.5_prec * (trxyz(4,i) + trxyz(4,i-1) )
	   endif

       do m=1,nl
          dql(m) = q_line(m,i-1,1) + c1 * d_qcr(m,i-1) + c2 * d_qcl(m,i-1) !* 用到的d_qcl：0~max
          dqr(m) = q_line(m,i  ,1) - c1 * d_qcl(m,i  ) - c2 * d_qcr(m,i  ) !* 用到的d_qcl：1~max+1
       enddo
       
       m = nl + 1
       gamaeql = q_line(m,i-1,1)
       gamaeqr = gamaeql

       call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)

       do m=1,nl
          fc(m,i) = f(m) !*TGH. 1~max+1 对应半节点1/2 ~ max+1/2
       enddo
    enddo
    
!--------------------------------------!
! 物面边界上不用虚点
!!-----------------------
!! 补充边界格式
!! 补充边界格式
    if(st_log)then !左边界
      if(euler_log)then
        do i=2,3
           nx = 0.5 * (trxyz(1,i) + trxyz(1,i-1) )
           ny = 0.5 * (trxyz(2,i) + trxyz(2,i-1) )
           nz = 0.5 * (trxyz(3,i) + trxyz(3,i-1) )
           nt = 0.5 * (trxyz(4,i) + trxyz(4,i-1) )
           do m=1,nl       
             dql(m) =  q_line(m,i-1,1)
             dqr(m) =  q_line(m,i  ,1)
           enddo
           call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)
           do m=1,nl
              fc(m,i) = f(m)
           enddo
        enddo
      else
		i=2 !对3/2点采用2阶中心和3阶迎风偏心插值混合
        nx = 0.5 * (trxyz(1,i) + trxyz(1,i-1) )
        ny = 0.5 * (trxyz(2,i) + trxyz(2,i-1) )
        nz = 0.5 * (trxyz(3,i) + trxyz(3,i-1) )
        nt = 0.5 * (trxyz(4,i) + trxyz(4,i-1) )
        do m=1,nl
          dql(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
!          dqr(m) = ( 3.0*q_line(m,i-1,1)+6.0*q_line(m,i  ,1) - q_line(m,i+1,1) )/8.0
          dqr(m) = ( 2.0*q_line(m,i-1,1)+5.0*q_line(m,i  ,1) - q_line(m,i+1,1) )/6.0
        enddo
        if(dqr(1) < 1.e-15 .or. dqr(5)<1.e-15)then
		   do m=1,nl
              dqr(m) = dql(m)
		   enddo
        endif
        call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)
        do m=1,nl
          fc(m,i) = f(m)
        enddo

		do i=3,5  !对5/2,7/2,8/2点采用3阶迎风偏心插值
           nx = 0.5 * (trxyz(1,i) + trxyz(1,i-1) )
           ny = 0.5 * (trxyz(2,i) + trxyz(2,i-1) )
           nz = 0.5 * (trxyz(3,i) + trxyz(3,i-1) )
           nt = 0.5 * (trxyz(4,i) + trxyz(4,i-1) )
           do m=1,nl
!             dql(m) = ( 9.0*(q_line(m,i  ,1)+q_line(m,i-1,1)) - (q_line(m,i+1,1) + q_line(m,i-2,1)) )/16.0
!			 dqr(m) = dql(m)

!             dql(m) = ( 3.0*q_line(m,i  ,1)+6.0*q_line(m,i-1,1) - q_line(m,i-2,1) )/8.0
!             dqr(m) = ( 3.0*q_line(m,i-1,1)+6.0*q_line(m,i  ,1) - q_line(m,i+1,1) )/8.0

!              dql(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
             dql(m) = ( 2.0*q_line(m,i  ,1)+5.0*q_line(m,i-1,1) - q_line(m,i-2,1) )/6.0
             dqr(m) = ( 2.0*q_line(m,i-1,1)+5.0*q_line(m,i  ,1) - q_line(m,i+1,1) )/6.0
           enddo

           if(dql(1) < 1.e-15 .or. dql(5)<1.e-15)then
               do m=1,nl
                 dql(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
			   enddo
           endif
           if(dqr(1) < 1.e-15 .or. dqr(5)<1.e-15)then
               do m=1,nl
                 dqr(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
			   enddo
           endif

           call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)
           do m=1,nl
             fc(m,i) = f(m)
           enddo
		enddo
      endif
    endif

    if(ed_log)then !右边界
      if(euler_log)then
        do i=ni-1,ni
           nx = 0.5 * (trxyz(1,i) + trxyz(1,i-1) )
           ny = 0.5 * (trxyz(2,i) + trxyz(2,i-1) )
           nz = 0.5 * (trxyz(3,i) + trxyz(3,i-1) )
           nt = 0.5 * (trxyz(4,i) + trxyz(4,i-1) )
           do m=1,nl
             dql(m) =  q_line(m,i-1,1)
             dqr(m) =  q_line(m,i  ,1)
           enddo
           call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)
           do m=1,nl
             fc(m,i) = f(m)
           enddo
        enddo
      else
		 i = ni !对ni-1/2点采用2阶中心和3阶迎风偏心插值混合
         nx = 0.5 * (trxyz(1,i) + trxyz(1,i-1) )
         ny = 0.5 * (trxyz(2,i) + trxyz(2,i-1) )
         nz = 0.5 * (trxyz(3,i) + trxyz(3,i-1) )
         nt = 0.5 * (trxyz(4,i) + trxyz(4,i-1) )
         do m=1,nl
!           dql(m) = ( 3.0*q_line(m,i  ,1)+6.0*q_line(m,i-1,1) - q_line(m,i-2,1) )/8.0
           dql(m) = ( 2.0*q_line(m,i  ,1)+5.0*q_line(m,i-1,1) - q_line(m,i-2,1) )/6.0
           dqr(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
         enddo

         if(dql(1) < 1.e-15 .or. dql(5)<1.e-15)then
            do m=1,nl
               dql(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
			enddo
         endif

         call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)
         do m=1,nl
           fc(m,i) = f(m)
         enddo

        do i=ni-3,ni-1 !对ni-3/2,ni-5/2,ni-7/2点采用3阶迎风偏心插值
           nx = 0.5 * (trxyz(1,i) + trxyz(1,i-1) )
           ny = 0.5 * (trxyz(2,i) + trxyz(2,i-1) )
           nz = 0.5 * (trxyz(3,i) + trxyz(3,i-1) )
           nt = 0.5 * (trxyz(4,i) + trxyz(4,i-1) )
           do m=1,nl
!             dql(m) = ( 9.0*(q_line(m,i  ,1)+q_line(m,i-1,1)) - (q_line(m,i+1,1) + q_line(m,i-2,1)) )/16.0
!			 dqr(m) = dql(m)

!             dql(m) = ( 3.0*q_line(m,i  ,1)+6.0*q_line(m,i-1,1) - q_line(m,i-2,1) )/8.0
!             dqr(m) = ( 3.0*q_line(m,i-1,1)+6.0*q_line(m,i  ,1) - q_line(m,i+1,1) )/8.0

!              dqr(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
             dql(m) = ( 2.0*q_line(m,i  ,1)+5.0*q_line(m,i-1,1) - q_line(m,i-2,1) )/6.0
             dqr(m) = ( 2.0*q_line(m,i-1,1)+5.0*q_line(m,i  ,1) - q_line(m,i+1,1) )/6.0
           enddo
           if(dql(1) < 1.e-15 .or. dql(5)<1.e-15)then
              do m=1,nl
                 dql(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
			  enddo
           endif
           if(dqr(1) < 1.e-15 .or. dqr(5)<1.e-15)then
              do m=1,nl
                 dqr(m) = 0.5*(q_line(m,i,1)+q_line(m,i-1,1))
			  enddo
           endif

           call flux_type(dql,dqr,nl,nx,ny,nz,nt,f,efix,gamaeql,gamaeqr)
           do m=1,nl
             fc(m,i) = f(m)
           enddo
        enddo
      endif
  
    endif
!!----- <end> 补充边界格式--------
! <end>  物面边界上不用虚点  
!--------------------------------------!


    do i=method,ni     !* 1~max
       do m=1,nl
          dfc(m,i) =  fc(m,i+1) - fc(m,i)
       enddo
    enddo

    return
end subroutine nnd_pv  
!_____________________________________________________________________!
subroutine mxdq1(prim,gykb,dql,dqr,fl,fr,efix)
    use global_const,only:nl,ns,ms1,beta1,tref,nchem,sml_sss
    implicit none
    integer :: m
    real :: hint,gama,ae,af,efix,eps
    real :: prim(nl),gykb(4),dql(nl),dqr(nl),fl(nl),fr(nl)
    real :: hs(ns),as(ns)
    real :: nx,ny,nz,nt,ct,cgm,cgm1,st
    real :: x1l,x2l,x1r,x2r
    real :: l1,l4,l5,l1fix,l4fix,l5fix,l1p,l4p,l5p,l1n,l4n,l5n
    real :: dhl,dcl,c2dcl,dhr,dcr,c2dcr
    real :: rm,um,vm,wm,pm,cm,c2,v2,tm,hm

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    call gettmhmgm(prim,tm,hint,gama)
    c2 = gama*pm/rm
    cm = sqrt(c2)
    v2 = um*um + vm*vm + wm*wm

    hm = hint + 0.5*v2

    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)

    ct = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
    eps = efix*efix*cgm*cgm

    l1 = ct
    l4 = ct + cm * cgm
    l5 = ct - cm * cgm

    l1fix = sqrt( l1*l1 + eps )
    l4fix = sqrt( l4*l4 + eps )
    l5fix = sqrt( l5*l5 + eps )

    l1p = 0.5 * ( l1 + l1fix )
    l4p = 0.5 * ( l4 + l4fix )
    l5p = 0.5 * ( l5 + l5fix )

    l1n = 0.5 * ( l1 - l1fix )
    l4n = 0.5 * ( l4 - l4fix )
    l5n = 0.5 * ( l5 - l5fix )

    x1l = ( 2.0*l1p - l4p - l5p )/( 2.0 * c2 )
    x2l = ( l4p - l5p )/( 2.0 * cm )

    x1r = ( 2.0*l1n - l4n - l5n )/( 2.0 * c2 )
    x2r = ( l4n - l5n )/( 2.0 * cm )

    cgm1 = 1.0/cgm
    nx = nx * cgm1
    ny = ny * cgm1
    nz = nz * cgm1
    st = (ct - nt )* cgm1 ! 减去nt ，在st中不包含nt 
    ae = gama - 1.0
    af = 0.5*ae*v2

    dcl = st * dql(1) - nx * dql(2) - ny * dql(3) - nz * dql(4)
    dhl = af * dql(1) - ae * ( um * dql(2) + vm * dql(3) + wm * dql(4) - dql(5) )
    c2dcl = c2 * dcl
    do m=6,nl
       dhl = dhl + as(m-5) * dql(m)
    enddo

    fl(1) = l1p * dql(1)              -    dhl   * x1l            -    dcl   * x2l
    fl(2) = l1p * dql(2) + ( nx*c2dcl - um*dhl ) * x1l + ( nx*dhl - um*dcl ) * x2l
    fl(3) = l1p * dql(3) + ( ny*c2dcl - vm*dhl ) * x1l + ( ny*dhl - vm*dcl ) * x2l
    fl(4) = l1p * dql(4) + ( nz*c2dcl - wm*dhl ) * x1l + ( nz*dhl - wm*dcl ) * x2l
    fl(5) = l1p * dql(5) + ( st*c2dcl - hm*dhl ) * x1l + ( st*dhl - hm*dcl ) * x2l
    do m=6,nl
       fl(m) = l1p * dql(m) - prim(m) * ( dhl * x1l + dcl * x2l )
    enddo

    dcr = st * dqr(1) - nx * dqr(2) - ny * dqr(3) - nz * dqr(4)
    dhr = af * dqr(1) - ae * ( um * dqr(2) + vm * dqr(3) + wm * dqr(4) - dqr(5) )
    c2dcr = c2 * dcr
    do m=6,nl
       dhr = dhr + as(m-5) * dqr(m)
    enddo

    fr(1) = l1n * dqr(1)              -    dhr   * x1r            -    dcr   * x2r
    fr(2) = l1n * dqr(2) + ( nx*c2dcr - um*dhr ) * x1r + ( nx*dhr - um*dcr ) * x2r
    fr(3) = l1n * dqr(3) + ( ny*c2dcr - vm*dhr ) * x1r + ( ny*dhr - vm*dcr ) * x2r
    fr(4) = l1n * dqr(4) + ( nz*c2dcr - wm*dhr ) * x1r + ( nz*dhr - wm*dcr ) * x2r
    fr(5) = l1n * dqr(5) + ( st*c2dcr - hm*dhr ) * x1r + ( st*dhr - hm*dcr ) * x2r
    do m=6,nl
       fr(m) = l1n * dqr(m) - prim(m) * ( dhr * x1r + dcr * x2r )
    enddo

    return
end subroutine mxdq1
!_____________________________________________________________________!
subroutine mxdq_roe(prim,gykb,dq,f,efix)
    use global_const,only:nl,ns,ms1,beta1,tref,nchem,sml_sss
    implicit none
    integer :: m
    real :: hint,gama,efix,eps
    real :: f(nl),prim(nl),gykb(4),dq(nl),hs(ns),as(ns)
    real :: nx,ny,nz,nt,ct,cgm,cgm1,st
    real :: l1,l4,l5,x1,x2
    real :: dh,dc,c2dc,ae,af
    real :: rm,um,vm,wm,pm,cm,c2,v2,tm,hm

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    call gettmhmgm(prim,tm,hint,gama)
    c2 = gama*pm/rm
    cm = sqrt(c2)
    v2 = um*um + vm*vm + wm*wm

    hm = hint + 0.5*v2

    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)

    ct = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)
	eps = ( efix * cm * cgm )**2

    l1 = abs( ct )
    l4 = abs( ct + cm * cgm )
    l5 = abs( ct - cm * cgm )

    l1 = sqrt( l1*l1 + eps )
    l4 = sqrt( l4*l4 + eps )
    l5 = sqrt( l5*l5 + eps )

    x1 = ( 2.0*l1 - l4 - l5 )/( 2.0 * c2 )
    x2 = ( l4 - l5 )/( 2.0 * cm )

    cgm1 = 1.0/cgm
    nx = nx * cgm1
    ny = ny * cgm1
    nz = nz * cgm1
    st = (ct - nt )* cgm1 ! 减去nt ，在矩阵中不包含nt 
    ae = gama - 1.0

    af = 0.5*ae*v2

    dc = st * dq(1) - nx * dq(2) - ny * dq(3) - nz * dq(4)
    dh = af * dq(1) - ae * ( um * dq(2) + vm * dq(3) + wm * dq(4) - dq(5) )

    do m=6,nl
       dh = dh + as(m-5) * dq(m)
    enddo

    c2dc = c2 * dc

    f(1) = l1 * dq(1)             -    dh   * x1           -    dc   * x2
    f(2) = l1 * dq(2) + ( nx*c2dc - um*dh ) * x1 + ( nx*dh - um*dc ) * x2
    f(3) = l1 * dq(3) + ( ny*c2dc - vm*dh ) * x1 + ( ny*dh - vm*dc ) * x2
    f(4) = l1 * dq(4) + ( nz*c2dc - wm*dh ) * x1 + ( nz*dh - wm*dc ) * x2
    f(5) = l1 * dq(5) + ( st*c2dc - hm*dh ) * x1 + ( st*dh - hm*dc ) * x2

    do m=6,nl
       f(m) = l1 * dq(m) - prim(m) * ( dh * x1 + dc * x2 )
    enddo

    return
end subroutine mxdq_roe
!_____________________________________________________________________!
subroutine flux(prim,gama,nl,nx,ny,nz,nt,f)  
    implicit none
    integer :: m,nl
    real :: rm,um,vm,wm,pm,em,hm,tm,gama
    real :: nx,ny,nz,nt,cit,prim(nl),f(nl)

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    call getem(prim,em,gama)

    cit = nx*um + ny*vm + nz*wm + nt

    f(1) = rm * cit
    f(2) = f(1) * um + nx * pm
    f(3) = f(1) * vm + ny * pm
    f(4) = f(1) * wm + nz * pm
    f(5) = ( em + pm ) * cit - pm * nt

    do m=6,nl
       f(m) = prim(m) * f(1)
    enddo

    return
end subroutine flux
!_____________________________________________________________________!
real function f_zero(x,y)
    implicit none
    real :: x,y
    f_zero = 0.0
    return
end function f_zero
!_____________________________________________________________________!
real function minmod(x,y)
    implicit none
    real :: x,y
    minmod = 0.5 * ( sign(1.0,x) + sign(1.0,y) )*min( abs(x) , abs(y) )
    return
end function minmod
!_____________________________________________________________________!
real function vanleer(x,y)
    implicit none
    real :: x,y,eps
    parameter(eps=1.0e-15)
    vanleer = ( sign(1.0,x) + sign(1.0,y) )*x*y/( abs(x) + abs(y) + eps )
    return
end function vanleer

!_____________________________________________________________________!
real function min_3q(x,y)
    implicit none
    real :: x,y,z,z1,z2,z3,eps
	parameter(eps = 1.0e-16 )
    z  = 0.5 * ( sign(1.0,x) + sign(1.0,y) )*min( abs(x) , abs(y) )
	z1 =1.0/( x*x + eps )**2
	z2 =1.0/( y*y + eps )**2
	z3 =abs(z1-z2)/(z1+z2)
	min_3q = (1.0-z3)* z + z3*( 2.0*x + y )/3.0
    return
end function min_3q
!_____________________________________________________________________!
real function min_3u(x,y)
    implicit none
    real :: x,y,z,z1,z2,z3
	z = 0.5 * ( sign(1.0,x) + sign(1.0,y) )
	if( abs(z) < 1.0e-10 ) then
	    min_3u = 0.0
	else   
	   z2 =( 2.0*x*y + 1.0e-12 )/( x*x + y*y + 1.0e-12 )
	   z3 = ( 2.0*x + y )/3.0
       min_3u = (1.0-z2)* z * min( abs(x) , abs(y) ) + z2 * z3
	endif
    return
end function min_3u
!_____________________________________________________________________!
real function min_van(x,y)
    implicit none
    real :: x,y,z,z1,z2,z3,eps
    parameter(eps=1.0e-15)
	z = 0.5 * ( sign(1.0,x) + sign(1.0,y) )
	if( abs(z) < 1.0e-10 ) then
	    min_van = 0.0
	else   
	   z2 =( 2.0*x*y + 1.0e-12 )/( x*x + y*y + 1.0e-12 )
	   z3 =  2.0*x*y/( abs(x) + abs(y) + eps )
       min_van = z * ( (1.0-z2)* min( abs(x) , abs(y) ) + z2 * z3 )
	endif
    return
end function min_van
!_____________________________________________________________________!
subroutine spectinv
    use global_variables
    implicit none
    integer :: i,j,k,cmethod,ii,jj,kk
    real :: um,vm,wm,cm
    real :: kx,ky,kz,ex,ey,ez,cx,cy,cz,kt,et,ct
    real :: cita,citb,citc,cgma,cgmb,cgmc

!    do k=1,nk-cmethod  !* TGH. 有限差分 *!
!       do j=1,nj-cmethod
!          do i=1,ni-cmethod
     do k=1,nk    !___cic   !*TGH. 注意，原始取法
       do j=1,nj
          do i=1,ni
             um = u(i,j,k)
             vm = v(i,j,k)
             wm = w(i,j,k)
             cm = c(i,j,k)

             call getr0kec_mml(i,j,k,kx,ky,kz,kt,1,0,0) !取网格点上的网格导数
             call getr0kec_mml(i,j,k,ex,ey,ez,et,0,1,0)
             call getr0kec_mml(i,j,k,cx,cy,cz,ct,0,0,1)

             cita = kx*um + ky*vm + kz*wm + kt
             citb = ex*um + ey*vm + ez*wm + et
             citc = cx*um + cy*vm + cz*wm + ct

             cgma = sqrt( kx*kx + ky*ky + kz*kz )
             cgmb = sqrt( ex*ex + ey*ey + ez*ez )
             cgmc = sqrt( cx*cx + cy*cy + cz*cz )

             sra(i,j,k) = abs( cita ) + cm * cgma
             srb(i,j,k) = abs( citb ) + cm * cgmb
             src(i,j,k) = abs( citc ) + cm * cgmc
          enddo
       enddo
    enddo

	do k=1,nk
	do j=1,nj
		do i=0,ni+1,ni+1  !*tgh. I方向的边界
             um = u(i,j,k)
             vm = v(i,j,k)
             wm = w(i,j,k)
             cm = c(i,j,k)
			 ii = max(1,min(ni,i))

             call getr0kec_mml(ii,j,k,kx,ky,kz,kt,1,0,0) !取网格点上的网格导数
             call getr0kec_mml(ii,j,k,ex,ey,ez,et,0,1,0)
             call getr0kec_mml(ii,j,k,cx,cy,cz,ct,0,0,1)

             cita = kx*um + ky*vm + kz*wm + kt
             citb = ex*um + ey*vm + ez*wm + et
             citc = cx*um + cy*vm + cz*wm + ct

             cgma = sqrt( kx*kx + ky*ky + kz*kz )
             cgmb = sqrt( ex*ex + ey*ey + ez*ez )
             cgmc = sqrt( cx*cx + cy*cy + cz*cz )

             sra(i,j,k) = abs( cita ) + cm * cgma
             srb(i,j,k) = abs( citb ) + cm * cgmb
             src(i,j,k) = abs( citc ) + cm * cgmc
		enddo
	enddo
	enddo

	do k=1,nk
	do i=1,ni
		do j=0,nj+1,nj+1  !*tgh. j方向的边界
             um = u(i,j,k)
             vm = v(i,j,k)
             wm = w(i,j,k)
             cm = c(i,j,k)
			 jj = max(1,min(nj,j))

             call getr0kec_mml(i,jj,k,kx,ky,kz,kt,1,0,0) !取网格点上的网格导数
             call getr0kec_mml(i,jj,k,ex,ey,ez,et,0,1,0)
             call getr0kec_mml(i,jj,k,cx,cy,cz,ct,0,0,1)

             cita = kx*um + ky*vm + kz*wm + kt
             citb = ex*um + ey*vm + ez*wm + et
             citc = cx*um + cy*vm + cz*wm + ct

             cgma = sqrt( kx*kx + ky*ky + kz*kz )
             cgmb = sqrt( ex*ex + ey*ey + ez*ez )
             cgmc = sqrt( cx*cx + cy*cy + cz*cz )

             sra(i,j,k) = abs( cita ) + cm * cgma
             srb(i,j,k) = abs( citb ) + cm * cgmb
             src(i,j,k) = abs( citc ) + cm * cgmc
		enddo
	enddo
	enddo

	do j=1,nj
	do i=1,ni
		do k=0,nk+1,nk+1  !*tgh. k方向的边界
             um = u(i,j,k)
             vm = v(i,j,k)
             wm = w(i,j,k)
             cm = c(i,j,k)
			 kk = max(1,min(nk,k))

             call getr0kec_mml(i,j,kk,kx,ky,kz,kt,1,0,0) !取网格点上的网格导数
             call getr0kec_mml(i,j,kk,ex,ey,ez,et,0,1,0)
             call getr0kec_mml(i,j,kk,cx,cy,cz,ct,0,0,1)

             cita = kx*um + ky*vm + kz*wm + kt
             citb = ex*um + ey*vm + ez*wm + et
             citc = cx*um + cy*vm + cz*wm + ct

             cgma = sqrt( kx*kx + ky*ky + kz*kz )
             cgmb = sqrt( ex*ex + ey*ey + ez*ez )
             cgmc = sqrt( cx*cx + cy*cy + cz*cz )

             sra(i,j,k) = abs( cita ) + cm * cgma
             srb(i,j,k) = abs( citb ) + cm * cgmb
             src(i,j,k) = abs( citc ) + cm * cgmc
		enddo
	enddo
	enddo

    return
end subroutine spectinv
!_____________________________________________________________________!
subroutine lxdq(prim,gykb,dq,f)
    use global_variables,only:nl,small,sml_sss
    implicit none
    real,parameter :: og2=0.707106781186548_8
!--  此子程序计算左特征矩阵与向量rq的乘积，结果为f向量
!--  og2为SQRT(2)的倒数
    integer :: m
    real    :: hint,gama
    real    :: f(nl),prim(nl),gykb(4),dq(nl)
    real    :: nx,ny,nz,nt,ct,cgm,cgm1
    real    :: dh,dc,ae,af
    real    :: rm,um,vm,wm,pm,cm,oc,c2,v2,tm,hm
    real    :: vv1,vv2,vv3,dvv1,dvv2,dvv3

    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    call gettmhmgm(prim,tm,hint,gama)
    c2 = gama*pm/rm
    if ( c2 < 0.0 ) then
       write(*,*)'lxdq(prim,gykb,dq,f)'
       write(*,*)'c2<0,gama,pm,rm',gama,pm,rm
    endif

    cm = sqrt(c2)
    oc = 1.0/( cm + small )
    v2 = um*um + vm*vm + wm*wm

    hm = hint + 0.5*v2

    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)

    ct = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

    cgm1 = 1.0/cgm
    nx = nx * cgm1
    ny = ny * cgm1
    nz = nz * cgm1
    ct = ct * cgm1
    ae = gama - 1.0
    af = 0.5*ae*v2

    dc = ct * dq(1) - nx * dq(2) - ny * dq(3) - nz * dq(4)
    dh = af * dq(1) - ae * ( um * dq(2) + vm * dq(3) + wm * dq(4) - dq(5) )

    vv1 = ny * wm - nz * vm
    vv2 = nz * um - nx * wm
    vv3 = nx * vm - ny * um

    dvv1 = ny * dq(4) - nz * dq(3)
    dvv2 = nz * dq(2) - nx * dq(4)
    dvv3 = nx * dq(3) - ny * dq(2)

    f(1) = oc * nx  * ( c2 * dq(1) - dh ) + ( dvv1 - vv1 * dq(1) )
    f(2) = oc * ny  * ( c2 * dq(1) - dh ) + ( dvv2 - vv2 * dq(1) )
    f(3) = oc * nz  * ( c2 * dq(1) - dh ) + ( dvv3 - vv3 * dq(1) )
    f(4) = oc * og2 * ( dh - cm * dc )
    f(5) = oc * og2 * ( dh + cm * dc )

    do m=6,nl
       f(m) = dq(m) - prim(m) * dq(1)
    enddo

    return
end subroutine lxdq
!_____________________________________________________________________!
subroutine rxdq(prim,gykb,dq,f)
    use global_variables,only:nl,small,sml_sss
    implicit none
    real,parameter :: og2=0.707106781186548_8
!--  此子程序计算右特征矩阵与向量rq的乘积，结果为f向量
!--  og2为SQRT(2)的倒数
    integer :: m
    real    :: hint,gama
    real    :: f(nl),prim(nl+1),gykb(4),dq(nl)
    real    :: nx,ny,nz,nt,ct,cgm,cgm1
    real    :: dc,ae,oae
    real    :: rm,um,vm,wm,pm,cm,oc,c2,v2,tm,hm
    real    :: vv1,vv2,vv3,dww1,dww2,dww3,dwwh,dq45
    
    rm = prim(1)
    um = prim(2)
    vm = prim(3)
    wm = prim(4)
    pm = prim(5)

    call gettmhmgm(prim,tm,hint,gama)
    c2 = gama*pm/rm
    if ( c2 < 0.0 ) then
       write(*,*)'rxdq(prim,gykb,dq,f)'
       write(*,*)'c2<0,gama,pm,rm',gama,pm,rm
    endif

    cm = sqrt(c2)
    oc = 1.0/( cm + small )
    v2 = um*um + vm*vm + wm*wm

    hm = hint + 0.5*v2

    nx = gykb(1)
    ny = gykb(2)
    nz = gykb(3)
    nt = gykb(4)

    ct = nx*um + ny*vm + nz*wm + nt
    cgm = max(sqrt(nx*nx + ny*ny + nz*nz),sml_sss)

    cgm1 = 1.0/cgm
    nx = nx * cgm1
    ny = ny * cgm1
    nz = nz * cgm1
    ct = ct * cgm1
    ae = gama - 1.0
    oae = 1.0/ae

    dc = nx * dq(1) + ny * dq(2) + nz * dq(3) + og2 * ( dq(4) + dq(5) )

    vv1 = ny * wm - nz * vm
    vv2 = nz * um - nx * wm
    vv3 = nx * vm - ny * um

    dww1 = nz * dq(2) - ny * dq(3)
    dww2 = nx * dq(3) - nz * dq(1)
    dww3 = ny * dq(1) - nx * dq(2)
    dwwh = vv1 * dq(1) + vv2 * dq(2) + vv3 * dq(3)
    dq45 = dq(4) - dq(5)

    f(1) = oc * dc
    f(2) = oc * um * dc + ( dww1 + og2 * nx * dq45 )
    f(3) = oc * vm * dc + ( dww2 + og2 * ny * dq45 )
    f(4) = oc * wm * dc + ( dww3 + og2 * nz * dq45 )
    f(5) = oc * hm * dc + ( dwwh + og2 * ct * dq45 ) - &
!!           oc * oae * ( nx*dq(1) + ny*dq(2) + nz*dq(3) )
           cm * oae * ( nx*dq(1) + ny*dq(2) + nz*dq(3) )

    do m=6,nl
       f(m) = oc * prim(m) * dc + dq(m)
    enddo

    return

end subroutine rxdq
!_____________________________________________________________________!
subroutine lxdq_p(c2,roc,nx,ny,nz,ddfv,ddfc)
    implicit none
    integer :: m
    real    :: c2,roc,nx,ny,nz,aa
    real    :: ddfv(5),ddfc(5)

	aa = ddfv(1) - ddfv(5)/c2
	ddfc(1) = nx*aa + nz*ddfv(3) - ny*ddfv(4) 
	ddfc(2) = ny*aa + nx*ddfv(4) - nz*ddfv(2)
	ddfc(3) = nz*aa + ny*ddfv(2) - nx*ddfv(3)
	 
	aa = nx*ddfv(2) + ny*ddfv(3) + nz*ddfv(4)
	ddfc(4) =  aa +  ddfv(5)/roc
	ddfc(5) = -aa +  ddfv(5)/roc
	return
end subroutine lxdq_p
!_____________________________________________________________________!
subroutine rxdq_p(rpc,roc,nx,ny,nz,ddfv,ddfc)
    implicit none
    integer :: m
    real    :: rpc,roc,nx,ny,nz
    real    :: ddfv(5),ddfc(5),aa,bb
	aa = (ddfv(4) - ddfv(5) )/2.0
	bb = (ddfv(4) + ddfv(5) )/2.0
	ddfc(1) =  bb*rpc + nx*ddfv(1) + ny*ddfv(2) + nz*ddfv(3)
	ddfc(2) =  aa*nx  + ny*ddfv(3) - nz*ddfv(2) 
	ddfc(3) =  aa*ny  + nz*ddfv(1) - nx*ddfv(3)
	ddfc(4) =  aa*nz  + nx*ddfv(2) - ny*ddfv(1)            
	ddfc(5) =  bb*roc                                     
	return
end subroutine rxdq_p
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

subroutine inviscd3d_cic( limiter,flux_line,flux_type )
!------------------------------------------------------------------------------!
!       在原inviscd3d的基础上进行了改动
!       改动后不使用虚点上的值，且计算了1~max的值（有限差分）
!       本程序和原程序都用了虚点上的值
!       对用到的外部子程序flux_line，目前只对nnd_pv 进行了相应的改动
!  By TU Guohua  2009.2 （还没有完成）
!------------------------------------------------------------------------------!
    use global_variables
    implicit none
    integer :: i,j,k,m,cmethd
    real,external :: limiter
    external  flux_line,flux_type
    real :: fc(1:nl,1:nmax),q_line(1:nl+1,-1:nmax+1,2),trxyz(5,nmax) !*TGH. NL=NM=5 *!

    if ( nchem == 0 .and. nchem_source == 0 ) then
       do i=-1,nmax+1
          do m=1,2
             q_line(nl+1,i,m) = gama
          enddo
       enddo
    endif

    cmethd = 0 + method
!    do k= cmethd,nk-1    !2,nk-1
!       do j= cmethd,nj-1      !2,nj-1
    do k= 1,nk         !___cic
       do j= 1,nj      !___cic
          do i= method-1,ni+1    !1,ni+1
             q_line(1,i,1) = r(i,j,k)
             q_line(2,i,1) = u(i,j,k)
             q_line(3,i,1) = v(i,j,k)
             q_line(4,i,1) = w(i,j,k)
             q_line(5,i,1) = p(i,j,k)
             do m=6,nl
                q_line(m,i,1) = fs(m-5,i,j,k) !* TGH. 化学反应组份 ??? *!
             enddo
             do m=1,nl
                q_line(m,i,2) = q(m,i,j,k)
             enddo
!             if ( nchem /= 0 .or. nchem_source /= 0) q_line(nl+1,i,1) = gm(i,j,k) 
             q_line(nl+1,i,2) = t(i,j,k)
          enddo
          do i=1,ni
             trxyz(1,i) = kcx(i,j,k)
             trxyz(2,i) = kcy(i,j,k)
             trxyz(3,i) = kcz(i,j,k)
             trxyz(4,i) = kct(i,j,k)
             trxyz(5,i) = vol(i,j,k)
          enddo
          call flux_line(xk,xb,ni,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz,q_line,fc) 
!          do i=cmethd,ni-1     !2,ni-1
          do i=1,ni     !___cic
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) +  fc(m,i)
             enddo
          enddo

       enddo
    enddo

!    do k=cmethd,nk-1        !2,nk-1
!       do i=cmethd,ni-1     !2,ni-1
    do k=1,nk        !___cic
       do i=1,ni     !___cic
          do j=-1,nj+1  !1,nj+1
             q_line(1,j,1) = r(i,j,k)
             q_line(2,j,1) = u(i,j,k)
             q_line(3,j,1) = v(i,j,k)
             q_line(4,j,1) = w(i,j,k)
             q_line(5,j,1) = p(i,j,k)
             do m=6,nl
                q_line(m,j,1) = fs(m-5,i,j,k)
             enddo
             do m=1,nl
                q_line(m,j,2) = q(m,i,j,k)
             enddo

!             if ( nchem /= 0 .or. nchem_source /= 0) q_line(nl+1,j,1) = gm(i,j,k) 
             q_line(nl+1,j,2) = t(i,j,k)
          enddo

          do j=1,nj
             trxyz(1,j) = etx(i,j,k)
             trxyz(2,j) = ety(i,j,k)
             trxyz(3,j) = etz(i,j,k)
             trxyz(4,j) = ett(i,j,k)
             trxyz(5,j) = vol(i,j,k)
          enddo

          call flux_line(xk,xb,nj,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz,q_line,fc)

!          do j=cmethd,nj-1     !2,nj-1
          do j=1,nj     !___cic
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) +  fc(m,j)
             enddo
          enddo
       enddo
    enddo

!    do j=cmethd,nj-1       !2,nj-1
!       do i=cmethd,ni-1    !2,ni-1
    do j=1,nj       !___cic
       do i=1,ni    !___cic
          do k= -1,nk+1    !1,nk+1
             q_line(1,k,1) = r(i,j,k)
             q_line(2,k,1) = u(i,j,k)
             q_line(3,k,1) = v(i,j,k)
             q_line(4,k,1) = w(i,j,k)
             q_line(5,k,1) = p(i,j,k)
             do m=6,nl
                q_line(m,k,1) = fs(m-5,i,j,k)
             enddo
             do m=1,nl
                q_line(m,k,2) = q(m,i,j,k)
             enddo
!             if ( nchem /= 0 .or. nchem_source /= 0) q_line(nl+1,k,1) = gm(i,j,k) 
             q_line(nl+1,k,2) = t(i,j,k)
          enddo
          do k=1,nk
             trxyz(1,k) = ctx(i,j,k)
             trxyz(2,k) = cty(i,j,k)
             trxyz(3,k) = ctz(i,j,k)
             trxyz(4,k) = ctt(i,j,k)
             trxyz(5,k) = vol(i,j,k)
          enddo
          call flux_line(xk,xb,nk,nmax,nl,method,cmethd,flux_type,limiter,efix,trxyz,q_line,fc)

!          do k=cmethd,nk-1      !2,nk-1
          do k=1,nk      !___cic
             do m=1,nl
                dq(m,i,j,k) = dq(m,i,j,k) + fc(m,k)
             enddo
          enddo
       enddo
    enddo

    return
end subroutine inviscd3d_cic

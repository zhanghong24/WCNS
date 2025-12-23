module define_precision_mod
!+--------------------------------------------------------------------+
!|   module define_precision_mod is developed to define single or     |
!| double real precision kind for selecting .                         |
!|   called by : globle share                                         |
!|   authors   : Deng Xiaobing                                        |
!|   input     : no                                                   |
!|   output    : no                                                   |
!|   modified history  :                                              |
!|     date       programmer    description                           |
!|   2003-06-17   Deng Xiaobing 1 to F90 free format                  |
!+--------------------------------------------------------------------+
   implicit none
!
   integer, parameter :: single_prec = 4  !!KIND(1.0e0)
   integer, parameter :: double_prec = 8  !!KIND(1.0d0)
!   integer, parameter ::        prec = single_prec
   integer, parameter ::        prec = double_prec
!
end module define_precision_mod
!_____________________________________________________________________!
module type_of_data
use define_precision_mod
!   һάʵ�Ͷ�̬����
    type real1d
        real,pointer::a1d( : )
    end type real1d
!   ��άʵ�Ͷ�̬����
    type real2d
        real,pointer::a2d( :,: )
    end type real2d
!   ��άʵ�Ͷ�̬����
    type real3d
        real,pointer::a3d( :,:,: )
    end type real3d
!   ��άʵ�Ͷ�̬����
    type real4d
        real,pointer::a4d( :,:,:,: )
    end type real4d
!   ��άʵ�Ͷ�̬����
    type real5d
        real,pointer::a5d( :,:,:,:,: )
    end type real5d
!   һά���Ͷ�̬����
    type int1d
        integer,pointer::a1d( : )
    end type int1d
!   ��ά���Ͷ�̬����
    type int2d
        integer,pointer::a2d( :,: )
    end type int2d
!   ��ά���Ͷ�̬����
    type int3d
        integer,pointer::a3d( :,:,: )
    end type int3d
!   ��ά���Ͷ�̬����
    type int4d
        integer,pointer::a4d( :,:,:,: )
    end type int4d

    type cbc
        integer :: bctype             !�߽�����
        integer :: s_nd,t_nd          !�߽��淽��:1,2,3��Ӧ��i,j,k
        integer :: s_lr,t_lr          !���ұ߽�-1,1��Ӧ�����ұ߽�
        integer :: nbs                !���
        integer :: nbt                !��Ӧ�ڿ����߽�����,ָ����ӦĿ��������Ϣ�ڵڼ�������
        integer :: s_fix,t_fix        !�̶�����(fixed_coor)
        integer :: s_st(3),t_st(3)    !��ʼ������(�ɶ�����)
        integer :: s_ed(3),t_ed(3)    !��ֹ������(�ɶ�����)
        integer :: s_lr3d(3),t_lr3d(3)
        integer,pointer :: image(:,:,: ),jmage(:,:,: ),kmage(:,:,: ) !��Ӧ�ڿɶ�����֮��ı任����
        integer :: bc_par_len         !ָ���߽��������
        real,pointer :: bc_par( : )   !ָ���߽����
        real,pointer :: bc_par1( : , : , : , : )  !ָ���߽����
       !  ƴ������߽�����
!        integer,pointer :: flag_ppt(:,:,:)  !��ֵ
!		integer,pointer :: fmage(:,:,: ) !��Ӧ�ڿɶ�����֮��ı任����
!		integer,pointer :: bmage(:,:,: ) !��Ӧ�ڿɶ�����֮��ı任����
!		integer,pointer :: imag1(:,:,: ),jmag1(:,:,: ),kmag1(:,:,: ) !��Ӧ�ڿɶ�����֮��ı任����
!		integer,pointer :: imag2(:,:,: ),jmag2(:,:,: ),kmag2(:,:,: ) !��Ӧ�ڿɶ�����֮��ı任����
!		integer,pointer :: imag3(:,:,: ),jmag3(:,:,: ),kmag3(:,:,: ) !��Ӧ�ڿɶ�����֮��ı任����
!        real   ,pointer :: qip_mml(:,:,:,:) ! ��ֵϵ��
!        real(prec) , pointer :: qip_mml(:,:,:,:) ! ��ֵϵ��

        !! ���м�����Ҫ
        integer :: ibcwin
        real,pointer :: qpv(:,:,:,:)
        real,pointer :: qpvpack(:,:,:,:)
        real,pointer :: dqv(:,:,:,:)
        real,pointer :: dqvpack(:,:,:,:)
        
        real,pointer :: vol(:,:,:)
        real,pointer :: volpack(:,:,:)
        real,pointer :: sxyz(:,:,:,:)
        real,pointer :: sxyzpack(:,:,:,:)
        real,pointer :: dis(:,:,:)
        real,pointer :: dispack(:,:,:)
        
        real,pointer :: qpv_t(:,:,:,:)
        real,pointer :: qpvpack_t(:,:,:,:)
        real,pointer :: dqv_t(:,:,:,:)
        real,pointer :: dqvpack_t(:,:,:,:)
        
        integer :: iused
    end type cbc
!

    type bcregion
        integer :: nregions
        type(cbc),pointer :: bc(:)
        
        integer,pointer :: bcindexs(:) 
    end type bcregion

!   �㷨�ṹ
 
    type group_type
        integer :: ngzone,ngstep
        integer,pointer :: s_nb(:),e_nb(:)
    end type group_type

    type fffff_type
        integer:: ntmst,nstart,ndisk,nimax,nreset
        integer:: nust,nvis,nlamtur,nchem,ntmodel,nwallfun,ntrans,nwtmax,ncmpcor  !�Ƕ�����ճ�ԣ���������ѧ��Ӧ�����¶�ģ��
        integer:: nlhs,nscheme,nlimiter,nsmooth,nflux
        integer:: nout,nmax
        integer:: ns,nf                            !ʵ�ʵ���Ԫ������Ӧ��
        integer:: nc                               !���ǻ�ѧ��Ӧ�����غ����ķ�������
        integer:: ne                               !���Ǽ��¶�ģ�����غ����ķ�������
        integer:: nchem_source,nchem_rad           !�Ƿ���㻯ѧ��ӦԴ��,Դ���װ뾶
        real:: ratio_mref,kmaxlim,xtrans                  !�ο�������(������)

        character(len=120) :: gasmodel,nameturb
        real :: gama,prl,prt
        integer :: num_variable                    ! ��������
        integer,pointer,dimension( : ) :: varseq   !��������Ӧ�����
        character(len=120),pointer,dimension( : ) :: varname !����������

        real,pointer,dimension( : ) :: ws,ms      !����Ԫ������(������),����Ԫ������(������)
        real,pointer,dimension( : ) :: ws1,ms1    !����Ԫ�������ĵ���(������),����Ԫ�������ĵ���(������)
        real,pointer,dimension( : ) :: cn_init    !����ְٷֱ�,����ת���ܺ͵�������

        real:: cfl,timedt,wholetime,phydt,tol,efix,csrv

!		real(prec):: phytime,phydtime  !*tgh. ����һ��phytime

        real:: cq_lam,cq_tur,visc
        real:: xk,xb,c_k2,c_k4


    end type fffff_type

!==================!
!*tgh. MaoMLΪ��ȫ�������Խӵ�����ӵ�Tpye
	type singularity
	    integer ::singularity
		integer,pointer ::singularities(:,:)
	end type singularity

	type match_point
	    integer :: npp               !һ�������ռ���ɼ��������������
		integer,pointer ::Nbijk(:,:) !��Ӧ�ڿɶ�����֮��ı任�����һ���ű��ǵ�����ţ���1��npp���ڶ����ű���1��4������1��Ӧ��ţ�
                                     !2��4��Ӧ��Ӧ�ĸ��ű���ţ�����Nbijk (:,1)��Nb��Nbijk (:,2)��i��Nbijk (:,3)��i��Nbijk (:,4)��k		
	end type match_point
!*tgh. end. MaoMLΪ��ȫ�������Խӵ�����ӵ�Tpye
!==================!


end module type_of_data
!+--------------------------------------------------------------------+
module global_constants_mod
!+--------------------------------------------------------------------+
!|   module global_constants_mod is developed to define and globally  |
!| share constants  .                                                 |
!|   called by : globle share                                         |
!|   authors   : Yuan XianXu                                          |
!|   input     : no                                                   |
!|   output    : no                                                   |
!|   modified history  :                                              |
!|     date       programmer    description                           |
!|   2003-06-17   Yuan XianXu   1 to F90 free format                  |
!+--------------------------------------------------------------------+
   implicit none
!
! constants for equation dimension
   integer,parameter :: md = 5
!
! constant used in sutherlan' law,=110.3 in NSMB5.0
   real, parameter :: Coo = 110.4
!
! constants for flowfield
   real, parameter :: pzero =  1.0e-6, dzero = 1.0e-6
   real, parameter :: zero = 0.0, nonzero = 1.0e-30
   real, parameter :: half = 0.5, one = 1.0, two = 2.0, big = 1.0e30
!
! constants for flowfield implicit residual smoothing
! pasi = [ 0.125 - 0.25 ]
   real, parameter :: dcr = 2.0 , pasi = 0.125
!
! not changed constants for Baldwin-Lomax turbulence model
   real, parameter :: Ap = 26.0 , Ccp = 1.6 , Ckleb = 0.3
   real, parameter :: Kinner = 0.4 , Kouter = 0.0168
   real, parameter :: peak_cut = 0.9 , mutexp = 1.0
   real, parameter :: CuserLim = 2.4
!
   real, parameter :: eii = 0.3 , ejj = 0.3 , ekk = 0.3
!
!
end module global_constants_mod
!_____________________________________________________________________!
module global_const
    use define_precision_mod
    implicit none
    
    real,parameter :: small=1.0e-38
    real,parameter :: large=1.0e+38
    
    real,parameter :: sml_sss=1.0e-30
    real,parameter :: sml_vol=1.0e-30
    
    real,parameter :: pole_vol=1.0e-24
    
    real,parameter :: pmin_limit=1.0e-8
    real,parameter :: pmax_limit=1.0e+3
    
    real,parameter :: rmin_limit=1.0e-8
    real,parameter :: rmax_limit=1.0e+3
    
    real,parameter :: tmin_limit=1.0e-8
    
    real :: rjmk,sc_l,sc_t,pai
    integer :: cyc(3,2),iii(3,3),nm,nl,nvst,nved,mb_varnumber

    parameter(nm=5,nl=5,nvst=5,nved = 5,mb_varnumber = 5)
    parameter(rjmk=8.31434)  !ͨ�����峣������λΪJ/(Mol*K)      

    parameter(sc_l=0.5,sc_t=0.5)
    parameter(pai=3.14159265)
    integer :: nwholefield                           !���������ֲ���(0,�볡;1,ȫ��)
    integer :: nblocks,nbgmax,nforce,nwerror,method  !���������������������������������������¼�в��,��ɢ����
    integer :: nmethod,ndim                          !�㷨������,�ռ�ά��
    integer :: ntmst,nstart,ndisk,nreset,nomax,ntemp_W,nplot
    integer :: nvis,nlamtur,nchem,ntmodel,nwallfun,ntrans,nwtmax,ncmpcor,ndes   !�Ƕ�����ճ�ԣ���������ѧ��Ӧ�����¶�ģ��
	real :: kmaxlim,xtrans
    integer :: nlhs,nscheme,nlimiter,newton,nsmooth,nflux
    integer :: nout,nmax
    integer :: ns,nf   !ʵ�ʵ���Ԫ������Ӧ��
!    integer :: nl      !ʵ�ʵĿ��ǻ�ѧ��Ӧ�����غ���ܷ�������
    integer :: nc      !���ǻ�ѧ��Ӧ�����غ����ķ�������
    integer :: ne      !���Ǽ��¶�ģ�����غ����ķ�������
    integer :: nchem_source,nchem_rad !�Ƿ���㻯ѧ��ӦԴ��,Դ���װ뾶
    real :: mref,mref1                !�ο�������(������),�ο�����������
    real :: beta1,beta2,beta3         !�����ٻ���ϵ��
    real :: gama,prl,prt
    character(len=120) :: gasmodel,nameturb

	real,dimension(20) :: Turb_const

    character(len=120),pointer,dimension( : ) :: varname !����������
    real,pointer,dimension( : ) :: ws,ms                 !����Ԫ������(������),����Ԫ������(������)
    real,pointer,dimension( : ) :: ws1,ms1               !����Ԫ�������ĵ���(������),����Ԫ�������ĵ���(������)
    real,pointer,dimension( : ) :: cn_init               !����ְٷֱ�,����ת���ܺ͵�������

    character(len=120) :: flowname,gridname,bcname,forcename,errhis,tecname
    real :: roo,uoo,voo,woo,poo,eoo,coo,hoo,too,moo,gmoo,height
    real :: koo,omgoo,goo,muoo,vtoo    ! ��������
    real :: q1oo,q2oo,q3oo,q4oo,q5oo
    real :: visloo
    real :: reynolds,re,rvl,rvl1,vl,vl1
    real :: attack,sideslip,tref,twall,pref,rref,vref
    real :: sref,lfref,lref,xref,yref,zref
    real :: cfl,timedt,wholetime,efix,csrv,CPU_TIME_USED
    real :: cfx,cfy,cfz,cmx,cmy,cmz,cl,cd,xcp
    real :: cq_lam,cq_tur,visc
    real :: pst,rst
    real :: xk,xb,c_k2,c_k4
!!!!!!!!!!!!!!!!!!!!
	integer :: nijk2d,nijk2nd  !*tgh ��ʱ�ӣ�<=nijk2d----��ά��<=nijk2nd----2��
	real :: timedt_rate  !�������ʱ�䲽��Ϊ��Сʱ�䲽����timedt_rate��
	real :: timedt_turb  !����ģ��ʱ�䲽��������ʱ�䲽��֮��
!!!!!!!!!!!!!!!!
end module global_const
!_____________________________________________________________________!
module global_variables
    use define_precision_mod
    use type_of_data
    use global_const
!   Ϊ�˲�ʹ�ӳ����ڲ����������Ķ���������ȫ�ֱ���ͳһ��ǰ׺Ϊmb_ 
    type(bcregion),pointer :: mb_bc(:)
    integer     ,pointer   :: mb_dim(:,:)                 !ά��
    type(real3d),pointer   :: mb_x(:),mb_y(:),mb_z(:)     !����
!    type(real3d),pointer   :: mb_x1(:),mb_y1(:),mb_z1(:)  !����
    type(real3d),pointer   :: mb_vol(:),mb_volt(:)        !���
    type(real3d),pointer   :: mb_kcx(:),mb_kcy(:),mb_kcz(:),mb_kct(:)
    type(real3d),pointer   :: mb_etx(:),mb_ety(:),mb_etz(:),mb_ett(:)
    type(real3d),pointer   :: mb_ctx(:),mb_cty(:),mb_ctz(:),mb_ctt(:)
    
    !!������ֵ��Ҫ(lhy)
    type(int3d),pointer    :: mb_flg(:,:)
    integer :: nbself
    
    type(real3d),pointer   :: mb_r(:),mb_u(:),mb_v(:),mb_w(:),mb_p(:)
    type(real4d),pointer   :: mb_fs(:)                 !��ѧ��Ӧ���

!*TGH. ����Ϊ����ӭ���ʽʱʹ�ü����غ���ʹ��
!*TGH. ��ʱ��û������
    type(real3d),pointer   :: mb_vol_l(:),mb_volt_l(:)        !���
    type(real3d),pointer   :: mb_vol_r(:),mb_volt_r(:)        !���
    type(real3d),pointer   :: mb_kcx_l(:),mb_kcy_l(:),mb_kcz_l(:),mb_kct_l(:)
    type(real3d),pointer   :: mb_kcx_r(:),mb_kcy_r(:),mb_kcz_r(:),mb_kct_r(:)
    type(real3d),pointer   :: mb_etx_l(:),mb_ety_l(:),mb_etz_l(:),mb_ett_l(:)
    type(real3d),pointer   :: mb_etx_r(:),mb_ety_r(:),mb_etz_r(:),mb_ett_r(:)
    type(real3d),pointer   :: mb_ctx_l(:),mb_cty_l(:),mb_ctz_l(:),mb_ctt_l(:)
    type(real3d),pointer   :: mb_ctx_r(:),mb_cty_r(:),mb_ctz_r(:),mb_ctt_r(:)

    real(prec),pointer,dimension( :,:,: ) :: vol_l,volt_l,vol_r,volt_r
    real(prec),pointer,dimension( :,:,: ) :: kcx_l,kcy_l,kcz_l,kct_l,kcx_r,kcy_r,kcz_r,kct_r
    real(prec),pointer,dimension( :,:,: ) :: etx_l,ety_l,etz_l,ett_l,etx_r,ety_r,etz_r,ett_r
    real(prec),pointer,dimension( :,:,: ) :: ctx_l,cty_l,ctz_l,ctt_l,ctx_r,cty_r,ctz_r,ctt_r
!*TGH. endΪ����ӭ���ʽʱʹ�ü����غ���ʹ��
    

    type(real3d),pointer   :: mb_t(:),mb_c(:)          !�¶�,����

    type(real3d),pointer   :: mb_sra(:),mb_srb(:),mb_src(:) !�װ뾶
    type(real3d),pointer   :: mb_srva(:),mb_srvb(:),mb_srvc(:) !�װ뾶
    type(real4d),pointer   :: mb_srs(:)                !��ѧ��ӦԴ���װ뾶
    type(real3d),pointer   :: mb_dtdt(:)               !ʱ�䲽��
    type(real3d),pointer   :: mb_visl(:),mb_vist(:)    !����������ճ��ϵ��
    type(real4d),pointer   :: mb_q(:),mb_dq(:)         !�غ������������
    type(real4d),pointer   :: mb_qke(:),mb_dqke(:)     !��������������
    type(real4d),pointer   :: mb_Turb_spec(:)		        
    type(real3d),pointer   :: mb_distance(:)           !  Two Equation Model Variables
    type(real3d),pointer   :: mb_dgrid(:)   
    type(real3d),pointer   :: mb_trancoe(:)            !  ����ת��ϵ��
    integer,pointer        :: mb_istran(:)             !  �����Ƿ���Ҫת��

    ! *** ˫ʱ�䲽���� ***
    type(real4d),pointer   :: mb_qnc(:),mb_qmc(:)
    integer :: ndualtst,nsubstep,nstopsub,nsubstmx
    real :: reftime,dtdts,ttdts,tolsub

    type(real4d),pointer   :: mb_rhs1(:)               !  �Ҷ���  cic

    type(group_type),pointer :: mb_group(:)
    type(fffff_type),pointer :: mb_control(:)
    integer,pointer :: mb_par_content(:)
    real,pointer :: qb_seq( :,: )                      !ÿ�ֿ��Ʒ��̼���ֵı���ֵ
!    integer :: mb_varnumber
    character(len=120) mb_varname(100)
    integer :: bc_par_len                              !ָ���߽��������
    real,pointer :: bc_par( : )                        !ָ���߽����

    integer:: ni,nj,nk
    real,pointer,dimension( :,:,: )   :: x,y,z
    real,pointer,dimension( :,:,: )   :: vol,volt
    real,pointer,dimension( :,:,: )   :: kcx,kcy,kcz,kct
    real,pointer,dimension( :,:,: )   :: etx,ety,etz,ett
    real,pointer,dimension( :,:,: )   :: ctx,cty,ctz,ctt
    real,pointer,dimension( :,:,: )   :: r,u,v,w,p
    real,pointer,dimension( :,:,: )   :: t,c
!    real,pointer,dimension( :,:,: )   :: h,gm
    real,pointer,dimension( :,:,: )   :: sra,srb,src,srva,srvb,srvc
    real,pointer,dimension( :,:,:,: ) :: srs
    real,pointer,dimension( :,:,: )   :: dtdt
    real,pointer,dimension( :,:,: )   :: visl
    real,pointer,dimension( :,:,: )   :: vist
    real,pointer,dimension( :,:,:,: ) :: fs
    real,pointer,dimension( :,:,:,: ) :: q,dq
    real,pointer,dimension( :,:,:,: ) :: qke,dqke
    real,pointer,dimension( :,:,: )   :: ds_turbulence !  Normal Distence from Solid Surface
    real,pointer,dimension( :,:,:,: ) :: spec
    real,pointer,dimension( :,:,:,:,:):: dmudxyz
    real,pointer,dimension( :,:,:,:)  :: dudxyz,dvdxyz,dwdxyz,drdxyz
    real,pointer,dimension( :,:,:,:)  :: viseq         !  ������Чճ��ϵ��
    real,pointer,dimension( :,:,: )   :: trancoe       !  ����ת�����
    real,pointer,dimension( :,:,: )   :: dgrid       !  DES����߶�

    real,pointer,dimension( :,:,:,: ) :: rhs1          !  �Ҷ���  cic

    real,pointer ::  gcl0(:,:,:,:)
    type(real4d),pointer ::  mb_gcl0(:)

	!˫ʱ�䲽����
    character(len=120) :: lastflow,meshbak
    integer::  nust,id_start,it_min,it_max,time_ac  !˫ʱ�䲽����������С������ӵ���������ʱ�侫��
    real   ::  phydt,cfl_all,tol   !����ʱ�䲽����ͳһCfl���������ݲ�
	real(prec)  :: phytime,phydtime    !*TGH. ע�⣬��������ʱ�������ʱ�䲽��
	
	!RoeMx
    real,pointer :: rpmin(:) !!lhy

	
	real,pointer,dimension( :,:,:,:,: ) :: q_unsteady ! ��������
    type(real5d),pointer   :: mb_q_unsteady(:)        ! ��������

	!������
    integer::  movetype,iv_start,moveaxes,movegrid
	  real   ::  location,amplitude,frequency,mass,i_mx,i_my,i_mz
	  real   ::  p_ang0,y_ang0,r_ang0,x_c0,y_c0,z_c0

	integer :: gcl,cic1,cic2,cic3,cic4,cic5 !*TGH. �����غ��ɣ������߽����
!*TGH. GCL-�������غ��㷨
!*TGH. CIC1-�����Խ�
!*TGH. CIC2-���������߽�
!*TGH. CIC3-�Գ��������߽�
!*TGH. CIC4-Զ�������߽�
!*TGH. CIC5-����Ԥ��������

!==================!
!*tgh. MaoMLΪ��ȫ�������Խӵ������
!   ȫ�ֱ����Ķ���
   	integer:: nptt                 !�������ռ俴���Խ����ϵ��ܵ���
	integer:: singles              ! ���ڱ���ͶԽ����ϵ���������Щ���еĿ鲢����ʾ���ڱ����ϡ�           
    type(match_point),pointer:: cdate(:)   !�ýṹ�ķ�ʽ���������λ�õĽǶȼ�¼�Խ����ϵ��λ��
	type(singularity),pointer:: single(:)  !���ڱ���ͶԽ����ϵ㣬����Щ���еĿ鲢����ʾ���ڱ����ϡ�
	integer:: connect_point,connect_order   !*TGH. ��ȡ/����Խ���ϸ��Ϣ
    character(len=120) :: conpointname
	real :: dis_tao
!*tgh. end. MaoMLΪ��ȫ�������Խӵ������
!==================!

end module global_variables
!+--------------------------------------------------------------------+
module global_matrixes__mod
!+--------------------------------------------------------------------+
!|   module global_matrixes__mod is developed to define and globally  |
!| share constant unit matrixes for select points and faces by integer|
!| constant nface.                                                    |
!|   called by : globle share                                         |
!|   authors   : Yuan XianXu                                          |
!|   input     : no                                                   |
!|   output    : no                                                   |
!|   modified history  :                                              |
!|     date       programmer    description                           |
!|   2003-06-17   Yuan XianXu   1 to F90 free format                  |
!+--------------------------------------------------------------------+
   use define_precision_mod
   implicit none
!
! get region direction decided by nafce
   real(prec),dimension(1:6) :: ndir_get = (/-1,1,-1,1,-1,1/)
!
! get dimension(1->i ; 2->j ; 3->k) decided by nafce
   integer,dimension(1:6) :: ndim_get = (/1,1,2,2,3,3/)
!
! 4 constant matrix used for integration aerodynamic force
! define integral region decided by nafce
   integer,dimension(1:3,1:6) :: ndim_add = reshape((/0,1,1,&
                                              0,1,1,&
                                              1,0,1,&
                                              1,0,1,&
                                              1,1,0,&
                                              1,1,0/), (/3,6/))
! i,j,k index for point 2 of integral quadrangle decided by nafce
   integer,dimension(1:3,1:6) :: ndim_p2  = reshape((/0,1,0,&
                                              0,1,0,&
                                              1,0,0,&
                                              1,0,0,&
                                              1,0,0,&
                                              1,0,0/), (/3,6/))
! i,j,k index for point 3 of integral quadrangle decided by nafce
   integer,dimension(1:3,1:6) :: ndim_p3  = reshape((/0,1,1,&
                                              0,1,1,&
                                              1,0,1,&
                                              1,0,1,&
                                              1,1,0,&
                                              1,1,0/), (/3,6/))
! i,j,k index for point 4 of integral quadrangle decided by nafce
   integer,dimension(1:3,1:6) :: ndim_p4  = reshape((/0,0,1,&
                                              0,0,1,&
                                              0,0,1,&
                                              0,0,1,&
                                              0,1,0,&
                                              0,1,0/), (/3,6/))
!
! 2 constant matrix used for boundary condition manipulate
! i,j,k index for left points near boundary
   integer,dimension(1:3,1:6) :: ndim_bL = reshape((/1,0,0,&
                                            -1,0,0,&
                                             0,1,0,&
                                             0,-1,0,&
                                             0,0,1,&
                                             0,0,-1/), (/3,6/))
! i,j,k index for right points near boundary
   integer,dimension(1:3,1:6) :: ndim_bR = reshape((/-1,0,0,&
                                              1,0,0,&
                                              0,-1,0,&
                                              0,1,0,&
                                              0,0,-1,&
                                              0,0,1/), (/3,6/))
!
! 4 constant matrix used in Baldwin-Lomax model
! i,j,k index transformation from l,m,n used in BL model
   integer,dimension(1:3,1:6) :: i_tm = reshape((/0,0,1,&
                                          0,0,1,&
                                          0,1,0,&
                                          0,1,0,&
                                          0,1,0,&
                                          0,1,0/), (/3,6/))
   integer,dimension(1:3,1:6) :: j_tm = reshape((/0,1,0,&
                                          0,1,0,&
                                          0,0,1,&
                                          0,0,1,&
                                          1,0,0,&
                                          1,0,0/), (/3,6/))
   integer,dimension(1:3,1:6) :: k_tm = reshape((/1,0,0,&
                                          1,0,0,&
                                          1,0,0,&
                                          1,0,0,&
                                          0,0,1,&
                                          0,0,1/), (/3,6/))
   integer,dimension(1:6) :: ijk_add = reshape((/1,-1,1,-1,1,-1/), (/6/))
!
! 3 constant matrix used in body_moving
! transformation matrix used for body ROLLing
   integer,dimension(1:3,1:3) :: move_dim = reshape((/1,0,0,&
                                             0,-1,0,&
                                             0,0,-1/), (/3,3/))
   integer,dimension(1:3,1:3) :: move_add = reshape((/0,0,1,&
                                              0,0,-1,&
                                              -1,0,0/), (/3,3/))
   real(prec),dimension(1:3,1:3)    :: move_deg = reshape((/0,0,0,&
                                                 90,270,0,&
                                                    0,0,0/), (/3,3/))
! interpolation coefficient matrix used by bc
! for singular axes
! for wall pressure and temperature
! for wall velocity and open bc
   real(prec),dimension(1:2,1:3) :: coe_interpolation = reshape((/1.0,0.0,&
                                                    2.0,-1.0,&
                                                    4.0/3.0,-1.0/3.0/), (/2,3/))
! 3 constant matrix used in patched-grid searching
   integer,dimension(1:3,1:6) :: mpg = reshape((/0,0,0,&
                                         1,0,1,&
                                         0,1,0,&
                                         0,1,1,&
                                         0,0,0,&
                                         1,0,0/), (/3,6/))
   integer,dimension(1:4) :: mquad = (/1,4,7,10/)
   integer,dimension(1:4) :: mnode = (/2,3,4,1/)
   integer,dimension(1:3) :: mfrom = (/1,2,4/)
   real(prec),dimension(1:2,1:4) :: quadxe = reshape((/0,0,&
                                               1,0,&
                                               1,1,&
                                               0,1/), (/2,4/))
!
end module global_matrixes__mod
!_____________________________________________________________________!
module stress_module
    real :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real :: txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
end module stress_module
!_____________________________________________________________________!
module resdual_module
    integer,pointer :: nres(:,:)
    real,pointer :: res(:),resmax(:)

    real :: dres_add,resmax_total
    integer :: n_error(5)
    real :: resall_dts0

end module resdual_module
!_____________________________________________________________________!
module duvwt_module
    real :: ukc,uet,uct,vkc,vet,vct,wkc,wet,wct,tkc,tet,tct
end module duvwt_module
!_____________________________________________________________________!
module geometry_module
    real :: kx,ky,kz,ex,ey,ez,cx,cy,cz,vjacob
end module geometry_module
!_____________________________________________________________________!
module sarea_module
    implicit none
    real :: skxp,skyp,skzp,skxm,skym,skzm
    real :: sexp,seyp,sezp,sexm,seym,sezm
    real :: scxp,scyp,sczp,scxm,scym,sczm
    real :: vvv1
end module sarea_module
!_____________________________________________________________________!

!_____________________________________________________________________!
module duvwt_all_field !����ճ����ʹ�õĴ�����
	USE DEFINE_PRECISION_MOD
	real,pointer,dimension(:,:,:,:) :: duvwt_mid,duvwt
END module duvwt_all_field

!_____________________________________________________________________!
module rk_3s_global
	use define_precision_mod
	use type_of_data
	implicit none 
	
	type(real4d),pointer :: mb_q00(:)
	real(prec),pointer :: q00(:,:,:,:)

end module rk_3s_global
!_____________________________________________________________________!

!=====================================================================!
!=====================================================================!
module store_GCL
    use define_precision_mod
    use type_of_data
    use global_variables, only : gcl0, mb_gcl0  ! <--- 关键：从上面引用进来
    implicit none
    public :: gcl0, mb_gcl0                     ! <--- 关键：设为 Public，相当于“转发”给别人用
end module store_GCL

!=====================================================================!
!=====================================================================!
!=============================================================================!

module store_rhs_for_GS
	use define_precision_mod
	use type_of_data
    implicit none
    real,pointer ::  rhs_nb(:,:,:,:)
    type(real4d),pointer ::  mb_rhs0(:)
end module store_rhs_for_GS
!=============================================================================!


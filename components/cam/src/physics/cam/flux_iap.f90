      module flux_iap
        !====================================================
        !这个module是大气所的通量算法
        !作者：马钰斌，程雪玲，吴琳，靳江波，李奇龙，曾庆存
        !2022-12-14
        !====================================================
        use shr_kind_mod , only : r8=> shr_kind_r8
        use physics_iap  , only : cal_grv,cal_psi_u,cal_psi_t,cal_psi_q,cal_qsat,cal_qsee
        implicit none 
        private 
        public cal_flux_iap 

       contains


        subroutine cal_flux_iap(z,u,v,t,q,u_s,v_s,t_s,q_s,lat,rho_a,tau,shf,lhf ) 
        !使用整体算法计算通量


        !输入变量
                real(r8) ,intent(in)        :: rho_a !空气密度

                real(r8) ,intent(in)        ::  z    ! 高度z

                real(r8) ,intent(in)        ::  u    ! 在高度z处的 风速u
                real(r8) ,intent(in)        ::  v    ! 在高度z处的 风速v
                real(r8) ,intent(in)        ::  t    ! 在高度z处的 温度
                real(r8) ,intent(in)        ::  q    ! 在高度z处的 比湿

                real(r8) ,intent(in)        ::  u_s    ! 表面的 风速u
                real(r8) ,intent(in)        ::  v_s    ! 表面的 风速v
                real(r8) ,intent(in)        ::  t_s    ! 表面的 温度
                real(r8) ,intent(in)        ::  q_s    ! 表面的 比湿
                real(r8) ,intent(in)        :: lat     !纬度（度）

        !中间变量
                integer         :: i           ! 循环索引

                real(r8)        ::  von         !冯-卡曼常数(一般取为0.4) 
                real(r8)        ::  z0         !动量粗糙度
                real(r8)        ::  z0t        !温度粗糙度
                real(r8)        ::  z0q        !比湿粗糙度
                real(r8)        :: charn        !计算粗糙度的系数charkon常数
                real(r8)        :: visa         !计算粗糙度的中间变量 
                real(r8)        :: u_min        !最小风速差
                real(r8)        :: u_flag       !风速差的方向
                real(r8)        :: delta_u      !u风速差
                real(r8)        :: delta_t      !温度差
                real(r8)        :: delta_q      !比湿差
                
                real(r8)        :: u_star      !摩擦速度
                real(r8)        :: t_star      !温度摩擦速度
                real(r8)        :: q_star      !比湿摩擦速度


                real(r8)        :: grv         !重力加速度
                
                real(r8)        :: beta        !\beta计算粗糙度理查森数系数
                
                real(r8)        :: ribu         !粗糙度理查森数
                
                real(r8)        :: ribcu         !比湿粗糙度
                real(r8)        :: pblh         !边界层高度
                real(r8)        :: nits         !循环次数
                real(r8)        :: zet         !稳定度参数 zet = z/L

                real(r8)        :: cc           !
                real(r8)        :: rr           !粗糙度雷诺数
                real(r8)        :: obklen        ! 奥布霍夫长度
                real(r8)        :: psi_u         !动量稳定度修正函数
                real(r8)        :: psi_t         !温度稳定度修正函数
                real(r8)        :: psi_q         !比湿稳定度修正函数

                real(r8)        :: t_v          ! 虚温
               
                real(r8)        :: t_v_star     ! 虚温特征量



        !输出变量

                real(r8),intent(out)        :: tau      !应力通量
                
                real(r8),intent(out)        :: shf      !感热通量
        
                real(r8),intent(out)        :: lhf      !潜热通量

!==========================================================================
                                     !计算过程
!==========================================================================                
!先定义常参数
                von = 0.4
                visa = 1.326e-5
                pblh = 600 !先定义一个数，后面再进行讨论
                call cal_grv(lat,grv)
!随风速增大的 charnock 数
                charn=0.011 
                if (u > 10) then
                  charn=0.011+(u-10)/(18-10)*(0.018-0.011)
                endif 
                if (u > 18) then
                  charn=0.018 
                endif 
! 由于一些时候，风速的差值很小，可能为0
! ，这样会导致计算不稳定，规定存在风速差最小值

                u_min  = 0.2    !这个数值可以调整，但是不能很大

!判断方向
               u_flag = (u- u_s) / abs( u-u_s )

!计算差值
                delta_u = sqrt( (u-u_s)**2  + u_min**2 )
                
!coare 的问题之一是：
!在处理速度的时候没有考虑方向，因此不能够计算出负的应力通量，之前的理论在计算陆地的通量是没问题的，因为陆地上速度总是为0，这样都是动量被地面消耗，但是海洋上，海流也是存在速度的，这样不分辨方向是致命的错误，可能是海洋向大气输送纬向动量。因此，在最后计算方向的时候乘上代表方向的单位量(u_flag = -1 或者 +1)! 因此加上方向  

!关于计算，计算只能是按照顺序的，前面不存在的变量将无法计算，如果仍需要计算，我们可以先预估一个数字，然后在后面不断逼近真值               

                delta_t = t - t_s 

                delta_q = q - q_s
                
!预估摩擦速度

                u_star = 0.035 *  delta_u * log(10/z)
 
!计算粗糙度z 
                z0  = charn * u_star / grv + 0.11 * visa / u_star
                z0t = 1.15e-4
                z0q = 1.15e-4
!先计算中性条件下 u_*, thv_*, q_*

               u_star   =( delta_u * von  )/(   log(z/z0)   )
               t_star   =( delta_t * von  )/(   log(z/z0t)  )  
               q_star   =( delta_q * von  )/(   log(z/z0q)  )

! 计算 L 
!没有通量没办法计算L， 只能通过理查森数来预估L
                beta = 1.2

               Ribcu=-z/pblh/.004/beta**3 
               Ribu=-grv*z/t*(delta_t+.61*t*delta_q)/delta_t**2 
              
            if (Ribu < 0) then 
                zet=CC*Ribu/(1+Ribu/Ribcu) 
            else 
                zet=CC*Ribu*(1+27/9*Ribu/CC)
            endif 

                obklen = z/zet

!定义迭代次数(如果非常稳定，就不迭代)
                if (zet>50)then
                        nits = 1
                elseif(zet<=50)then
                        nits = 5
                endif

!计算非中性条件下的u_*,t_*,q_*
                call cal_psi_u(zet,psi_u)
                call cal_psi_t(zet,psi_t)
                call cal_psi_q(zet,psi_q)

               u_star   =( delta_u * von  )/(   log(z/z0) -psi_u )
               t_star   =( delta_t * von  )/(   log(z/z0t)-psi_t )
               q_star   =( delta_q * von  )/(   log(z/z0q)-psi_q )

                do i = 1, nits 

!首先计算奥布霍夫长度,稳定度参数                        
                        t_v = t * ( 1 + 0.61 * q )  !单位再确认下

                        t_v_star = t_star*(1+0.61*q)+0.61*t*q_star

                        obklen = u_star**2 / ( von * (grv/t_v) * t_v_star )
                        zet = z / obklen 
 
!计算表面粗糙度
                        rr = z0*u_star/visa
                        z0q = min(1.15e-4,5.5e-5/rr**0.6 )
                        z0t = 1.15e-4

!计算u_*，t_*,q_*
                        call cal_psi_u(zet,psi_u)
                        call cal_psi_t(zet,psi_t)
                        call cal_psi_q(zet,psi_q)

                        u_star   =( delta_u * von  )/(   log(z/z0) -psi_u )
                        t_star   =( delta_t * von  )/(   log(z/z0t)-psi_t )
                        q_star   =( delta_q * von  )/(   log(z/z0q)-psi_q )
                end do

                tau =   rho_a * u_star * u_star * u_flag
                
                shf = - rho_a * u_star * t_star 

                lhf = - rho_a * u_star * q_star 
                       
                
        contains

        end subroutine cal_flux_iap 

      end module flux_iap

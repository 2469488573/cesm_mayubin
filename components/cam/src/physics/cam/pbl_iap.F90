        module pbl_iap

        !=========================================================
        !这个module是大气所的边界层参数化方案，主要计算边界层    !
        !表面层通量，                                            !
        !边界层高度，                                            !
        !云量和雾，                                              !
        !其中深度学习的方法被应用                                !
        !                                                        !
        !作者： 马钰斌、程雪玲、吴琳、靳江波、李奇龙、曾庆存     !
        !创建时间： 2022年12月8日                                !
        !=========================================================


        use shr_kind_mod, only: r8 => shr_kind_r8

        implicit none

        private

        public calc_pbl_h_vector
        public calc_pbl_h_vector_no_iterate
        public get_biaozhungaodu
        public calc_pbl_h_deeplearn
        public yunliang
        public calc_pblh_tbf

        contains

      subroutine   get_biaozhungaodu(ncol,pver,p,h_pver)
         ! get mei ge mo shi ceng de biao zhun de gao du (unit:m) 

         integer,intent(in)    ::  pver
         integer,intent(in)    ::  ncol
         real(r8),intent(in)   ::  p(ncol,pver)  ! unit : hPa
         integer               ::  i , j
         real(r8)              ::  h_pver_sum 
         real(r8)              ::  h_ncol_pver(ncol,pver)
         real(r8),intent(out)  ::  h_pver(pver)

         h_ncol_pver = 44300_r8 * ( 1_r8 - ( p/1013.25_r8 )**( 1_r8/5.256_r8 ) )

        !dui ncol jinxing pingjun ,dedao h_pver(pver)
         do j = 1,pver
                 h_pver_sum = 0
                 do i = 1,ncol 
                         h_pver_sum = h_pver_sum + h_ncol_pver(i,j)
                 end do 
                 h_pver(j) = h_pver_sum / ncol
         end do 
      end subroutine get_biaozhungaodu
      
        subroutine calc_pbl_h_vector(h_ncol_pver,pbl_h_old,obklen,latvals,pver,n,T_c,tau_star,hsb_star,g_in,n2,PBL_H)
        !--------------------------------------!
        !purpose:calculate the height of PBL   !
        !mothod from the zili(2005)  paper     !
        !                                      !
        !mayubin 2022-3-25                     !
        !--------------------------------------!

        !=================!
        !  input variable !
        !=================!
        
        integer, intent(in)  :: n
        integer, intent(in)  :: pver
        real(r8),intent(in)  :: pbl_h_old(n)
        real(r8),intent(in)  :: obklen(n)
        real(r8),intent(in)  :: g_in  
        real(r8),intent(in)  :: T_c(n)
        real(r8),intent(in)  :: tau_star(n) 
        real(r8),intent(in)  :: hsb_star(n)      
        real(r8),intent(in)  :: n2(n,pver)
        real(r8),intent(in)  :: latvals(n)     !lat in radians (rad)
        real(r8),intent(in)  :: h_ncol_pver(n,pver)
        !=================!
        !  local variable !
        !=================!

        real(r8)             :: C_R =0.6_r8
        real(r8)             :: C_CN =1.36_r8
        real(r8)             :: C_NS =0.51_r8
        real(r8)             :: NN=0.01_r8     !ping lv
        real(r8)             :: f(n)           !defult 0.00005_r8 
        real(r8)             :: h_1(n)
        real(r8)             :: h_2(n)
        real(r8)             :: h_3(n)
        real(r8)             :: beta_b(n)
        real(r8)             :: h_pver_sum 
        integer              :: i
        integer              :: j
        integer              :: k
        real(r8)             :: n2_vertical_mean(n)
        real(r8)             :: pbl_h_stable(n)
        integer              :: diedaiindex =2
        integer,parameter    :: zuidadiedaicishu = 10 !die dai zui da ci shu 
        real(r8)             :: pbl_h_array(n,zuidadiedaicishu)
        integer              :: pingjuncengshu = 1
        real(r8)             :: h_pver(pver)
        integer              :: deltamin,deltamin_index
        real(r8)             :: delta(pver) 
        real(r8)             :: hsb_star_nocp(n)
        real(r8)             :: h_2_temp
        real(r8)             :: pbl_h_init(n)
        real(r8)             :: tau_star_2(n)
        !==================!
        !  output variable !
        !==================!

        real(r8),intent(out) :: PBL_H(n)

        !===================!
        ! ji suan guo cheng !
        !===================!

!!!!!====zheli xian ceshi yixia nayi ceng suanchulaide zuizhengchang ==
        ! brunt vaisaila frequency chuizhifangxiang pingjun    
        do i=1,n
                n2_vertical_mean(i) = 0
                do k =1,pver
                        n2_vertical_mean(i) = n2_vertical_mean(i)+n2(i,k)
                end do 
                n2_vertical_mean(i) = n2_vertical_mean(i)/pver
        end do 


        f=2_r8*7.292_r8*0.00001_r8*sin(latvals)

        beta_b = g_in/T_c

  
        ! hsb_star  ke neng duo cheng le yi ge C_p ,yao bu yao chu diao ,zheyang hui
        ! xian zhu zeng jia bian jie ceng gao du  
       !tau_star he u_star

        tau_star_2 = tau_star**2

        hsb_star_nocp  = hsb_star 

        h_1 = (f**2)/((C_R**2)*tau_star)

        h_2 = (sqrt(abs(n2_vertical_mean))*abs(f))/(tau_star*(C_CN**2))

        h_3 = abs(f*hsb_star_nocp*beta_b)/(C_CN**2*tau_star**2)

        pbl_h_stable = sqrt(1/(h_1+h_2+h_3))
        pbl_h_init = pbl_h_stable

        !you yu dui fu li ping lv ji fen shi bu zhi dao gao du ,yin ci xu yao die dai
        !ji suan 
        
!        pbl_h_array(:,1)= pbl_h_stable 

        !dui ncols jin xin pingjun deidao h_pver(pver)
!         do j = 1,pver
!                 h_pver_sum = 0
!                 do i = 1,n 
!                         h_pver_sum = h_pver_sum + h_ncol_pver(i,j)
!                 end do 
!                 h_pver(j) = h_pver_sum / n
!                 print*,h_pver(j)
!         end do 

!        do i = 1,pver 
!                print*, "di" ,i , "ceng,biaohzungaodu=" ,h_pver(i)
!        end do 


        
!        do i = 1,pver 
!                do j = 1,n
!                   print*," i = ",i,", j = ",j,", biaozhungaodu =",h_ncol_pver(j,i)
!                end do 
!        end do 
      do k = 1,n
      diedaiindex = 2
      pbl_h_array(k,1)= pbl_h_stable(k) 
        do while( abs(pbl_h_array(k,diedaiindex)-pbl_h_array(k,diedaiindex-1)) > 100&
                 .and. diedaiindex<zuidadiedaicishu )
                !pan duan pbl_h he na ge pver zui jie jin 
                do i = 1 , pver
                        delta(i) = abs(pbl_h_array(k,diedaiindex) - h_ncol_pver(k,i) )
                end do 
                
                !jisuan zuixiaozhi bing jilu  index 
                deltamin = delta(1)
                deltamin_index = 1

                do i = 1,pver
                        if( delta(i)<deltamin) then 
                                deltamin = delta(i)
                                deltamin_index = i 
                        end if 

                end do 
                         
                !jisuan zuijiejin de na yi ceng yi ji shangmian ji ceng N_mean

                
                       ! n2_vertical_mean(i) = 0_r8

                      !  do j =deltamin_index-pingjuncengshu,deltamin_index
                      !          n2_vertical_mean(i) = n2_vertical_mean(i)+&
                      !          n2(i,j)
                      !  end do 
                      !  n2_vertical_mean(i) = n2_vertical_mean(i)/(pingjuncengshu+1)
                        n2_vertical_mean(k)=n2(k,deltamin_index-1)

                !ji suan xin de bian jie ceng gao du  
                h_2_temp = (sqrt(abs(n2_vertical_mean(k)))*abs(f(k)))/(tau_star(k)*(C_CN**2))
                
                pbl_h_array(k,diedaiindex) = sqrt(1/(h_1(k)+h_2_temp+h_3(k)))


                !index j add 1 
                diedaiindex=diedaiindex+1
        end do 
          pbl_h_stable(k) = pbl_h_array(k,diedaiindex-1)

     end do 


        !dui wen ding du jin xing pan duan 

        do  i = 1, n
                if ( obklen(i) > 0 ) then 
                        PBL_H(i) = pbl_h_stable(i) 
                else
                        PBL_H(i) = pbl_h_stable(i)
                end if
        end do
 


       !!!test

        !===============================
        !PBL_H = tau_star**2/g_in
        !===============================
        !PBL_H = g_in/(f**2)
        !===============================
        !PBL_H = abs(tau_star)/f
        !===============================
        !PBL_H = ((tau_star**4)*(f**2))/(g_in**3)
        !===============================
        !PBL_H = g_in/n2_vertical_mean   
        !===============================
        !PBL_H = (g_in*tau_star**2)**3  
        !===============================   
        !PBL_H = T_c/hsb_star 
        !===============================

      end subroutine calc_pbl_h_vector

 subroutine calc_pbl_h_vector_no_iterate(ri,h_ncol_pver,pbl_h_old,obklen,latvals,pver,n,T_c,tau_star,hsb_star,g_in,n2,PBL_H)
        !--------------------------------------!
        !purpose:calculate the height of PBL   !
        !mothod from the zili(2005)  paper     !
        !                                      !
        !mayubin 2022-5-19                     !
        !--------------------------------------!

        !=================!
        !  input variable !
        !=================!
        
        integer, intent(in)  :: n
        integer, intent(in)  :: pver
        real(r8),intent(in)  :: pbl_h_old(n)
        real(r8),intent(in)  :: obklen(n)
        real(r8),intent(in)  :: g_in  
        real(r8),intent(in)  :: T_c(n)
        real(r8),intent(in)  :: tau_star(n) 
        real(r8),intent(in)  :: hsb_star(n)      
        real(r8),intent(in)  :: n2(n,pver),ri(n,pver)
        real(r8),intent(in)  :: latvals(n)     !lat in radians (rad)
        real(r8),intent(in)  :: h_ncol_pver(n,pver)
        !=================!
        !  local variable !
        !=================!

        real(r8)             :: C_R =0.6_r8
        real(r8)             :: C_CN =1.36_r8
        real(r8)             :: C_NS =0.51_r8
        real(r8)             :: NN=0.01_r8     !ping lv
        real(r8)             :: f(n)           !defult 0.00005_r8 
        real(r8)             :: h_1(n)
        real(r8)             :: h_2(n)
        real(r8)             :: h_3(n)
        real(r8)             :: beta_b(n)
        real(r8)             :: h_pver_sum 
        integer              :: i
        integer              :: j
        integer              :: k
        real(r8)             :: n2_vertical_mean(n)
        real(r8)             :: pbl_h_stable(n)
        integer              :: diedaiindex =2
        integer,parameter    :: zuidadiedaicishu = 10 !die dai zui da ci shu 
        real(r8)             :: pbl_h_array(n,zuidadiedaicishu)
        integer              :: pingjuncengshu = 4
        real(r8)             :: h_pver(pver)
        integer              :: deltamin,deltamin_index
        real(r8)             :: delta(pver) 
        real(r8)             :: hsb_star_nocp(n)
        real(r8)             :: h_2_temp
        real(r8)             :: pbl_h_init(n)
        real(r8)             :: tau_star_2(n)
        real(r8)             :: f_min
        !==================!
        !  output variable !
        !==================!

        real(r8),intent(out) :: PBL_H(n)

        !===================!
        ! ji suan guo cheng !
        !===================!

        !!!!!====zheli xian ceshi yixia nayi ceng suanchulaide zuizhengchang ==
        ! brunt vaisaila frequency chuizhifangxiang pingjun    
        do i=1,n
                n2_vertical_mean(i) = 0
                do k =24,28
                        n2_vertical_mean(i) = n2_vertical_mean(i)+n2(i,k)
                end do 
                n2_vertical_mean(i) = n2_vertical_mean(i)/5
        end do 
do i = 1,n

        if(latvals(i)<0.and.latvals(i)>-15*3.14/180)then
                f(i) = -15*3.14/180
        end if
        if(latvals(i) >= 0 .and. latvals(i) <15* 3.14/180)then
                f(i) = 15*3.14/180
        end if

end do


        f=2_r8*7.292_r8*0.00001_r8*sin(latvals)
        
!        do i = 1,n
!                if(f(i)<f_min)then
!                        f(i) = f_min
!                end if
!        end do

        beta_b = g_in/T_c
  
        ! hsb_star  ke neng duo cheng le yi ge C_p ,yao bu yao chu diao ,zheyang hui
        ! xian zhu zeng jia bian jie ceng gao du  
        !tau_star he u_star
        hsb_star_nocp  = hsb_star 

        h_1 = (f**2)/((C_R**2)*tau_star)
        h_2 = (sqrt(abs(n2_vertical_mean))*abs(f))/(tau_star*(C_CN**2))
        h_3 = abs(f*hsb_star_nocp*beta_b)/(C_CN**2*tau_star**2)
        pbl_h_stable = sqrt(1/(h_1+h_2+h_3))

!dui xinfangan de maxmin h jinxing xianzhi 

        do i = 1,n
                if(pbl_h_stable(i)>3000)then
                        pbl_h_stable(i) = 3000
                end if
                if(pbl_h_stable(i)<50)then 
                        pbl_h_stable = 50
                end if  
        end do
        !chakan ri
!        do i = 1,n
!                do j = 1,pver
!                        print*,i,j,ri(i,j)
!                end do
!        end do 
        !dui wen ding du jin xing pan duan 

        do  i = 1, n
                
                if ( obklen(i) > 0 .and. ri(i,31)>0.25_r8) then
!                        print*,"di",i,"ncol" ,ri(i,31) 
                        PBL_H(i) = pbl_h_stable(i) 
                else
                        PBL_H(i) = pbl_h_old(i)
                end if
        end do

      end subroutine calc_pbl_h_vector_no_iterate



      subroutine calc_flux()
        !================mayubin==============
        !
        !purpose : cal flux of surface 
        !moment flux 
        !sensible heat flux 
        !latent heat flux 
        !=================================





      end subroutine calc_flux
     
       subroutine calc_pbl_h_deeplearn_three_factor( ncol, tau ,shf ,lhf ,pbl_h_dp )
        !mayubin 2022/6/24 
        !input
        integer, intent(in)   :: ncol 
        real(r8),intent(in)   :: tau(ncol)
        real(r8),intent(in)   :: shf(ncol)
        real(r8),intent(in)   :: lhf(ncol)
        
        !local 
        real(r8)             :: weight00(64,3)
        real(r8)             :: weight02(64,1)
        real(r8)             :: bias00(64,1)
        real(r8)             :: bias02 = 7.667402
        integer              :: i,j
        real(r8)             :: wxb1(64,ncol) , wx2(64,ncol), reluwxb1(64,ncol)
        real(r8)             :: wx2sum
        real(r8)             :: pingjun(ncol),biaozhuncha(ncol),tau_biaozhunhua(ncol),&
                                shf_biaozhunhua(ncol),lhf_biaozhunhua(ncol),shf_w(ncol),lhf_w(ncol),&
                                tau_std(ncol),shf_std(ncol),lhf_std(ncol)
        real(r8) tau_ncol_sum,shf_ncol_sum,lhf_ncol_sum,&
                 tau_pingjun,shf_pingjun,lhf_pingjun,&
                 tau_std_sum,shf_std_sum,lhf_std_sum,&
                 tau_std_sum_mean,shf_std_sum_mean,lhf_std_sum_mean




        !output 
        real(r8),intent(out)  :: pbl_h_dp(ncol)
        ! danweizhuanhua J he W.h
        lhf_w = lhf
        shf_w = shf
        !nn output parameter are constant 
        data weight00/1.073475,-0.544754,0.750596,-0.264516,-0.490043,-0.467062,1.462132,1.504056,-0.621389,0.187883,0.395087,-0.833701,0.615115,-0.581890,-0.746298,0.830845,-0.777304,-0.497521,0.542383,0.305653,-1.905314,0.473748,1.681869,-0.724897,1.840948,-0.600823,-0.026229,-0.698829,0.095268,0.342287,0.799918,-0.951834,-0.614076,0.509474,-0.910255,0.870776,-0.757339,-0.422914,2.024335,2.005388,0.701373,-0.439053,0.035714,0.193852,-0.063280,0.001888,-0.700263,-0.637911,0.345497,0.290306,1.742877,-0.564842,-0.821523,0.428326,-0.639363,0.740396,0.836748,0.529985,-0.677184,0.586683,1.516016,-0.055916,-1.665320,-0.708855,0.532824,0.468415,0.068437,0.895945,0.261880,0.585849,1.147956,0.502968,0.685725,-0.484371,-0.081515,-1.602120,-0.041258,0.504835,-0.456619,0.391319,-0.305853,0.301241,-0.576034,-0.262759,-2.097549,0.856788,-0.011481,0.010287,1.098459,0.832180,0.068230,-0.448985,-0.054990,-0.071903,1.489364,-1.968611,0.630247,-0.326224,-0.342264,1.671796,-0.355720,0.203099,0.450155,1.085203,-0.263546,-0.725751,0.780647,-0.608344,0.944698,-0.014309,0.260791,-0.477309,-0.071403,-0.450961,1.145328,0.718276,-0.405357,-0.764062,0.346342,-0.534813,-0.263367,0.319319,0.459508,1.266264,0.881500,0.213432,-1.878592,-0.085977,0.177805,0.666727,0.920415,0.871609,0.619972,0.341816,0.314589,-1.052345,0.183469,-0.492031,-0.719630,1.013233,-0.526999,0.788828,-0.267416,0.573676,-0.351517,0.574789,-0.057034,-0.665392,1.111592,0.022094,-0.194048,-0.459517,-0.533533,-0.029879,-0.070508,-0.306848,-0.199895,-0.628550,0.271560,1.242666,0.173730,-0.551856,-0.206569,0.262017,-0.033655,0.786774,-0.208505,-0.452115,-0.422749,-0.489955,0.816471,-0.158786,0.758739,0.007972,0.755914,-0.136933,-0.736746,-0.225053,-0.456523,0.524465,-0.041105,-0.066913,0.446612,-0.044642,-0.299893,0.623790,0.341658,0.804704,-1.227542,-0.488892,0.796174,0.888756/

       data weight02/2.028993,1.997387,2.081208,1.830254,2.096154,2.223962,1.325882,2.358625,1.991714,2.090166,2.198164,0.972804,2.254223,1.846843,2.006792,1.977402,2.056325,2.222537,2.309003,2.042206,1.332056,2.102832,1.832685,2.239376,1.590063,2.251463,-0.026031,2.073334,-0.063853,2.283603,1.799824,1.103766,2.086398,1.980928,2.020166,1.758786,2.485721,2.154185,1.611543,1.472993,2.070014,2.021032,1.959863,2.319531,1.886000,-0.025161,1.912920,2.407532,2.125054,2.496415,1.454528,1.873072,2.294024,2.120467,2.067123,2.171715,2.206501,2.228948,2.056305,2.168296,2.366634,-0.055226,0.526147,1.970670/
       data bias00/2.903134,2.941154,2.885986,2.965427,2.522568,2.511720,1.871562,2.403379,2.918393,3.084845,2.896975,0.266564,2.812622,3.211048,3.085637,3.106334,3.143068,2.541197,2.724314,3.253774,0.088564,2.857390,2.215092,2.855591,2.032742,2.802315,-0.034591,3.009327,-0.556008,2.537859,2.149551,0.402572,2.822366,3.043907,3.314974,2.360861,2.658409,2.432288,1.876585,1.833542,3.109980,2.874225,2.509058,2.641164,2.729093,-0.013446,3.276453,2.619589,2.999193,2.524535,2.235947,3.114990,2.891423,3.055188,2.875186,2.865269,2.694156,2.524986,3.006970,2.793157,2.107884,-0.523572,0.485908,3.076036/
       bias02 = 1.493356 

       !biaozhunhua
        
        tau_ncol_sum = 0.0
        shf_ncol_sum = 0.0
        lhf_ncol_sum = 0.0

     do j =1,ncol
        tau_ncol_sum = tau_ncol_sum + tau(j)
        shf_ncol_sum = shf_ncol_sum + shf(j)
        lhf_ncol_sum = lhf_ncol_sum + lhf(j)
     end do
        tau_pingjun  = tau_ncol_sum/ncol
        shf_pingjun  = shf_ncol_sum/ncol
        lhf_pingjun  = lhf_ncol_sum/ncol
        
        tau_std_sum  = 0.0
        shf_std_sum  = 0.0 
        lhf_std_sum  = 0.0

        do j = 1,ncol
          tau_std(j) = (tau(j) - tau_pingjun)**2
          tau_std_sum = tau_std_sum + tau_std(j)
          shf_std(j) = (shf(j)-shf_pingjun)**2
          shf_std_sum = shf_std_sum + shf_std(j)
          lhf_std(j) =(lhf(j)-lhf_pingjun)**2
          lhf_std_sum = shf_std_sum + shf_std(j)
                
        end do 
        tau_std_sum_mean = tau_std_sum/ncol
        shf_std_sum_mean = shf_std_sum/ncol
        lhf_std_sum_mean = lhf_std_sum/ncol 

        do j = 1,ncol 
        tau_biaozhunhua(j) = (tau(j) - tau_pingjun)/ tau_std_sum_mean
        shf_biaozhunhua(j) = (shf(j) - shf_pingjun)/ shf_std_sum_mean
        lhf_biaozhunhua(j) = (lhf(j) - lhf_pingjun)/ lhf_std_sum_mean
        enddo
        
     do j  = 1,ncol
        do i = 1,64
                wxb1(i,j) = bias00(i,1)+weight00(i,1)*tau_biaozhunhua(j)+weight00(i,2)*lhf_biaozhunhua(j)&
                +weight00(i,3)*shf_biaozhunhua(j)
        end do            

        do i =1,64

                if(wxb1(i,j)<0) then
                        reluwxb1(i,j) = 0
                else 
                        reluwxb1(i,j) = wxb1(i,j)
                end if 
        end do

        do i =1,64
                wx2(i,j) = weight02(i,1)*reluwxb1(i,j)
        end do 
        wx2sum = 0
        do i = 1,64
                wx2sum = wx2sum + wx2(i,j)
        end do
        pbl_h_dp(j) = wx2sum + bias02 
     end do
        contains
        subroutine  biaozhunhua(x,n,x_mean,x_std_sum_mean,x_biaozhunhua)
          real,intent(in)    :: x(n)
          integer,intent(in) :: n
          real,intent(out)   :: x_mean,x_std_sum_mean
          real(r8)            :: x_sum,x_std_sum 
          real(r8)            :: x_std(n)
          real(r8),intent(out) :: x_biaozhunhua(n)

          integer            :: i 
          x_sum = 0.0
          do i =1,n
            x_sum = x_sum + x(i)
          enddo 
          x_mean = x_sum/n

          x_std = 0.0
          do i = 1,n
            x_std(i) = x(i) - x_mean
          end do 
          x_std_sum = 0.0
          do i = 1,n
            x_std_sum = x_std_sum + x_std(i)
          enddo 
          x_std_sum_mean = x_std_sum/n
          do i = 1,n
            x_biaozhunhua(i) =( x(i) -x_mean)/x_std_sum_mean
          end do 
        end subroutine biaozhunhua
      end subroutine calc_pbl_h_deeplearn_three_factor
      
      subroutine calc_pblh_tbf(ncol,t_c,tau,shf,lhf,pblh_tbf)
!torch_bridge_fortran 被写成了module，方便在其他程序中调用
!其中使用者可以把它当成黑箱，深度学习训练好的模型被封装在其中
!需要输入的是一个二维数组，x_array(cesm_m,ncol),cesm_m 是因子个数
!ncol和cesm的ncol是一致的。输出的y_cesm就是深度模型的计算量
!x1, i --|              ______
!x2, i --|             |     |     
!x3, i --|---------->  | box | ------->  y_cesm(i)
!x4, i --|             |_____|
!x5, i --|              

!这是cesm中调用TBF的一个实例程序，在bridge程序中，计算被进行
!参数传递的过程在bridge进行，
!参数优化的过程则是在python文件夹下进行
!计算过程则主要写在calculation.f90中，这一部分不建议修改，需要对深度模型很了解
!文件读取在file_io.f90中进行
!----------------------------------------------------------------

!在fortran模型代码中，我们主要需要做下面几件事


!1.use module 

        use bridge , only: tbf
        use shr_kind_mod , only:r8 => shr_kind_r8  

!-----------------------------------------------------------------

!2.申明变量和变量的维度

        implicit none 
        integer, intent(in)   :: ncol 
        real(r8),intent(in)   :: tau(ncol)
        real(r8),intent(in)   :: shf(ncol)
        real(r8),intent(in)   :: lhf(ncol)
        
        real(r8),intent(in)   :: t_c(ncol)
        !local 
        integer              :: i,j
 
        integer              :: m = 4
        real(r8),allocatable :: x_array(:,:)
        real(r8)             :: pblh_tbf(ncol)


!-----------------------------------------------------------------


!3.明确自变量（可以通过其他函数传递）!如果多个变量向量则需要进行数组拼接

        allocate(x_array(m,ncol))


                x_array(1,:) = t_c
                x_array(2,:) = tau
                x_array(3,:) = shf
                x_array(4,:) = lhf            


!------------------------------------------------------------------

!4.调用子程序tbf( 数组长度，因子个数，自变量数组, 预报量数组)

        call tbf(ncol,m,x_array,pblh_tbf)

!------------------------------------------------------------------

!5.检查和传递计算结果

        print*,pblh_tbf 

      end subroutine calc_pblh_tbf


      subroutine calc_pbl_h_deeplearn( ncol, tau ,shf ,lhf ,pbl_h_dp )
        !mayubin 2022/6/24 
        !input
        integer, intent(in)   :: ncol 
        real(r8),intent(in)   :: tau(ncol)
        real(r8),intent(in)   :: shf(ncol)
        real(r8),intent(in)   :: lhf(ncol)
        
        !local 
        real(r8)             :: weight00(64,3)
        real(r8)             :: weight02(64,1)
        real(r8)             :: bias00(64,1)
        real(r8)             :: bias02 = 7.667402
        integer              :: i,j
        real(r8)             :: wxb1(64,ncol) , wx2(64,ncol), reluwxb1(64,ncol)
        real(r8)             :: wx2sum
        real(r8)             :: pingjun(ncol),biaozhuncha(ncol),tau_biaozhunhua(ncol),&
                                shf_biaozhunhua(ncol),lhf_biaozhunhua(ncol),shf_w(ncol),lhf_w(ncol),&
                                tau_std(ncol),shf_std(ncol),lhf_std(ncol)
        real(r8) tau_ncol_sum,shf_ncol_sum,lhf_ncol_sum,&
                 tau_pingjun,shf_pingjun,lhf_pingjun,&
                 tau_std_sum,shf_std_sum,lhf_std_sum,&
                 tau_std_sum_mean,shf_std_sum_mean,lhf_std_sum_mean




        !output 
        real(r8),intent(out)  :: pbl_h_dp(ncol)
        ! danweizhuanhua J he W.h
        lhf_w = lhf
        shf_w = shf
        !nn output parameter are constant 
        data weight00/1.073475,-0.544754,0.750596,-0.264516,-0.490043,-0.467062,1.462132,1.504056,-0.621389,0.187883,0.395087,-0.833701,0.615115,-0.581890,-0.746298,0.830845,-0.777304,-0.497521,0.542383,0.305653,-1.905314,0.473748,1.681869,-0.724897,1.840948,-0.600823,-0.026229,-0.698829,0.095268,0.342287,0.799918,-0.951834,-0.614076,0.509474,-0.910255,0.870776,-0.757339,-0.422914,2.024335,2.005388,0.701373,-0.439053,0.035714,0.193852,-0.063280,0.001888,-0.700263,-0.637911,0.345497,0.290306,1.742877,-0.564842,-0.821523,0.428326,-0.639363,0.740396,0.836748,0.529985,-0.677184,0.586683,1.516016,-0.055916,-1.665320,-0.708855,0.532824,0.468415,0.068437,0.895945,0.261880,0.585849,1.147956,0.502968,0.685725,-0.484371,-0.081515,-1.602120,-0.041258,0.504835,-0.456619,0.391319,-0.305853,0.301241,-0.576034,-0.262759,-2.097549,0.856788,-0.011481,0.010287,1.098459,0.832180,0.068230,-0.448985,-0.054990,-0.071903,1.489364,-1.968611,0.630247,-0.326224,-0.342264,1.671796,-0.355720,0.203099,0.450155,1.085203,-0.263546,-0.725751,0.780647,-0.608344,0.944698,-0.014309,0.260791,-0.477309,-0.071403,-0.450961,1.145328,0.718276,-0.405357,-0.764062,0.346342,-0.534813,-0.263367,0.319319,0.459508,1.266264,0.881500,0.213432,-1.878592,-0.085977,0.177805,0.666727,0.920415,0.871609,0.619972,0.341816,0.314589,-1.052345,0.183469,-0.492031,-0.719630,1.013233,-0.526999,0.788828,-0.267416,0.573676,-0.351517,0.574789,-0.057034,-0.665392,1.111592,0.022094,-0.194048,-0.459517,-0.533533,-0.029879,-0.070508,-0.306848,-0.199895,-0.628550,0.271560,1.242666,0.173730,-0.551856,-0.206569,0.262017,-0.033655,0.786774,-0.208505,-0.452115,-0.422749,-0.489955,0.816471,-0.158786,0.758739,0.007972,0.755914,-0.136933,-0.736746,-0.225053,-0.456523,0.524465,-0.041105,-0.066913,0.446612,-0.044642,-0.299893,0.623790,0.341658,0.804704,-1.227542,-0.488892,0.796174,0.888756/

       data weight02/2.028993,1.997387,2.081208,1.830254,2.096154,2.223962,1.325882,2.358625,1.991714,2.090166,2.198164,0.972804,2.254223,1.846843,2.006792,1.977402,2.056325,2.222537,2.309003,2.042206,1.332056,2.102832,1.832685,2.239376,1.590063,2.251463,-0.026031,2.073334,-0.063853,2.283603,1.799824,1.103766,2.086398,1.980928,2.020166,1.758786,2.485721,2.154185,1.611543,1.472993,2.070014,2.021032,1.959863,2.319531,1.886000,-0.025161,1.912920,2.407532,2.125054,2.496415,1.454528,1.873072,2.294024,2.120467,2.067123,2.171715,2.206501,2.228948,2.056305,2.168296,2.366634,-0.055226,0.526147,1.970670/
       data bias00/2.903134,2.941154,2.885986,2.965427,2.522568,2.511720,1.871562,2.403379,2.918393,3.084845,2.896975,0.266564,2.812622,3.211048,3.085637,3.106334,3.143068,2.541197,2.724314,3.253774,0.088564,2.857390,2.215092,2.855591,2.032742,2.802315,-0.034591,3.009327,-0.556008,2.537859,2.149551,0.402572,2.822366,3.043907,3.314974,2.360861,2.658409,2.432288,1.876585,1.833542,3.109980,2.874225,2.509058,2.641164,2.729093,-0.013446,3.276453,2.619589,2.999193,2.524535,2.235947,3.114990,2.891423,3.055188,2.875186,2.865269,2.694156,2.524986,3.006970,2.793157,2.107884,-0.523572,0.485908,3.076036/
       bias02 = 1.493356 

       !biaozhunhua
        
        tau_ncol_sum = 0.0
        shf_ncol_sum = 0.0
        lhf_ncol_sum = 0.0

     do j =1,ncol
        tau_ncol_sum = tau_ncol_sum + tau(j)
        shf_ncol_sum = shf_ncol_sum + shf(j)
        lhf_ncol_sum = lhf_ncol_sum + lhf(j)
     end do
        tau_pingjun  = tau_ncol_sum/ncol
        shf_pingjun  = shf_ncol_sum/ncol
        lhf_pingjun  = lhf_ncol_sum/ncol
        
        tau_std_sum  = 0.0
        shf_std_sum  = 0.0 
        lhf_std_sum  = 0.0

        do j = 1,ncol
          tau_std(j) = (tau(j) - tau_pingjun)**2
          tau_std_sum = tau_std_sum + tau_std(j)
          shf_std(j) = (shf(j)-shf_pingjun)**2
          shf_std_sum = shf_std_sum + shf_std(j)
          lhf_std(j) =(lhf(j)-lhf_pingjun)**2
          lhf_std_sum = shf_std_sum + shf_std(j)
                
        end do 
        tau_std_sum_mean = tau_std_sum/ncol
        shf_std_sum_mean = shf_std_sum/ncol
        lhf_std_sum_mean = lhf_std_sum/ncol 

        do j = 1,ncol 
        tau_biaozhunhua(j) = (tau(j) - tau_pingjun)/ tau_std_sum_mean
        shf_biaozhunhua(j) = (shf(j) - shf_pingjun)/ shf_std_sum_mean
        lhf_biaozhunhua(j) = (lhf(j) - lhf_pingjun)/ lhf_std_sum_mean
        enddo
        
     do j  = 1,ncol
        do i = 1,64
                wxb1(i,j) = bias00(i,1)+weight00(i,1)*tau_biaozhunhua(j)+weight00(i,2)*lhf_biaozhunhua(j)&
                +weight00(i,3)*shf_biaozhunhua(j)
        end do            

        do i =1,64

                if(wxb1(i,j)<0) then
                        reluwxb1(i,j) = 0
                else 
                        reluwxb1(i,j) = wxb1(i,j)
                end if 
        end do

        do i =1,64
                wx2(i,j) = weight02(i,1)*reluwxb1(i,j)
        end do 
        wx2sum = 0
        do i = 1,64
                wx2sum = wx2sum + wx2(i,j)
        end do
        pbl_h_dp(j) = wx2sum + bias02 
     end do
        contains
        subroutine  biaozhunhua(x,n,x_mean,x_std_sum_mean,x_biaozhunhua)
          real,intent(in)    :: x(n)
          integer,intent(in) :: n
          real,intent(out)   :: x_mean,x_std_sum_mean
          real(r8)            :: x_sum,x_std_sum 
          real(r8)            :: x_std(n)
          real(r8),intent(out) :: x_biaozhunhua(n)

          integer            :: i 
          x_sum = 0.0
          do i =1,n
            x_sum = x_sum + x(i)
          enddo 
          x_mean = x_sum/n

          x_std = 0.0
          do i = 1,n
            x_std(i) = x(i) - x_mean
          end do 
          x_std_sum = 0.0
          do i = 1,n
            x_std_sum = x_std_sum + x_std(i)
          enddo 
          x_std_sum_mean = x_std_sum/n
          do i = 1,n
            x_biaozhunhua(i) =( x(i) -x_mean)/x_std_sum_mean
          end do 
        end subroutine biaozhunhua
      end subroutine calc_pbl_h_deeplearn
     
        subroutine yunliang(ncol,obklen , rh , n_l)
                integer,intent(in)  :: ncol
                real(r8)            :: rh(ncol) 
                real(r8),intent(in) :: obklen(ncol)
                integer             :: i
                real(r8)            :: n_l_stable(ncol)
                real(r8)            :: n_l_unstable(ncol)

                real(r8),intent(out) :: n_l(ncol)
                !rh >100 jiu jiashe wei 100 
                do i = 1,ncol 
                        if(rh(i)>100_r8)then 
                                rh(i) = 100_r8
                        end if
                end do 

                !stable 
                do i = 1, ncol 
                        if(rh(i)>= 80_r8 )then 
                                n_l_stable(i)  = ((rh(i) - 80_r8)**2_r8) / 4_r8
                        else if(rh(i)<80_r8)then
                                n_l_stable(i)  = 0_r8
                        end if
                end do


                !unstable 
                do  i = 1,ncol
                        if(rh(i)>=57_r8 )then  
                                 n_l_unstable(i)  =-0.74_r8 + -1.257_r8*rh(i) + 0.022644_r8*rh(i)*rh(i)
                        else if (rh(i)<57_r8)then 
                                n_l_unstable(i)  = 0_r8  
                        end if 
                end do 

                !panduan 
                do i = 1 , ncol
                        if(obklen(i)>=0_r8) then 
                                n_l(i) = n_l_stable(i)
                        else if( obklen(i)< 0_r8 )then 
                                n_l(i) = n_l_unstable(i)
                        end if 
                end do 
        
                !test 
!                do i = 1,ncol
!                        n_l(i) = obklen(i)
!                end do 
                        
        end subroutine yunliang





        end module pbl_iap

        module  physics_iap 
     !这个module主要放物理定律的代码
        use shr_kind_mod , only : r8 => shr_kind_r8
  
        implicit none 
        private 
        public cal_grv 
        public cal_psi_u 
        public cal_psi_t 
        public cal_psi_q 
        public cal_qsat 
        public cal_qsee

        contains        

        subroutine  cal_grv(lat,grv) 
     ! 使用纬度计算当地的重力加速度，相关的研究看“重力加速度和维度的关系”文章
     ! 将重力加速度设为和纬度相关的三角函数
     ! lat 的单位是度,其中w 约等于2*pi/180,grv 的单位是m/s^2

              real(r8),intent(in) :: lat 
              real(r8)   ::  g_0 , a_1 , b_1 , w
              real(r8),intent(out)   :: grv

              g_0  =  9.806
              a_1  = -0.02575 
              b_1  = -0.0008862
              w    =  0.03542
      
              grv = g_0 + a_1 * cos( w * lat ) + b_1 * sin( w * lat )      

      end subroutine cal_grv
      
      
      subroutine cal_psi_u(zet,psi_u)
          real(r8) :: x, psik,psic,f,zet,c,psi_u
          x=(1.-15.*zet)**.25 
          psik=2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.) 
          x=(1.-10.15*zet)**.3333 
          psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
          f=zet*zet/(1+zet*zet) 
          psi_u=(1-f)*psik+f*psic                                                
          if(zet>0)then 
            c=min(50.,.35*zet) 
            psi_u=-((1+1.0*zet)**1.0+.667*(zet-14.28)/exp(c)+8.525)
          endif 
      end subroutine cal_psi_u 
      
      subroutine cal_psi_t(zet,psi_t)
  
          real(r8) :: x, psik,psic,f,zet,c,psi_t
          x=(1.-(15*zet))**.5 
          psik=2*log((1+x)/2) 
          x=(1.-(34.15*zet))**.3333 
          psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.) 
          f=zet*zet/(1+zet*zet) 
          psi_t=(1-f)*psik+f*psic   
         
          if(zet>0)then 
            c=min(50.,.35*zet) 
            psi_t=-((1.+2./3.*zet)**1.5+.6667*(zet-14.28)/exp(c)+8.525)
         endif
      end subroutine cal_psi_t
  
  
  
      
       subroutine cal_psi_q(zet,psi_q)
  
          real(r8) :: a,b,zet_0_q ,obklen,z0q,z,zet,psi_q
          a = 0.63
          b = 1
          z = 10
          z0q = 10**-6 
          zet_0_q =0               ! z0q/obklen
  
!          psi_q =a*(zet_0_q - zet) !+ (1-b)*log(z/z0q) 
          psi_q = -a*0.01*zet
  
      end subroutine cal_psi_q    
  
  
  
  
       subroutine cal_qsat(t,p,es,qsat )
     ! 计算饱和水汽压和饱和比湿的函数
      real(r8) :: t,p,es ,qsat
  
  
      es=6.112*exp(17.62*t/(t+243.12))*(1.006+3.15e-6*p-0.0074*1/p)
  
      qsat=1000*es*0.62198/(p-.378*es)
  
      end subroutine cal_qsat
      
      subroutine cal_qsee(t,p,es,qsee)
  
     !计算海面的饱和比湿
  
      real(r8) :: t,p,es,qsee
  
      
      
      es=6.112*exp(17.62*t/(t+243.12))*.98*(1.006+3.15e-6*p - 0.0074/p)
  
      qsee=1000*es*62198/(p-.378*es)     !1000 *() 是为了转化单位，从 kg/kg 到 g/kg
  
      end subroutine cal_qsee
     
  
        end module  physics_iap 

        module file_io
! 这个mdule 主要解决从文件中输入输出参数，分为单值和数组
        use shr_kind_mod , only : r8=>shr_kind_r8
                implicit none 

                private

                public danzhi
                public shuzu
                public filelinenum        
                public weiducanshu
                public array_1d
                public array_2d
                public array_b_dense
                public array_w_dense
        contains

        subroutine weiducanshu(filename,n,m,o)
!从文件中读取m,n,o
                character(len = 100),intent(in) ::filename
                integer,intent(inout)      ::  n,m,o 
                open(unit=13,file=filename)
                read(13,*) 
                read(13,*)
                read(13,*)
                read(13,*), m
                read(13,*)
                read(13,*), n
                read(13,*)
                read(13,*), o
 

        end subroutine weiducanshu

        subroutine danzhi(filename,var)
                !从文件filename中，读取一个数到var中
                character(len = 100),intent(in) ::filename
                real(r8)      ::  var 
                open(unit=10,file=filename)
                read(10,100) var
100             format(F8.4) 
!                print*,trim(filename)
!                print*,var
        end subroutine danzhi

        subroutine array_1d(filename,n,var)

                integer,intent(in)  :: n
                character(len = 100),intent(in)  :: filename 
                integer             :: i
                real(r8)                :: var(n)
                open(unit =15,file = filename)
                read(15,*) (var(i),i = 1,n)
!                print*,trim(filename)
!                print*,var

        end subroutine array_1d

        subroutine array_2d(filename,n,m,var)

                integer,intent(in)  :: n,m
                character(len = 100),intent(in)  :: filename 
                integer             :: i,j
                real(r8)                :: var(n,m)
                open(unit =16,file = filename)
                read(16,*) ((var(i,j),j=1,m),i = 1,n)
!                print*,trim(filename)
!                print*,var

        end subroutine array_2d

        subroutine array_b_dense(filename,n,o,var)

                integer,intent(in)  :: n,o
                character(len = 100),intent(in)  :: filename 
                integer             :: i,j
                real(r8)                :: var(n,o)
                open(unit =17,file = filename)

                do i = 1,o
                        do j = 1, n 
                                read(17,*) var(j,i)
                        enddo
                enddo

!                print*,trim(filename)
!                print*,var

        end subroutine array_b_dense

        subroutine array_w_dense(filename,n,m,o,var)

                integer,intent(in)  :: n,m,o
                character(len = 100),intent(in)  :: filename 
                integer             :: i,j,k
                real(r8)                :: var(n,m,o)
                open(unit =18,file = filename)

                do k = 1,o
                        do i = 1,n
                                read(18,*) (var(i,j,k),j=1,m)
                        enddo              
                enddo

!                print*,trim(filename)
!                print*,var

        end subroutine array_w_dense




        subroutine shuzu(filename,hangshu,lieshu,var)
        !从文件中读取多维数组,读取到数组var中
        
                integer,intent(in)  :: hangshu,lieshu
                character(len = 100),intent(in)  :: filename
                integer             :: i,j 
                real(r8)                :: var(hangshu,lieshu)
                open(unit =11,file = filename)
                !read(11,100) (var(i,j),i = 1,hangshu)
                do i =1,hangshu
                        do j =1,lieshu
                               read(11,100)  var(i,j)
                        enddo
                enddo
100             format(F8.4)
                print*,'在数组子程序中输出var = ',var


        end subroutine shuzu

        integer function filelinenum(a)
                integer ios
                character a*100
                open(22,file=trim(a))
                filelinenum=0
                do
                        read(22,*,iostat=ios)
                        if(ios/=0)exit
                        filelinenum=filelinenum+1
                end do
                close(22)

        end function filelinenum


      end module file_io

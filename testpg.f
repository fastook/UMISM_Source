      ipage=38
      iline=30
      do j=1,1000000
        do i=1,iline
          print *,i
        enddo
        do i=1,ipage-iline
          print *
        enddo
        call wait(100000)
      enddo
      end
      subroutine wait(n)
      do i=1,n
        x=sin(real(i))
      enddo
      end


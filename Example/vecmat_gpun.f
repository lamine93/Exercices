      program vecmat
      implicit real*8 (a-h,o-z)
c     parameter (ndim=10000)
c     parameter (ndim=100)
      allocatable a(:,:),v(:),vnew(:)
      read(5,*)ndim
      allocate (a(ndim,ndim),v(ndim),vnew(ndim))
      a(1,1)=1.0d0
      call init(a,ndim)
      v(1)=0.5d0
      v(2)=0.5d0
      v(3)=0.5d0
      v(4)=0.5d0
c     call prev(a,v,ndim,ndim,ndim)
      tol=1.0d-3
      error=1.0d0
      iter=0
      t1=omp_get_wtime()
!$omp target data map(to:a,v,vnew) 
      do while (error.gt.tol)
c .and.iter.lt.150)
!$omp target teams distribute parallel do if (ndim>100)
         do i=1,ndim
            vnew(i)=0.0d0
            do j=1,ndim
               vnew(i)=vnew(i)+a(j,i)*v(j)
            enddo
         enddo
!$omp end target teams distribute parallel do 
         dnorm=0.0d0
!$omp target teams distribute parallel do reduction(+:dnorm) map(dnorm)
         do i=1,ndim
            dnorm=dnorm+vnew(i)*vnew(i)
         enddo
!$omp end target teams distribute parallel do
         val=dsqrt(dnorm)
!$omp target teams distribute parallel do map(to:val) 
         do i=1,ndim
            vnew(i)=vnew(i)/val
         enddo
!$omp end target teams distribute parallel do
         error=0.0d0
!$omp target teams distribute parallel do reduction(+:error) map(error)
         do i=1,ndim
            error=error+(v(i)-vnew(i))**2
            v(i)=vnew(i)
         enddo
!$omp end target teams distribute parallel do
         iter=iter+1
c        write(6,'(I5,2F10.5)')iter,val,error
         if (mod(iter,100).eq.0)write(6,'(I5,2F10.5)')iter,val,error
      enddo
!$omp end target data
      t2=omp_get_wtime()
      write(6,'(A,2I10,F10.5,F10.2)')'Convergence after : ',ndim,iter,
     &val,t2-t1
      end

      subroutine init(a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n,n)
      do i=1,n-1
         j=i+1
         a(i,j)= 1.0d0
         a(j,i)= 1.0d0
      enddo
      return
      end

      subroutine prev(v,e,m,n,ndim)
c
c     ----- print out e and v-matrices
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(ndim,*),e(*)
      iwr=6
      max = 5
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,9008)
      write (iwr,9068) (e(i),i = imin,imax)
      write (iwr,9008)
      write (iwr,9028) (i,i = imin,imax)
      write (iwr,9008)
      do j = 1,n
      write (iwr,9048) j,(v(j,i),i = imin,imax)
      enddo
 140  if (imax .lt. m) go to 100
      return
 9008 format(/)
 8028 format(17x,10(3x,i3,3x))
 8048 format(i5,12x,10f9.4)
 8068 format(17x,10f9.4)
 9028 format(17x,5(9x,i3,9x))
 9048 format(i5,2x,10x,5f21.10)
 9068 format(17x,5f21.10)
      end


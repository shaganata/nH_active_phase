!-----------------------------------------------------------------------
!           Fitting EGRESS/INGRESS nH during active phases
!      by polynomial expression: n_H(b) = n1*b^{-1} + nK*b^{-K},
!                                b is impact parameter
!             input file columns: fi, nH, nHhigh, nHlow
!             output file columns: n_1, n_K, K, chi_red, max. rel. error
!                                  and xi for 50 best fits
! Natalia Shagatova 2020
!
! GNU General Public License v3.0
! https://opensource.org/licenses/GPL-3.0
!-----------------------------------------------------------------------

program nHinversionFit
!-----------------------------------------------------------------------

implicit none               
                              
! mH = mass of H atom, mi = mean molecular weight                                              
real(8), parameter :: pi=3.14159265358979d0, mH=1.67372d-24, mi=1.4d0         
! p = separation of binary components, Rgiant = giant radius
! Rp = Rgiant/p, pR = p/Rgiant
real(8), parameter :: p=2.784d13, Rgiant=6.96d12, Rp=0.25d0, pR=4.d0 

! coefficients in n_H(b) as elements of array: n(1)=n1, n(2)=nK, n(3)=K
! nH = column density, CHIred = CHI^2_red
! er = rel.error, ermax = max. rel. error, maxCHI = max. CHI^2_red
! A =  constant in nH integral, xi = parameter in velocity profile v(r)
real(8) :: n(3), nH, CHIred, er, ermax, maxCHI, A, xi     
! arrays to load observed values of nH        
real(8), allocatable :: fi(:), nHobserved(:), nHhigh(:), nHlow(:)
! arrays to model column densities, b is impact parameter
real(8), allocatable :: b(:), nHtot(:), deltanH(:)
! arrays for saving 50 best values of n1, nK, K, CHIred, xi, and K
! nN(1,i)=n1,nN(2,i)=nK,nN(3,i)=CHIred,nN(4,i)=ermax,nN(5,i)=xi,kN(i)=K
real(8) :: nN(5,50)            
integer ::  kN(50), K

! variables to generate random numbers        
integer ::  seed(12), clock, j, m

character(40) :: input, output, head

! numlines = number of lines; integers for loops: j2, j3, j4
! maxi = number of models generated is 5*maxi
integer :: IOstatus, numlines, j2, j3, j4, maxi

! i = orbital inclination
! nHwd = contribution to the nH from the hot component
! lambda(50) = array for eigenvalues of Abel operator
! and the geometry variables (Shagatova et al. 2016, A&A 588, A83)
real(8) :: i, nHwd, lambda(50), theta, l, l2    
! variables for numerical integration
real(8) :: suma, suma0, k1, k2, k3, step, y      

!--------------------input and output-----------------------------------

! enter input file, i and nHwd:
input='all-nh_fi_EGRESS.dat'  
i=90.0d0 
nHwd=1.5d22      

!-------------------------------------------

write(*,*) 'Input file is:', input  
write(*,*) 'Enter the name of output file:'
read(*,*) output

i=i*pi/180.d0      

open(1,file=input)
open(2,file=output, status='new')

!---------------seed to generate random numbers-------------------------  

m=12                             
   
CALL RANDOM_SEED(size = m)
 
CALL SYSTEM_CLOCK(COUNT=clock)
 
seed = clock + 37 * (/ (j - 1, j = 1, m) /)

CALL RANDOM_SEED(PUT = seed)

!-----------------loading values to arrays------------------------------

read(1,*) head                                  !reading the header line
 
allocate(fi(1))
allocate(nHobserved(1))
allocate(nHhigh(1))   
allocate(nHlow(1))
numlines=0         

! counting the number of lines = length of arrays
do
    read(1,*,iostat=IOstatus) fi(1), nHobserved(1), nHhigh(1), nHlow(1)	    
    if (IOstatus==0) then   
	  numlines=numlines+1
	end if 
	
    if (IOstatus>0) then
      write(*,*) '...error in data...'
      exit
    end if
   
    if (IOstatus<0) then
	  exit
    end if
end do

deallocate(fi)
deallocate(nHobserved)
deallocate(nHhigh)
deallocate(nHlow)

rewind(1)                                               
allocate(fi(numlines))
allocate(nHobserved(numlines))
allocate(nHhigh(numlines))
allocate(nHlow(numlines))
allocate(deltanH(numlines))
allocate(b(numlines))
allocate(nHtot(numlines))

read(1,*) head                                 ! reading the header line
do j=1,numlines,1
  read(1,*) fi(j), nHobserved(j), nHhigh(j), nHlow(j)
  nHobserved(j)=nHobserved(j)-nHwd
  fi(j)=fi(j)*2.d0*pi
  b(j)=pR*dsqrt( (dcos(i)**2.d0)+((dsin(fi(j))*dsin(i))**2.d0) )                                         
  deltanH(j)=(nHhigh(j)-nHlow(j))
end do 
        
deallocate(nHhigh)
deallocate(nHlow)

!reading the eigenvalues of Abel operator
open(3,file='eigenvaluesAbel.txt')    
do j=1,50,1
   read(3,*) lambda(j)
end do

!------assigning the values of parameters by Monte Carlo method---------

write(2,*)'# nH = n1*b^{-1} + nK*b^{-K}'    
write(2,*)'#    n1          nK       K    CHIred      ermax         xi'
nN=0.d0
kN=0
do j=1,50,1                ! array is filled by large values of chi2_red
    nN(3,j)=1.d100
end do


do j4=1,5,1    ! loop over j4 aid to see the fraction of iterations done
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maxi=10000                        ! total number of iterations is 5*maxi

! random numbers assignment to fit. parameters within given limits
do j2=1,maxi,1     

   call FindRandomNumbers
   
   A=n(1)/lambda(1)
   xi=n(2)*lambda(1)/(lambda(K)*n(1))
   
   !------------------evaluation of nHtot(i,phi)------------------------
   
   do j3=1,numlines,1
   
      if ((fi(j3)>=(3.d0*pi/2.d0)) .or. (fi(j3)<=(pi/2.d0))) then
         theta=dasin(Rgiant*b(j3)/p)                
      else
         theta=((2*dasin(dsqrt((dcos(i)**2.d0)+&
         ((dsin(pi/2.d0)*dsin(i))**2.d0))))-dasin(Rgiant*b(j3)/p))
      end if   
    
      if (fi(j3)>pi) then
         theta=-theta        
      end if  
            
     !...............................................................
     step=1.d12     
     suma=0.d0
     suma0=1.d0          ! to fulfill the while cycle condition at start  
     
     ! white dwarf position = starting point of integration
     if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then            
	     l=dsqrt((p**2.d0)-((Rgiant*b(j3))**2.d0))  
     else      
	     l=-dsqrt((p**2.d0)-((Rgiant*b(j3))**2.d0)) 
     end if    
     
     ! numerical integration - Runge Kutta 4th order                                                  
     do while (dabs((suma-suma0)/suma0)>=0.0000001d0)                        
        suma0=suma
        l2=l
        k1=nHtotF(l,j3)
        l=l-(0.5d0*step)
	    k2=nHtotF(l,j3)
	    l=l-(0.5d0*step)
	    k3=nHtotF(l,j3)
        suma=suma+((step/6.d0)*(k1+(4.d0*k2)+k3))	  
     end do
   
     y=Abs(suma-suma0)
 
     !......................................................
     
     suma=suma*A*0.5d0*Rgiant                           
     nHtot(j3)=suma
   end do
     
   
   !------------------------------------------
   call EvaluateCHIred
   
   ermax=0.d0                 !computation of the maximal relative error
   do  j=1,numlines,1
	     er=dabs((nHtot(j))-nHobserved(j))/nHobserved(j)
	     if(er>ermax) then
	         ermax=er
	     end if    
   end do

   
   maxCHI=0.d0      
   do j=1,50,1                     ! finding the fit with highest CHIred
      if ((nN(3,j)>maxCHI) ) then  
         maxCHI=nN(3,j)
	     m=j
	  end if
   end do						
			
! substitution of the worst fit parameters by better ones  							
   if ((CHIred<maxCHI) ) then       
	   nN(1,m)=n(1)                 ! 50 best values are saved in arrays
	   nN(2,m)=n(2)
	   nN(3,m)=CHIred
	   nN(4,m)=ermax              
	   nN(5,m)=xi
	   kN(m)=K
   end if
   
 
   
end do

write(*,'(A28,I2,A4)') "Fraction of iterations done:", j4, " / 5"
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end do

! ordering best fits by CHIred value
do j=1,50,1            
   maxCHI=0.d0
   do j2=1,50,1     
      if (nN(3,j2)>maxCHI) then
         maxCHI=nN(3,j2)
	     m=j2
	  end if
   end do
   write (2,'(2ES12.4,I4,3ES12.4)') nN(1,m), nN(2,m), kN(m), nN(3,m),&
   nN(4,m), nN(5,m)    
   nN(3,m)=0.d0
end do

close(1)
close(2)
close(3) 

!-----------------------------------------------------------------------
!                       subroutines
!-----------------------------------------------------------------------
contains
     
	subroutine FindRandomNumbers
         real(8) :: n1l, n1h, n1delta, nKl, nKh, nKdelta, Kl, Kh, Kdelta
         ! the ranges of fitting parameters values
         n1l=dlog10(1.0d23)
         n1h=dlog10(1.0d24)        
         nKl=dlog10(1.0d23)
         nKh=dlog10(1.0d27)         
         Kl=1.d0
         Kh=20.d0
         
         n1delta=n1h-n1l
         nKdelta=nKh-nKl
         Kdelta=Kh-Kl
         CALL RANDOM_NUMBER(n)		 
	     n(1)=10**((n(1)*n1delta)+n1l)
	   	 n(2)=10**((n(2)*nKdelta)+nKl)
	   	 n(3)=(n(3)*Kdelta)+Kl    
		 K=nint(n(3))                                  ! nearest integer
		 !K=7                  
  	end subroutine

	subroutine EvaluateCHIred
	    integer :: ii
		CHIred=0.d0
		do ii=1,numlines,1
		CHIred=CHIred+(((nHtot(ii)&
		-nHobserved(ii))/deltanH(ii))**2.d0)
		end do
		CHIred=CHIred/(numlines-4.d0)  ! DoF for 3 parameters: n1,nK,K
	end subroutine
	
! sub-integral function for nH computation, v(r) from Abel. inversion
    function nHtotF (ll, ii)       
        real(8) :: ll, nHtotF
		integer :: ii
        nHtotF=(1.d0/(((b(ii)*Rgiant)**2.d0)+(ll**2.d0)))&
		*(1.d0+(xi*(Rgiant**(dble(K)-1.d0))&    !
        /((((b(ii)*Rgiant)**2.d0)+(ll**2.d0))**((dble(K)-1.d0)/2.d0)))) 
    end function


!--------------------------------------------------------------
end program

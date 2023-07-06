!-----------------------------------------------------------------------
!              Model of nH(fi,i) during active phases
!- egress/ingress v(r) from inversion of Abel operator 
!  (Knill et al. 1993, A&A 274, 1002)
!- smooth transitional v(r) in 2 regions
!  (Shagatova et al. 2017, A&A 602, A71)          
!
!       input: orb. inclination, E/IN v(r), nHwd, psiE/IN
!       output file columns: phi, nH, nH_error, transition
!                                
! Natalia Shagatova 2020
!
! GNU General Public License v3.0
! https://opensource.org/licenses/GPL-3.0
!-----------------------------------------------------------------------

program columndensityall
!-----------------------------------------------------------------------

implicit none 

! mH = mass of H atom, mi = mean molecular weight                                              
real(8), parameter :: pi=3.14159265358979d0, mH=1.67372d-24, mi=1.4d0         
! p = separation of binary components, Rgiant = giant radius
! Rp = Rgiant/p
real(8), parameter :: p=2.784d13, Rgiant=6.96d12, Rp=0.25d0

! n - integer for loop in orb. phase, dof=degrees of freedom
integer :: n=0, dof, IOstatus
                       
! fi = orbital phase, b = impact parameter, v = wind velocity
! Mdot = mass loss rate, nHobs = measured column density
! nHmin, nHplus = -/+ errors of measured column density
! nHwd = contribution to the nH from white dwarf wind
! CHIred= Chi^2 reduced, dof = degrees of freedom
real(8) :: fi, b, v, Mdot, nHobs, nHmin, nHplus, nHwd, CHIred  
! transition = variable to write the v(r) used at the end of integration
!              at given fi
! header = variable to read the header line in file
character(7) :: transition, header        

! parameters of unified velocity profile, vterE/IN=terminal egr./ingr. v
! (Shagatova et al. 2017, A&A 602, A71) 
! psi is the angle in the coordin. system in the plane of observations
! centered in the giant position
real(8) :: vterE, vterIN, xiE, Ke, xiIN, Kin, n1e, n1in, C
real(8) :: psi, psiE, psiIN, psiE1, psiE2, psiE3, psiIN1, psiIN2, psiIN3          

! ii = orbital inclination [°], i=orb. incl. [rad], i90=(90.d0-ii)[rad]
! and the geometry variables (Shagatova et al. 2016, A&A 588, A83)
real(8) :: i, i90, ii, theta, l, l2, x, xp, y, x3, y3   
! variables for numerical integration
real(8) :: suma, suma0, k1, k2, k3, step, nHerr  

!-----INPUT: orb. inclination, E/IN v(r), nHwd ; OUTPUT file name-------

ii=90.0d0 

! EGRESS column density model and v(r)    
xiE=3.56d1
Ke=3.d0
n1e=3.84d23

! INGRESS column density model and v(r) 
xiIN=1.86d2
Kin=4.d0
n1in=1.79d23

! the value of nHwd have to be the same as used for egress/ingress nH 
! fitting in nH_inversion_fit_5_2.f90 code
nHwd=1.5d22
         
open(1,file="nH_act_phase_i90_wd1p5e22.txt",status="replace")

!................................   

! values of angle psi defining the extent of regions with "transitional"
! velocity profile
psiE=1.1d0      !1.1d0  
psiIN=0.8d0     !0.8d0 
 
vterE=2.0d6     
vterIN=(n1e/n1in)*vterE   ! condition ensuring constant Mdot

i90=(90.d0-ii)*pi/180.d0                               
i=ii*pi/180.d0                                   

psiE1=psiE
psiE2=pi-psiE
psiE3=(2.d0*pi)+psiE 
psiIN2=pi+psiIN
psiIN1=-psiIN
psiIN3=(2.d0*pi)-psiIN

Mdot=2.d0*pi*mi*mH*(1.d0/1.57d0)*Rgiant*n1e*vterE

C=pi/( 2.d0*( psiIN1-psiE1 ) )


!----------------------------------------------------------------------- 
!                       evaluation of Chi^2_red
!-----------------------------------------------------------------------

CHIred=0.d0
dof=0
open(2,file="nH_act_phase_measured.txt",status="old")
read(2,*) header

do
  read(2,*,iostat=IOstatus) fi, nHobs, nHplus, nHmin	  
     
  if (IOstatus==0) then   
   !-----------------------evaluation of nH(fi)------------------------
   fi=fi*2.d0*pi 
   dof=dof+1
   
   b=p*dsqrt( (dcos(i)**2.d0)+((dsin(fi)*dsin(i))**2.d0) )
   
   if ((fi>=(3.d0*pi/2.d0)) .or. (fi<=(pi/2.d0))) then
      theta=dasin(b/p)                
   else
     theta=((2*dasin(dsqrt((dcos(i)**2.d0)+&
     ((dsin(pi/2.d0)*dsin(i))**2.d0))))-dasin(b/p))
   end if   
    
   if (fi>pi) then
      theta=-theta           
   end if
     
   suma=0.d0
   suma0=0.d0
      
   step=1.d12
   
   ! white dwarf position = starting point of integration
   if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then             
	  l=dsqrt((p**2.d0)-(b**2.d0))  
   else      
	  l=-dsqrt((p**2.d0)-(b**2.d0)) 
   end if    
   
   ! during inferior conjunction of the giant, starting point of 
   ! integration = giant surface
   if ((b<Rgiant) .and. (theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then
      l=-dsqrt((Rgiant**2.d0)-(b**2.d0))
   end if
     
   suma=0.d0
   suma0=1.d0            ! to fulfill the while cycle condition at start                       
  
   do while (dabs((suma-suma0)/suma0)>=0.0000001d0)   

      suma0=suma
      l2=l

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~k1
      
      ! radial distance from the white dwarf
      if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then           
    	x=dsqrt(1.d0-((b/p)**2.d0))-(l/p)   
      else      
    	x=-dsqrt(1.d0-((b/p)**2.d0))-(l/p)  
      end if    
  	
      xp=x*dcos(theta)
      y=x*dsin(theta)
  	
  	  ! coordinate system in the plane of observations
  	  ! centered in the giant position
      x3= ((x*dcos(fi)-1.d0)*dcos(i90))&
  	     + (x*dsin(fi)*dsin(i90))
      y3=-((x*dcos(fi)-1.d0)*dsin(i90))&
  	     + (x*dsin(fi)*dcos(i90))
  	    
      if ( (x3>=0.d0) .and. (y3>=0.d0) ) then ! Quadrant I
  	    psi=datan(y3/x3)
      else if (x3<0.d0) then   ! Quadrant II and III
  	    psi=datan(y3/x3)+pi
      else
  	    psi=datan(y3/x3)+(2.d0*pi)  ! Quadrant IV
      end if
	
	
      if (psi<=psiIN2)  then
         if (psi<psiE1) then
            k1=1.d0/(((l**2.d0)+(b**2.d0))*&
               ( (Vegress(x,theta)*(dcos(C*(psi-psiE1))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(-psiIN1+psi))**2.d0)) )  )
		    transition="    v_1"     !trans. velocity profile in Quadrant I
         else if (psi<=psiE2) then      
		    k1=NHe(l)/vterE
	        transition=" egress"                                  
         else	 
		    k1=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE2))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(psiIN2-psi))**2.d0)) ) )
		    transition="   v_23"  !tr. velocity profile in Quadrants II-III
         end if
      end if	

      if (psi>psiIN2) then	
         if (psi<psiIN3) then
            k1=NHin(l)/vterIN                                 
	        transition="ingress"
         else
		    k1=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE3))**2.d0))+& 
               (Vingress(x,theta)*(dcos(C*(psiIN3-psi))**2.d0)) )  )
		    transition="    v_4"    !trans. velocity profile in Quadrant IV
         end if 
      end if 
	
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~k2	
	
      l=l-(0.5d0*step)
	
      if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then            
	    x=dsqrt(1.d0-((b/p)**2.d0))-(l/p)   
      else      
	    x=-dsqrt(1.d0-((b/p)**2.d0))-(l/p)  
      end if    
	    
      xp=x*dcos(theta)
      y=x*dsin(theta)
	    
      x3= ((x*dcos(fi)-1.d0)*dcos(i90))&
	         + (x*dsin(fi)*dsin(i90))
      y3=-((x*dcos(fi)-1.d0)*dsin(i90))&
	         + (x*dsin(fi)*dcos(i90))
	     
	  if ( (x3>=0.d0) .and. (y3>=0.d0) ) then   
	      psi=datan(y3/x3)
	  else if (x3<0.d0) then   
	      psi=datan(y3/x3)+pi
	  else
	      psi=datan(y3/x3)+(2.d0*pi)  
	  end if
	
      if (psi<=psiIN2)  then
         if (psi<psiE1) then
            k2=1.d0/(((l**2.d0)+(b**2.d0))*&
               ( (Vegress(x,theta)*(dcos(C*(psi-psiE1))**2.d0))+&    
               (Vingress(x,theta)*(dcos(C*(-psiIN1+psi))**2.d0)) )  )
	        transition="    v_1"		       
         else if (psi<=psiE2) then      
		    k2=NHe(l)/vterE
		    transition=" egress"
         else	 
	        k2=1.d0/(((l**2.d0)+(b**2.d0))*&
	           ( (Vegress(x,theta)*(dcos(C*(psi-psiE2))**2.d0))+& 
               (Vingress(x,theta)*(dcos(C*(psiIN2-psi))**2.d0)) ) )
		    transition="   v_23"
         end if
      end if	

      if (psi>psiIN2) then	
         if (psi<psiIN3) then
            k2=NHin(l)/vterIN 
		    transition="ingress"
         else
		    k2=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE3))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(psiIN3-psi))**2.d0)) )  )
		    transition="    v_4"
         end if 
      end if 	 
	 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~k3	
	
      l=l-(0.5d0*step)
   	
      if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then            
    	x=dsqrt(1.d0-((b/p)**2.d0))-(l/p)   
      else      
    	x=-dsqrt(1.d0-((b/p)**2.d0))-(l/p)  
      end if    
   	
      xp=x*dcos(theta)
      y=x*dsin(theta)
   	
      x3= ((x*dcos(fi)-1.d0)*dcos(i90))&
   	     + (x*dsin(fi)*dsin(i90))
      y3=-((x*dcos(fi)-1.d0)*dsin(i90))&
   	     + (x*dsin(fi)*dcos(i90))
   
      if ( (x3>=0.d0) .and. (y3>=0.d0) ) then    
   	    psi=datan(y3/x3)
      else if (x3<0.d0) then   
   	    psi=datan(y3/x3)+pi
      else
   	    psi=datan(y3/x3)+(2.d0*pi)  
   	  end if
	
      if (psi<=psiIN2)  then
         if (psi<psiE1) then
            k3=1.d0/(((l**2.d0)+(b**2.d0))*&
               ( (Vegress(x,theta)*(dcos(C*(psi-psiE1))**2.d0))+&  
               (Vingress(x,theta)*(dcos(C*(-psiIN1+psi))**2.d0)) )  ) 
	        transition="    v_1" 
         else if (psi<=psiE2) then      
		    k3=NHe(l)/vterE                                  
	        transition=" egress"
         else	 
	        k3=1.d0/(((l**2.d0)+(b**2.d0))*&
	           ( (Vegress(x,theta)*(dcos(C*(psi-psiE2))**2.d0))+&      
               (Vingress(x,theta)*(dcos(C*(psiIN2-psi))**2.d0)) ) )
	        transition="   v_23"
         end if
      end if	

      if (psi>psiIN2) then	
         if (psi<psiIN3) then
            k3=NHin(l)/vterIN 
		    transition="ingress"
         else
		    k3=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE3))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(psiIN3-psi))**2.d0))  ) )  
	        transition="    v_4"
         end if 
      end if 		 
		 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	 				 
      suma=suma+((step/6.d0)*(k1+(4.d0*k2)+k3))
      
   end do  
   suma=suma*Mdot/(4.d0*pi*mi*mH)   
   suma=suma+nHwd
          
   !-------------------------------------------------------------------  
   CHIred=CHIred+(((suma-nHobs)/(nHmin+nHplus))**2.d0)  
   !-------------------------------------------------------------------	  
  end if 
	
  if (IOstatus>0) then
    write(*,*) '...error in data reading...'
    exit
  end if
   
  if (IOstatus<0) then
    exit
  end if
end do

dof=dof-7-1   !parameters: n1e, n1in, nKe, nKin, Ke, Kin, nHwd
CHIred=CHIred/dof	

write(*,"(A11,F8.3,A15,I5)") "Chi^2_red =", CHIred, "DoF =", dof 

!-----------------------  OUTPUT FILE: header --------------------------

write(1,"(A4,F4.1,A7,ES8.2E2,A7,ES8.2E2,A7,F2.0,A7,ES8.2E2,A7,ES8.2E2,&
      A7,F2.0,A7,ES8.2E2,A8,F4.2,A9,F4.2)") "# i=", ii, "n1e=", n1e,&
       "xiE=", xiE, "Ke=", Ke, "n1in=", n1in, "xiIN=", xiIN, "Kin=", &
       Kin, "nHwd=", nHwd, "psiE=", psiE, "psiIN=", psiIN   
write(1,"(A13,F8.3,A15,I5)") "# Chi^2_red =", CHIred, "DoF =", dof     
write(1,"(A53)") "#   phi           nH          nH_error     transition" 
 
 
!----------------------------------------------------------------------- 
!                  evaluation of nH(fi) for fi <0,1>
!-----------------------------------------------------------------------

do n=0,360,1
   fi=real(n)*pi/180.0d0
   b=p*dsqrt( (dcos(i)**2.d0)+((dsin(fi)*dsin(i))**2.d0) )
   
   if ((fi>=(3.d0*pi/2.d0)) .or. (fi<=(pi/2.d0))) then
      theta=dasin(b/p)                
   else
     theta=((2*dasin(dsqrt((dcos(i)**2.d0)+&
     ((dsin(pi/2.d0)*dsin(i))**2.d0))))-dasin(b/p))
   end if   
    
   if (fi>pi) then
      theta=-theta           
   end if
     
   suma=0.d0
   suma0=0.d0
      
   step=1.d12
   
   ! white dwarf position = starting point of integration
   if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then             
	  l=dsqrt((p**2.d0)-(b**2.d0))  
   else      
	  l=-dsqrt((p**2.d0)-(b**2.d0)) 
   end if    
   
   ! during inferior conjunction of the giant, starting point of 
   ! integration = giant surface
   if ((b<Rgiant) .and. (theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then
      l=-dsqrt((Rgiant**2.d0)-(b**2.d0))
   end if
     
   suma=0.d0
   suma0=1.d0            ! to fulfill the while cycle condition at start                       
  
   do while (dabs((suma-suma0)/suma0)>=0.0000001d0)   

      suma0=suma
      l2=l

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~k1
      
      ! radial distance from the white dwarf
      if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then           
    	x=dsqrt(1.d0-((b/p)**2.d0))-(l/p)   
      else      
    	x=-dsqrt(1.d0-((b/p)**2.d0))-(l/p)  
      end if    
  	
      xp=x*dcos(theta)
      y=x*dsin(theta)
  	
  	  ! coordinate system in the plane of observations
  	  ! centered in the giant position
      x3= ((x*dcos(fi)-1.d0)*dcos(i90))&
  	     + (x*dsin(fi)*dsin(i90))
      y3=-((x*dcos(fi)-1.d0)*dsin(i90))&
  	     + (x*dsin(fi)*dcos(i90))
  	    
      if ( (x3>=0.d0) .and. (y3>=0.d0) ) then ! Quadrant I
  	    psi=datan(y3/x3)
      else if (x3<0.d0) then   ! Quadrant II and III
  	    psi=datan(y3/x3)+pi
      else
  	    psi=datan(y3/x3)+(2.d0*pi)  ! Quadrant IV
      end if
	
	
      if (psi<=psiIN2)  then
         if (psi<psiE1) then
            k1=1.d0/(((l**2.d0)+(b**2.d0))*&
               ( (Vegress(x,theta)*(dcos(C*(psi-psiE1))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(-psiIN1+psi))**2.d0)) )  )
		    transition="    v_1"     !trans. velocity profile in Quadrant I
         else if (psi<=psiE2) then      
		    k1=NHe(l)/vterE
	        transition=" egress"                                  
         else	 
		    k1=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE2))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(psiIN2-psi))**2.d0)) ) )
		    transition="   v_23"  !tr. velocity profile in Quadrants II-III
         end if
      end if	

      if (psi>psiIN2) then	
         if (psi<psiIN3) then
            k1=NHin(l)/vterIN                                 
	        transition="ingress"
         else
		    k1=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE3))**2.d0))+& 
               (Vingress(x,theta)*(dcos(C*(psiIN3-psi))**2.d0)) )  )
		    transition="    v_4"    !trans. velocity profile in Quadrant IV
         end if 
      end if 
	
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~k2	
	
      l=l-(0.5d0*step)
	
      if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then            
	    x=dsqrt(1.d0-((b/p)**2.d0))-(l/p)   
      else      
	    x=-dsqrt(1.d0-((b/p)**2.d0))-(l/p)  
      end if    
	    
      xp=x*dcos(theta)
      y=x*dsin(theta)
	    
      x3= ((x*dcos(fi)-1.d0)*dcos(i90))&
	         + (x*dsin(fi)*dsin(i90))
      y3=-((x*dcos(fi)-1.d0)*dsin(i90))&
	         + (x*dsin(fi)*dcos(i90))
	     
	  if ( (x3>=0.d0) .and. (y3>=0.d0) ) then   
	      psi=datan(y3/x3)
	  else if (x3<0.d0) then   
	      psi=datan(y3/x3)+pi
	  else
	      psi=datan(y3/x3)+(2.d0*pi)  
	  end if
	
      if (psi<=psiIN2)  then
         if (psi<psiE1) then
            k2=1.d0/(((l**2.d0)+(b**2.d0))*&
               ( (Vegress(x,theta)*(dcos(C*(psi-psiE1))**2.d0))+&    
               (Vingress(x,theta)*(dcos(C*(-psiIN1+psi))**2.d0)) )  )
	        transition="    v_1"		       
         else if (psi<=psiE2) then      
		    k2=NHe(l)/vterE
		    transition=" egress"
         else	 
	        k2=1.d0/(((l**2.d0)+(b**2.d0))*&
	           ( (Vegress(x,theta)*(dcos(C*(psi-psiE2))**2.d0))+& 
               (Vingress(x,theta)*(dcos(C*(psiIN2-psi))**2.d0)) ) )
		    transition="   v_23"
         end if
      end if	

      if (psi>psiIN2) then	
         if (psi<psiIN3) then
            k2=NHin(l)/vterIN 
		    transition="ingress"
         else
		    k2=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE3))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(psiIN3-psi))**2.d0)) )  )
		    transition="    v_4"
         end if 
      end if 	 
	 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~k3	
	
      l=l-(0.5d0*step)
   	
      if ((theta>=-(pi/2.d0)) .and. (theta<=(pi/2.d0))) then            
    	x=dsqrt(1.d0-((b/p)**2.d0))-(l/p)   
      else      
    	x=-dsqrt(1.d0-((b/p)**2.d0))-(l/p)  
      end if    
   	
      xp=x*dcos(theta)
      y=x*dsin(theta)
   	
      x3= ((x*dcos(fi)-1.d0)*dcos(i90))&
   	     + (x*dsin(fi)*dsin(i90))
      y3=-((x*dcos(fi)-1.d0)*dsin(i90))&
   	     + (x*dsin(fi)*dcos(i90))
   
      if ( (x3>=0.d0) .and. (y3>=0.d0) ) then    
   	    psi=datan(y3/x3)
      else if (x3<0.d0) then   
   	    psi=datan(y3/x3)+pi
      else
   	    psi=datan(y3/x3)+(2.d0*pi)  
   	  end if
	
      if (psi<=psiIN2)  then
         if (psi<psiE1) then
            k3=1.d0/(((l**2.d0)+(b**2.d0))*&
               ( (Vegress(x,theta)*(dcos(C*(psi-psiE1))**2.d0))+&  
               (Vingress(x,theta)*(dcos(C*(-psiIN1+psi))**2.d0)) )  ) 
	        transition="    v_1" 
         else if (psi<=psiE2) then      
		    k3=NHe(l)/vterE                                  
	        transition=" egress"
         else	 
	        k3=1.d0/(((l**2.d0)+(b**2.d0))*&
	           ( (Vegress(x,theta)*(dcos(C*(psi-psiE2))**2.d0))+&      
               (Vingress(x,theta)*(dcos(C*(psiIN2-psi))**2.d0)) ) )
	        transition="   v_23"
         end if
      end if	

      if (psi>psiIN2) then	
         if (psi<psiIN3) then
            k3=NHin(l)/vterIN 
		    transition="ingress"
         else
		    k3=1.d0/(((l**2.d0)+(b**2.d0))*&
		       ( (Vegress(x,theta)*(dcos(C*(psi-psiE3))**2.d0))+&
               (Vingress(x,theta)*(dcos(C*(psiIN3-psi))**2.d0))  ) )  
	        transition="    v_4"
         end if 
      end if 		 
		 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	 				 
      suma=suma+((step/6.d0)*(k1+(4.d0*k2)+k3))
      
   end do  
   suma=suma*Mdot/(4.d0*pi*mi*mH)   

   suma0=suma0*Mdot/(4.d0*pi*mi*mH) 
   nHerr=Abs(suma-suma0)
   !...................................................... 
   
   fi=(fi/(2.d0*pi))                                    
   write(1,"(F9.5)", advance="no") fi         
   write(1,"(ES14.3E2)", advance="no") suma
   write(1,"(ES16.3E2)", advance="no") nHerr
   write(1,"(A10)") transition
   
end do


!-----------------------------------------------------------------------

close(1)
close(2)

!-----------------------------------------------------------------------
!                           subroutines
!-----------------------------------------------------------------------
contains

! sub-integral function for nH computation with EGRESS velocity profile
function NHe (ll)                                           
real(8) :: ll, NHe
NHe=(1.d0/((b**2.d0)+(ll**2.d0)))*(1.d0+(xiE*(Rgiant**(Ke-1.d0))&    
/((((b**2.d0)+(ll**2.d0))**((Ke-1.d0)/2.d0))))) 
end function

! sub-integral function for nH computation with INGRESS velocity profile
function NHin (ll)                                           
real(8) :: ll, NHin
NHin=(1.d0/((b**2.d0)+(ll**2.d0)))*(1.d0+(xiIN*(Rgiant**(Kin-1.d0))&    
/((((b**2.d0)+(ll**2.d0))**((Kin-1.d0)/2.d0))))) 
end function

! egress velocity profile
function Vegress (xx, tt)
real(8) :: xx, tt, Vegress
Vegress=vterE/(1.d0 + (xiE*((Rp/dsqrt((xx**2.d0)+1.d0-&
        (2.d0*xx*dcos(tt))))**(Ke-1.d0)) ) )
end function

! ingress velocity profile
function Vingress (xx, tt)
real(8) :: xx, tt, Vingress
Vingress=vterIN/(1.d0 + (xiIN*((Rp/dsqrt((xx**2.d0)+1.d0-&
         (2.d0*xx*dcos(tt))))**(Kin-1.d0)) ) )
end function

!-----------------------------------------------------------------------
end program

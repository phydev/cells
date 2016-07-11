module misc_m

  use global_m

  implicit none

  private

  public ::  hfield_calc, hfunc, deltak,&
       heaviside, ran2, output_aux, spherical_surface, gen_cell_points

  contains

    subroutine hfield_calc(cell, aux, r, Lsize, lxyz, lxyz_inv, lxyz_part, lxyz_inv_part, ncell, tcell, ntype, np )

      type(mesh_t), allocatable, intent(in) :: cell(:,:)
      integer, allocatable, intent(in) :: r(:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:), lxyz_part(:,:), lxyz_inv_part(:,:), ncell(:)
      type(mesh_t), allocatable, intent(inout) :: aux(:,:)
      integer, intent(in) ::  ntype, np, tcell, Lsize(2)
      integer :: icell, itype, b(2), ip2, ip, rmin(2)

      aux(:,:)%phi = 0.d0
      do ip=1, np
         
         do itype = 1, ntype
          
            do icell=1 + (itype-1)*(ncell(itype-1)), ncell(itype) + (itype-1)*(ncell(itype-1))

               b(1:2) = abs(lxyz(ip,1:2) - lxyz(r(icell),1:2)) ! r position of the icell
               rmin(1:2) = min(b(1:2),2*Lsize(1:2)-b(1:2)-1)
               ip2 = lxyz_inv_part(rmin(1),rmin(2))
               aux(ip,itype)%phi = aux(ip,itype)%phi + hfunc(cell(ip2,icell)%phi)
               
            end do
         end do
      end do

    end subroutine hfield_calc


    function hfunc(x)
      
      real :: hfunc, x,y
	y = 0.5*(x+1)
	hfunc = y**2*(3-2*y)
      !hfunc = (0.25)*(1.d0+x)*(1.d0-x)
    end function hfunc

    function deltak(l,m)

      implicit none

      real :: deltak
      integer :: l, m

      if(l.eq.m) then
         deltak = 1.d0
      else
         deltak = 0.d0
      end if

    end function deltak

    function heaviside(x)
      
      implicit none
      
      real :: x, heaviside
      
      if ( x<0.d0) then
         heaviside = 0.d0
      else
         heaviside = 1.d0
      end if
    end function heaviside

    function ran2(idum) 
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      real :: ran2,AM,EPS,RNMX 
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, & 
           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, & 
           NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS) 
      integer :: idum2,j,k,iv(NTAB),iy 
      save iv,iy,idum2 
      data idum2/123456789/, iv/NTAB*0/, iy/0/ 
      if (idum.le.0) then 
         idum=max(-idum,1) 
         idum2=idum 
         do 11 j=NTAB+8,1,-1 
            
            k=idum/IQ1 
            idum=IA1*(idum-k*IQ1)-k*IR1 
            if (idum.lt.0) idum=idum+IM1 
            if (j.le.NTAB) iv(j)=idum 
11          continue 
            iy=iv(1) 
         endif 
         k=idum/IQ1 
         idum=IA1*(idum-k*IQ1)-k*IR1 
         if (idum.lt.0) idum=idum+IM1 
         k=idum2/IQ2 
         idum2=IA2*(idum2-k*IQ2)-k*IR2 
         if (idum2.lt.0) idum2=idum2+IM2 
         j=1+iy/NDIV 
         iy=iv(j)-idum2 
         iv(j)=idum 
         if(iy.lt.1)iy=iy+IMM1 
         ran2=min(AM*iy,RNMX) 
         return 
       end function ran2 


    subroutine output_aux(aux, id, file_name, dir_name, itype, ntypes, np, lxyz)

      implicit none
      type(mesh_t), intent(in) :: aux(:,:)
      integer, intent(in) :: id, itype, ntypes, np, lxyz(:,:)
      character(3), intent(in) :: dir_name
      character(5), intent(in) :: file_name
      integer :: itypes, ip

      OPEN (UNIT=id,FILE=trim(dir_name//'/sc'//file_name//'.xyz'))
      do ip=1, np
         do itypes = itype, ntypes
            !if(aux(ip,itypes)%phi>=0.d0) then
               write(id,'(I10,I10,F10.2,I10)') lxyz(ip,1:2), aux(ip,itypes)%phi, itypes
            !end if
         end do
      end do
      close(id)
    end subroutine output_aux


    subroutine spherical_surface(Rc, R, ndim, delta)
      
      implicit none
      
      real, intent(in) :: Rc, delta
      integer, intent(inout) :: ndim
      integer, allocatable, intent(inout) :: R(:,:)
      real :: theta, phi, M_Pi
      character(len=10) :: dr, ds, sdx, sdy
      integer :: i, dx, dy, dz
      M_Pi = 3.14159265359
      ! Note:
      ! No intrinsic exists to convert between a numeric value and a formatted character 
      ! string representation  for instance, given the CHARACTER value '154',
      ! obtaining an INTEGER or REAL value with the value 154, or vice versa.
      ! Instead, this functionality is provided by internal-file I/O, 
      ! as in the following example: 
    
      ! Convert a string to a numeric value
      ! read (string,'(I10)') value
      ! print *, value
          
      ! Convert a value to a formatted string
      ! write (string2,'(I10)') value
      ! print *, string2


      i = 0
      theta = 0.0
      
      
      phi = 0.0
      do while (phi<=2.d0*M_Pi)
         dx =  int(anint(Rc * cos(phi)) )
         dy =  int(anint(Rc * sin(phi) ))
         
         write(sdx,'(I2)') abs(dx)
         write(sdy,'(I2)') abs(dy)
         
         
         dr = trim(sdx)//trim(sdy)
         
         if (dr .ne. ds) then
            i = i + 1
            
            R(i,1) = dx 
            R(i,2) = dy  
            
            ds = dr
         end if
         phi = phi + delta
      end do
      theta = theta + delta
      
      
      i = i +1
      R(i,1) = 0
      R(i,2) = 0
      
      ndim = i

    end subroutine spherical_surface
    


    subroutine gen_cell_points(Rc,R,ndim)
    
      implicit none
      
      real, intent(in) :: Rc
      integer, intent(inout) :: ndim
      integer, intent(inout) :: R(:,:)
      real :: dx, i, j,  M_Pi
      integer :: counter
      M_Pi = 3.14159265359
      counter = 0
      i = -Rc
      do while(i<=Rc)
         j =-Rc
         do while (j<=Rc)
  
               dx = sqrt(i**2 + j**2 )
               if(dx <= Rc) then
                  counter = counter + 1
                  R(counter,1) = int(anint(i))
                  R(counter,2) = int(anint(j))
                    
                  
               end if

            j = j +1.d0
         end do
         i = i + 1.d0
      end do

      ndim = counter

    end subroutine gen_cell_points



end module misc_m
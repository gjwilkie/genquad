program quad_test
use genquad, only: get_quadrature_rule
implicit none
  logical:: diagnostics_on = .true.
  integer:: N_scheme = 4
  integer:: N_integrand = 4
  integer:: N_res = 59
  real:: vmax = 4.0
  real:: pi

  real,dimension(:,:,:), allocatable:: errors, errors_gs2
  real,dimension(:,:), allocatable:: wgt_table_test, analytic_result
  real,dimension(:), allocatable:: v_table_test, abscissae, int_wgts, integrand
  real,dimension(:), allocatable:: gs2_wgts, gs2_vgrid, gs2_integrand, uniform_vgrid
  real,dimension(:), allocatable:: simpson_wgts, vgrid_overresolved, wgts_overresolved
  real:: delta_v , numerical_result, zero = 0.0, gs2_result, maxerr_v, maxerr_w
  integer::  iv,ischeme,ires,iint, Ntot
  integer, dimension(:),allocatable:: Nv_res
  real,external:: w1,w2,w3,w4
  real,dimension(10):: hermitenodes,hermiteweights

  allocate(Nv_res(N_res))

  pi = 3.14159265358979323846264338327950288419716939938
  
  do ires = 1,N_res
    Nv_res(ires) =ires+1
  end do

  open(unit=22,file="results.out", status="replace")
  open(unit=23,file="wgt1.plt", status="replace")
  open(unit=24,file="wgt2.plt", status="replace")
  open(unit=25,file="wgt3.plt", status="replace")
  open(unit=26,file="wgt4.plt", status="replace")
 
  allocate(analytic_result(N_scheme,N_integrand))
  allocate(abscissae(1:maxval(Nv_res(:))))
  allocate(int_wgts(1:maxval(Nv_res(:))))
  allocate(integrand(1:maxval(Nv_res(:))))
  allocate(gs2_integrand(1:maxval(Nv_res(:))))
  allocate(gs2_wgts(1:maxval(Nv_res(:))))
  allocate(gs2_vgrid(1:maxval(Nv_res(:))))
  allocate(errors(N_scheme,N_integrand,N_res))
  allocate(errors_gs2(N_scheme,N_integrand,N_res))

  ! Define analytic integrals
  ! first index is weight scheme
  ! second index is integrand
  ! Calculated from mathematica
  analytic_result(1,1) = 23.48501352449811
  analytic_result(1,2) = 0.03246362468013194
  analytic_result(1,3) = 0.801369616941368
  analytic_result(1,4) = 27.61348974269745
 
  analytic_result(2,1) = 619.0
  analytic_result(2,2) = 0.05882352941176441
  analytic_result(2,3) = 0.949590264857363
  analytic_result(2,4) = 12.64928118705395
  
  analytic_result(3,1) = 59.22503999731359
  analytic_result(3,2) = 3.39405260298183
  analytic_result(3,3) = 3.604857301284351
  analytic_result(3,4) = 58.20370297469304
   
  analytic_result(4,1) = 302.4775729802816
  analytic_result(4,2) = 0.02621488484812071
  analytic_result(4,3) = 1.267867002167052
  analytic_result(4,4) = 16.38784893315581
   
  do ischeme = 1,4
    write(22,*) "*****************************"
    write(22,'(A27,I2)') " Quadrature weight Scheme ", ischeme
    write(22,*) "*****************************"
    do ires = 1,N_res
      write(22,'(A5,I2)') "Nv = ", Nv_res(ires)
!      write(*,'(A5,I2)') "Nv = ", Nv_res(ires)
  
      if (ischeme .EQ. 1 ) then
        call get_quadrature_rule(w1,Nv_res(ires),0.0,1.0,abscissae,int_wgts,.true.,.true.)
      else if (ischeme .EQ. 2) then
        call get_quadrature_rule(w2,Nv_res(ires),0.0,1.0,abscissae,int_wgts,.false.,.true.)
      else if (ischeme .EQ. 3) then
        call get_quadrature_rule(w3,Nv_res(ires),-1.0,4.0,abscissae,int_wgts,.false.,.false.)
      else if (ischeme .EQ. 4) then
        call get_quadrature_rule(w4,Nv_res(ires),0.0,4.0,abscissae,int_wgts,.false.,.false.)
      end if

      write(22,*) "Calculated quadrature abscissae and weights:"
      do iv = 1,Nv_res(ires)
         write(22,*) abscissae(iv), int_wgts(iv)
      end do
      do iint = 1,N_integrand
        call set_integrand(iint,Nv_res(ires),abscissae(1:Nv_res(ires)),integrand(1:Nv_res(ires)))
        call perform_integral(Nv_res(ires),integrand(1:Nv_res(ires)),int_wgts(1:Nv_res(ires)),numerical_result)

        write(22,*) "Analytic result:  ", analytic_result(ischeme,iint)
        write(22,*) "Numerical result: ", numerical_result

        write(22,*) "---------------------"
        write(22,*) "Error = ", abs( numerical_result - analytic_result(ischeme,iint))/abs(numerical_result)
        write(22,*) "---------------------"
        write(22,*) " "
        errors(ischeme,iint,ires) = abs(numerical_result - analytic_result(ischeme,iint))/abs(analytic_result(ischeme,iint))
      end do
    end do
  end do

  close(22)

     do ires = 1,N_res
        write(23,*) Nv_res(ires), (errors(1,iint,ires),iint=1,N_integrand),(errors_gs2(1,iint,ires),iint=1,N_integrand)
        write(24,*) Nv_res(ires), (errors(2,iint,ires),iint=1,N_integrand),(errors_gs2(2,iint,ires),iint=1,N_integrand)
        write(25,*) Nv_res(ires), (errors(3,iint,ires),iint=1,N_integrand),(errors_gs2(3,iint,ires),iint=1,N_integrand)
        write(26,*) Nv_res(ires), (errors(4,iint,ires),iint=1,N_integrand),(errors_gs2(4,iint,ires),iint=1,N_integrand)
     end do

  close(23)
  close(24)
  close(25)
  close(26)

  open(unit=27,file="hermitecomp.dat", status="replace")
  call get_quadrature_rule(w1,20,0.0,1.0,abscissae,int_wgts,.true.,.true.)

  ! Actual Hermite quadrature rule for N=20, from Abramowitz and Stegun
  hermitenodes = (/ 0.2453407083009, 0.7374737285454, 1.2340762153953, 1.7385377121166, 2.2549740020893, 2.7888060584281, &
       3.3478545673832, 3.9447640401156, 4.6036824495507, 5.3874808900112  /)

  hermiteweights = (/4.622436696006e-1, 2.866755053628e-1, 1.090172060200e-1, 2.481052088746e-2, 3.243773342238e-3, & 
      2.283386360163e-4, 7.802556478532e-6, 1.086069370769e-7, 4.399340992273e-10, 2.229393645534e-13 /)

  maxerr_v = 0.0
  maxerr_w = 0.0
  do iv=1,10
    write(27,*) iv, " & ", abscissae(iv), " & ", -hermitenodes(10-iv+1), " & ", int_wgts(iv), " & ", hermiteweights(10-iv+1) &
      , " \\ "
    maxerr_v = max(maxerr_v, abs( (abscissae(iv)+hermitenodes(10-iv+1))/hermitenodes(10-iv+1)))
    maxerr_w = max(maxerr_w, abs( (int_wgts(iv)-hermiteweights(10-iv+1))/hermiteweights(10-iv+1)))
  end do
  do iv=11,20
    write(27,*) iv, " & ", abscissae(iv), " & ", hermitenodes(iv-10), " & ", int_wgts(iv), " & ", hermiteweights(iv-10)," \\ "
    maxerr_v = max(maxerr_v, abs( (abscissae(iv)-hermitenodes(iv-10))/hermitenodes(iv-10)))
    maxerr_w = max(maxerr_w, abs( (int_wgts(iv)-hermiteweights(iv-10))/hermiteweights(iv-10)))
  end do
  write(*,*) "Max error in Hermite abscissae = ", maxerr_v
  write(*,*) "Max error in Hermite weights = ", maxerr_w
    
  close(27)

contains

subroutine set_integrand(opt,N,v,f)
implicit none
  integer,intent(in):: opt, N
  real,dimension(1:), intent(in):: v
  real,dimension(1:), intent(inout):: f
  integer:: iv
!  real,external:: g1,g2,g3,g4

  do iv =1,N
    select case (opt)
      case(1) 
        f(iv) = g1(v(iv))
      case(2) 
        f(iv) = g2(v(iv))
      case(3) 
        f(iv) = g3(v(iv))
      case(4) 
        f(iv) = g4(v(iv))
    end select
  end do

end subroutine set_integrand


  ! Integrand 1:
  ! Polynomial
real function g1(v)
implicit none
  real,intent(in):: v
  g1 = 9.0 + v*(-10.0 + v*(7.0 + v*(-3.0 + v*(1.0 + v*5.0))))
  return
end function

  ! Integrand 2:
  ! Triognometric
real function g2(v)
implicit none
  real,intent(in):: v
  g2 = cos(4.0*v)
  return
end function

  ! Integrand 2:
  ! Hyperbolic
real function g3(v)
implicit none
  real,intent(in):: v
  g3 = tanh(2.0*v+1.0)
  return
end function
  ! Integrand 4:
  ! Exponential
real function g4(v)
implicit none
  real,intent(in):: v
  g4 = 10.0*exp(-v**4+2.0*v**2)
  return
end function

subroutine perform_integral(N,f,w,answer)
implicit none
   real,dimension(1:), intent(in):: f,w 
   real, intent(out):: answer
   integer, intent(in):: N
   real:: sum
   integer:: iv
  
   sum = 0.0
   do iv= 1,N
     sum = sum + w(iv)*f(iv)
   end do

   answer = sum

end subroutine perform_integral


end program

 
  ! Scheme 1: w(v) =  exp(-v)
  ! Exponential
real function w1(v)
implicit none
  real,intent(in):: v
  w1 = exp(-v**2)
  return
end function

real function w2(v)
implicit none
  real,intent(in):: v
  w2 = exp(-v)
  return
end function

real function w3(v)
implicit none
  real,intent(in):: v
  w3 = 1.0/(0.1+10.0*abs(v)**3)
!  w3 = 1.0
  return
end function

real function w4(v)
implicit none
  real,intent(in):: v
  w4 = 1.0/(1.0+v**2)
  return
end function



module depac_wrapper_test
       ! Note: the name of the module and the subroutine should not be the same (excluding _c)
       ! That's why I put _test at the end of module's name
 use, intrinsic :: iso_c_binding, only: c_char, c_int, c_float, c_null_char
 use LE_DryDepos_Gas_DEPAC, only : DryDepos_Gas_DEPAC
 implicit none
contains
 ! MODIFIED: Added new parameters to extract all resistance components and compensation points
 subroutine depac_wrapper_c(compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, &
                            lai, sai, nwet, lu, iratns, &
                            rc_tot, ccomp_tot, hlaw, react, &
                            status, c_ave_prev_nh3, ra, rb, catm, rc_eff, &
                            gsoil_eff_out, rsoil_eff_out, p, &
                            gw_out, gstom_out, cw_out, cstom_out, csoil_out) bind(C, name="depac_wrapper")
!     character(len=4), intent(in) :: compnam(*)
!     integer, intent(in) :: day_of_year, nwet, lu, iratns
!     real, intent(in) :: lat, t, ust, glrad, sinphi, rh, lai, sai, hlaw, react, c_ave_prev_nh3
!     real, intent(out) :: rc_tot, ccomp_tot
!     integer, intent(out) :: status

   character(kind=c_char), intent(in) :: compnam(*)

   ! Note: The "value" attribute in Fortran is used in interoperable subroutines (those with the bind(C) attribute) 
   ! to indicate that the argument should be passed by value rather than by reference.
   integer(c_int), value, intent(in) :: day_of_year
   integer(c_int), value, intent(in) :: nwet
   integer(c_int), value, intent(in) :: lu
   integer(c_int), value, intent(in) :: iratns 

   real(c_float), value, intent(in) :: lat
   real(c_float), value, intent(in) :: t
   real(c_float), value, intent(in) :: ust
   real(c_float), value, intent(in) :: glrad
   real(c_float), value, intent(in) :: sinphi
   real(c_float), value, intent(in) :: rh
   real(c_float), value, intent(in) :: lai
   real(c_float), value, intent(in) :: sai
   real(c_float), value, intent(in) :: hlaw
   real(c_float), value, intent(in) :: react
   real(c_float), value, intent(in) :: p                ! ADDED: pressure (Pa)
   !real(c_float), value, intent(in) :: tsea             ! sea surface temperature (K)
   !real(c_float), value, intent(in) :: smi              ! soil moisture index (there is still a
                                                         ! problem with the soil moisture over sea =0 
   real(c_float), value, intent(in) :: c_ave_prev_nh3
   !real(c_float), value, intent(in) :: c_ave_prev_so2
   real(c_float), value, intent(in) :: catm             ! actual atmospheric concentration (ug/m3)
   !real(c_float), value, intent(in) :: gamma_soil_water ! gammawater from gammawater.nc file
   real(c_float), value, intent(in) :: ra              ! aerodynamic resistance (s/m)
   real(c_float), value, intent(in) :: rb              ! boundary layer resistance (s/m)

   real(c_float), intent(out) :: rc_tot
   real(c_float), intent(out) :: ccomp_tot
   real(c_float), intent(out) :: rc_eff
   ! ADDED: Output parameters for soil conductance and resistance
   real(c_float), intent(out) :: gsoil_eff_out  ! Effective soil conductance (m/s)
   real(c_float), intent(out) :: rsoil_eff_out  ! Effective soil resistance (s/m)
   
   ! NEW: Added output parameters for additional resistances
   real(c_float), intent(out) :: gw_out         ! External leaf conductance (m/s)
   real(c_float), intent(out) :: gstom_out      ! Stomatal conductance (m/s)
   real(c_float), intent(out) :: cw_out         ! External leaf compensation point (ug/m3)
   real(c_float), intent(out) :: cstom_out      ! Stomatal compensation point (ug/m3)
   real(c_float), intent(out) :: csoil_out      ! Soil compensation point (ug/m3)
   
   integer(c_int), intent(out) :: status
   character(len=6) :: f_compnam

   call c_f_string(compnam, f_compnam)
   
   ! MODIFIED: Added all additional output parameters to get full details from DEPAC
   call DryDepos_Gas_DEPAC(f_compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, &
                           lai, sai, nwet, lu, iratns, &
                           rc_tot, ccomp_tot, hlaw, react, &
                           status, p=p, c_ave_prev_nh3=c_ave_prev_nh3, ra=ra, rb=rb, &
                           catm=catm, rc_eff=rc_eff, &
                           gw_out=gw_out, gstom_out=gstom_out, gsoil_eff_out=gsoil_eff_out, &
                           cw_out=cw_out, cstom_out=cstom_out, csoil_out=csoil_out)

   ! ADDED: Calculate effective soil resistance from conductance
   ! Use -9999.0 as error value when conductance is 0 or negative
   if (gsoil_eff_out > 0.0) then
     rsoil_eff_out = 1.0/gsoil_eff_out
   else
     rsoil_eff_out = -9999.0
   endif
   
   ! NOTE: We don't need to calculate cw and cstom here since 
   ! they are provided directly by DEPAC through cw_out and cstom_out
 end subroutine depac_wrapper_c

! This subroutine copies characters from a C-style null-terminated string (c_string)
! into a Fortran string (f_string), stopping when it encounters a null character or
! reaches the end of f_string.   
  subroutine c_f_string(c_string, f_string)
    character(kind=c_char), intent(in) :: c_string(*)
    character(len=*), intent(out) :: f_string
    integer :: i
    f_string = ''
    do i = 1, len(f_string)
      if (c_string(i) == c_null_char) exit
      f_string(i:i) = c_string(i)
    end do
  end subroutine c_f_string
end module depac_wrapper_test

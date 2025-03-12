module depac_wrapper_test
       ! Note: the name of the module and the subroutine should not be the same (excluding _c)
       ! That's why I put _test at the end of module's name
 ! Import necessary components from iso_c_binding for C-Fortran interoperability
 use, intrinsic :: iso_c_binding, only: c_char, c_int, c_float, c_null_char, c_bool
 use LE_DryDepos_Gas_DEPAC, only : DryDepos_Gas_DEPAC
 implicit none
contains
 ! MODIFIED: Added parameter to accept ccomp_tot as input when needed and to use override values
 ! This wrapper allows MicroHH to either use DEPAC's calculated compensation points
 ! or to override them with user-provided values.
 subroutine depac_wrapper_c(compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, &
                            lai, sai, nwet, lu, iratns, &
                            rc_tot, ccomp_tot, hlaw, react, &
                            status, c_ave_prev_nh3, ra, rb, catm, rc_eff, &
                            gsoil_eff_out, rsoil_eff_out, p, &
                            gw_out, gstom_out, cw_out, cstom_out, csoil_out, &
                            use_input_ccomp) bind(C, name="depac_wrapper")

   ! Input parameters: Component name (e.g., "NH3")
   character(kind=c_char), intent(in) :: compnam(*)

   ! Input parameters: Time and location
   integer(c_int), value, intent(in) :: day_of_year   ! Day of year (1-366)
   integer(c_int), value, intent(in) :: nwet          ! Surface wetness (0=dry, 1=wet, 9=snow)
   integer(c_int), value, intent(in) :: lu            ! Land use type
   integer(c_int), value, intent(in) :: iratns        ! NH3/SO2 ratio regime

   ! Input parameters: Meteorological and surface conditions
   real(c_float), value, intent(in) :: lat            ! Latitude (degrees)
   real(c_float), value, intent(in) :: t              ! Temperature (°C)
   real(c_float), value, intent(in) :: ust            ! Friction velocity (m/s)
   real(c_float), value, intent(in) :: glrad          ! Global radiation (W/m²)
   real(c_float), value, intent(in) :: sinphi         ! Sine of solar elevation angle
   real(c_float), value, intent(in) :: rh             ! Relative humidity (%)
   real(c_float), value, intent(in) :: lai            ! Leaf area index
   real(c_float), value, intent(in) :: sai            ! Surface area index
   real(c_float), value, intent(in) :: hlaw           ! Henry's law constant
   real(c_float), value, intent(in) :: react          ! Reactivity factor
   real(c_float), value, intent(in) :: p              ! Pressure (Pa)
   real(c_float), value, intent(in) :: c_ave_prev_nh3 ! Previous NH3 concentration
   real(c_float), value, intent(in) :: catm           ! Atmospheric concentration (µg/m³)
   real(c_float), value, intent(in) :: ra             ! Aerodynamic resistance (s/m)
   real(c_float), value, intent(in) :: rb             ! Boundary layer resistance (s/m)

   ! Output parameters: Resistances and compensation points
   real(c_float), intent(out) :: rc_tot               ! Total canopy resistance (s/m)
   real(c_float), intent(inout) :: ccomp_tot          ! Compensation point (µg/m³), can be input or output
   real(c_float), intent(out) :: rc_eff               ! Effective total resistance (s/m)
   real(c_float), intent(out) :: gsoil_eff_out        ! Effective soil conductance (m/s)
   real(c_float), intent(out) :: rsoil_eff_out        ! Effective soil resistance (s/m)
   real(c_float), intent(out) :: gw_out               ! External leaf conductance (m/s)
   real(c_float), intent(out) :: gstom_out            ! Stomatal conductance (m/s)
   real(c_float), intent(out) :: cw_out               ! External leaf compensation point (µg/m³)
   real(c_float), intent(out) :: cstom_out            ! Stomatal compensation point (µg/m³)
   real(c_float), intent(out) :: csoil_out            ! Soil compensation point (µg/m³)
   integer(c_int), intent(out) :: status              ! Status code (0 = success)

   ! Input flag: Whether to use the provided ccomp_tot value instead of calculating it
   integer(c_int), value, intent(in) :: use_input_ccomp ! 0 = calculate, non-zero = use input value

   ! Local variables
   character(len=6) :: f_compnam
   real(c_float) :: original_ccomp

   ! Convert C string to Fortran string
   call c_f_string(compnam, f_compnam)
   
   ! Store the input compensation point if we're going to use it
   if (use_input_ccomp /= 0) then
     original_ccomp = ccomp_tot
   end if
   
   ! Always call DEPAC to calculate all variables (resistances, conductances, etc.)
   call DryDepos_Gas_DEPAC(f_compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, &
                           lai, sai, nwet, lu, iratns, &
                           rc_tot, ccomp_tot, hlaw, react, &
                           status, p=p, c_ave_prev_nh3=c_ave_prev_nh3, ra=ra, rb=rb, &
                           catm=catm, rc_eff=rc_eff, &
                           gw_out=gw_out, gstom_out=gstom_out, gsoil_eff_out=gsoil_eff_out, &
                           cw_out=cw_out, cstom_out=cstom_out, csoil_out=csoil_out)

   ! If we're using the input compensation point, restore it and recalculate rc_eff
   if (use_input_ccomp /= 0) then
     ! Restore the original ccomp_tot value that was passed in
     ccomp_tot = original_ccomp
     
     ! Special case: If ccomp_tot is zero (or very close to zero), set rc_eff = rc_tot
     ! This is for the traditional deposition-only model with no compensation point
     if (abs(ccomp_tot) < 1.0e-12) then
       rc_eff = rc_tot
     else
       ! For non-zero compensation points, calculate rc_eff using the formula:
       ! rc_eff = ((ra + rb)*ccomp_tot + rc_tot*catm)/(catm-ccomp_tot)
       
       ! Check to prevent division by zero or negative values
       if (abs(catm - ccomp_tot) < 1.0e-12) then
         ! Surface concentration equal to atmospheric concentration - no exchange
         rc_eff = 9999999.0
       else
         rc_eff = ((ra + rb)*ccomp_tot + rc_tot*catm)/(catm-ccomp_tot)
       end if
     end if
   end if

   ! Calculate effective soil resistance from conductance
   ! Use -9999.0 as error value when conductance is 0 or negative
   if (gsoil_eff_out > 0.0) then
     rsoil_eff_out = 1.0/gsoil_eff_out
   else
     rsoil_eff_out = -9999.0
   endif
   
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

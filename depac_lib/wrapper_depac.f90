module depac_wrapper_test
       ! Note: the name of the module and the subroutine should not be the same (excluding _c)
       ! That's why I put _test at the end of module's name
 !==============================================================================
 ! DEPAC Wrapper for MicroHH - Explanation of Compensation Point Handling
 !==============================================================================
 ! This wrapper provides two modes of operation:
 !
 ! 1. Normal Mode (use_input_ccomp = 0):
 !    - DEPAC calculates all resistances, conductances, and compensation points
 !    - All calculated values are returned to MicroHH unchanged
 !
 ! 2. Override Mode (use_input_ccomp = 1):
 !    - The ccomp_tot value is provided by MicroHH (from ccomp_override_value)
 !    - All individual compensation points (cw_out, cstom_out, csoil_out) are 
 !      forced to ZERO regardless of the ccomp_tot value
 !    - This allows for custom bidirectional exchange modeling with controlled
 !      compensation points
 !
 ! The override mode is useful for:
 !   - Traditional deposition-only models (all comp. points = 0)
 !   - Testing specific compensation point values
 !   - Implementing alternative bidirectional exchange models
 !==============================================================================

 use, intrinsic :: iso_c_binding, only: c_char, c_int, c_float, c_null_char, c_bool
 use LE_DryDepos_Gas_DEPAC, only : DryDepos_Gas_DEPAC
 implicit none
contains
 ! MODIFIED: Completely revised wrapper to handle compensation point overrides properly
 subroutine depac_wrapper_c(compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, &
                            lai, sai, nwet, lu, iratns, &
                            rc_tot, ccomp_tot, hlaw, react, &
                            status, c_ave_prev_nh3, ra, rb, catm, rc_eff, &
                            gsoil_eff_out, rsoil_eff_out, p, &
                            gw_out, gstom_out, cw_out, cstom_out, csoil_out, &
                            use_input_ccomp, tsea, c_ave_prev_so2) bind(C, name="depac_wrapper")   ! CHANGED: added tsea, then c_ave_prev_so2 as the last argument (must match the extern "C" prototype in deposition.cxx)

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
   ! ADDED: water surface temperature (K). DEPAC reads tsea without a present() check for water
   !        land-use classes (le_drydepos_gas_depac.f90 l. 1753-1755). Before, it was not passed,
   !        so DEPAC read an absent optional argument (undefined behaviour, likely segfault).
   real(c_float), value, intent(in) :: tsea
   ! ADDED: long-term mean SO2 concentration C_ave,prev,SO2 [ug m-3]. DEPAC uses it only in the NH3 co-deposition
   !        factor of the external leaf resistance. The co-deposition block requires c_ave_prev_nh3 > 0 and c_ave_prev_so2 > 0,
   !        so c_ave_prev_so2 = 0 skips it and reproduces the results without this argument.
   real(c_float), value, intent(in) :: c_ave_prev_so2

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
   real(c_float) :: temp_cw_out, temp_cstom_out, temp_csoil_out
   logical, save :: so2_debug_done = .false.   ! ADDED: DEBUG (temporary). Print the long-term means once per rank
   integer :: env_status                       ! ADDED: DEBUG (temporary)

   ! Convert C string to Fortran string
   call c_f_string(compnam, f_compnam)

   ! ADDED: DEBUG (temporary, remove after verification). With the environment variable DEPAC_DEBUG_SO2 set,
   !        print the values that reach DEPAC on the first call of this rank.
   if (.not. so2_debug_done) then
     so2_debug_done = .true.
     call get_environment_variable("DEPAC_DEBUG_SO2", status=env_status)
     if (env_status == 0) then
       write(*,'(a,es12.5,a,es12.5,a,i0)') 'DEPAC_DEBUG_SO2 c_ave_prev_nh3=', c_ave_prev_nh3, &
                                          ' c_ave_prev_so2=', c_ave_prev_so2, ' iratns=', iratns
     end if
   end if
   
   ! Store the input compensation point if we're going to use it
   if (use_input_ccomp /= 0) then
     original_ccomp = ccomp_tot
   end if
   
   ! Always call DEPAC to calculate all variables (resistances, conductances, etc.)
   call DryDepos_Gas_DEPAC(f_compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, &
                          lai, sai, nwet, lu, iratns, &
                          rc_tot, ccomp_tot, hlaw, react, &
                          status, p=p, tsea=tsea, c_ave_prev_nh3=c_ave_prev_nh3, &
                          c_ave_prev_so2=c_ave_prev_so2, ra=ra, rb=rb, &   ! CHANGED: added c_ave_prev_so2=c_ave_prev_so2 (keyword argument, order irrelevant here)
                          catm=catm, rc_eff=rc_eff, &
                          gw_out=gw_out, gstom_out=gstom_out, gsoil_eff_out=gsoil_eff_out, &
                          cw_out=temp_cw_out, cstom_out=temp_cstom_out, csoil_out=temp_csoil_out)   ! CHANGED: added tsea=tsea

   ! CRITICAL CHANGE: Handle override mode more explicitly
   if (use_input_ccomp /= 0) then
     ! Restore the original ccomp_tot value that was passed in
     ccomp_tot = original_ccomp
     
     ! Force ALL individual compensation points to ZERO regardless of ccomp_tot value
     cw_out = 0.0
     cstom_out = 0.0
     csoil_out = 0.0
     
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
   else
     ! If not using override mode, use the values calculated by DEPAC
     cw_out = temp_cw_out
     cstom_out = temp_cstom_out
     csoil_out = temp_csoil_out
   end if

   ! Calculate effective soil resistance from conductance
   ! Use -9999.0 as error value when conductance is 0 or negative
   ! (deposition.cxx now treats rsoil_eff_out < 0 as "no soil path", i.e. v_e = 0)
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

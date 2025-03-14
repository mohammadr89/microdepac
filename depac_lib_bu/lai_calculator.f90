module lai_calculator
  implicit none

  type laitype
    integer     :: sgs50, egs50, s_lai_len, e_lai_len
    real        :: dsgs, degs, laimin, laimax
  end type laitype

  type(laitype), parameter :: lai_par(18) = (/ &
    laitype(0, 366, 140, 135, 0.0, 0.0, 2.0, 3.5), &  ! 1. grass
    laitype(130, 250, 35, 65, 0.0, 0.0, 0.0, 4.2), &  ! 2. arable land
    laitype(130, 250, 35, 65, 0.0, 0.0, 0.0, 4.2), &  ! 3. permanent crops
    laitype(0, 366, 1, 1, 0.0, 0.0, 5.0, 5.0), &      ! 4. coniferous forest
    laitype(100, 307, 20, 30, 1.5, -2.0, 0.0, 4.0), & ! 5. deciduous forest
    laitype(-999, -999, -999, -999, -999.0, -999.0, -999.0, -999.0), & ! 6. water
    laitype(-999, -999, -999, -999, -999.0, -999.0, -999.0, -999.0), & ! 7. urban
    laitype(0, 366, 140, 135, 0.0, 0.0, 2.0, 3.5), &  ! 8. other
    laitype(-999, -999, -999, -999, -999.0, -999.0, -999.0, -999.0), & ! 9. desert
    laitype(-999, -999, -999, -999, -999.0, -999.0, -999.0, -999.0), & ! 10. ice
    laitype(0, 366, 140, 135, 0.0, 0.0, 2.0, 3.5), &  ! 11. savanna
    laitype(100, 307, 20, 30, 1.5, -2.0, 0.0, 4.0), & ! 12. tropical forest
    laitype(-999, -999, -999, -999, -999.0, -999.0, -999.0, -999.0), & ! 13. water_inland
    laitype(130, 250, 35, 65, 0.0, 0.0, 0.0, 4.2), &  ! 14. mediterrean scrub
    laitype(0, 366, 140, 135, 0.0, 0.0, 2.0, 3.5), &  ! 15. semi-natural grassland
    laitype(130, 250, 35, 65, 0.0, 0.0, 0.0, 4.2), &  ! 16. wheat
    laitype(100, 307, 20, 30, 1.5, -2.0, 0.0, 4.0), & ! 17. beech
    laitype(0, 366, 1, 1, 0.0, 0.0, 5.0, 5.0) &       ! 18. spruce
  /)

contains

  subroutine calc_lai(day_of_year, lat, lu, lai, sai, laimax)
    integer, intent(in) :: day_of_year, lu
    real, intent(in) :: lat
    real, intent(out) :: lai, sai, laimax
    
    integer :: sgs, egs
    type(laitype) :: lai_params

    lai_params = lai_par(lu)

    sgs = nint(lai_params%sgs50 + lai_params%dsgs * (lat-50.0))
    egs = nint(lai_params%egs50 + lai_params%degs * (lat-50.0))

    if (day_of_year < sgs .or. day_of_year > egs) then
        lai = 0.0
    else if (day_of_year <= sgs + lai_params%s_lai_len) then
        lai = lai_params%laimin + (lai_params%laimax-lai_params%laimin) * &
              (day_of_year-sgs) / lai_params%s_lai_len
    else if (day_of_year >= egs - lai_params%e_lai_len) then
        lai = lai_params%laimin + (lai_params%laimax-lai_params%laimin) * &
              (egs-day_of_year) / lai_params%e_lai_len
    else
        lai = lai_params%laimax
    end if

    ! Calculate SAI based on land use type
    select case(lu)
    case(4, 5, 17, 18)  ! Coniferous forest, Deciduous forest, Beech, Spruce
        sai = lai + 1.0
    case(3)  ! Permanent crops
        sai = lai + 0.5
    case(2, 16)  ! Arable land, Wheat
        if (day_of_year < sgs .or. day_of_year > egs) then
            sai = lai + 0.5
        else
            sai = lai + min(1.5, lai/lai_params%laimax * 1.5)
        end if
    case default
        sai = lai
    end select

    laimax = lai_params%laimax
  end subroutine calc_lai

end module lai_calculator

module lai_wrapper_mod
    use lai_calculator
    implicit none
contains
    subroutine calc_lai_wrapper(day_of_year, lat, lu, lai, sai, laimax) bind(C, name="calc_lai_wrapper")
        integer, value, intent(in) :: day_of_year, lu
        real, value, intent(in) :: lat
        real, intent(out) :: lai, sai, laimax
        
        call calc_lai(day_of_year, lat, lu, lai, sai, laimax)
    end subroutine calc_lai_wrapper
end module lai_wrapper_mod

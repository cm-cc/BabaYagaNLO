program momenta_correction

    use recola
    implicit none
    real(kind=8) :: psp(4,4),r1(2), r2(2), mm

    call use_gfermi_scheme_rcl(a=7.5552553808086852D-003)
    call set_pole_mass_muon_rcl(20.565799999989277d0, 0d0)
    call set_light_fermions_rcl(1d-3)


    call define_process_rcl(1, 'A A -> tau+ tau-', 'LO')

    call set_dynamic_settings_rcl(1)

    call generate_processes_rcl()

    psp(:,1) = [ 69.885301257628583d0,  0d0, 0d0, 69.885301257628583d0]
    psp(:,2) = [ 69.885301257628726d0, 0d0, 0d0, -69.885301257628726d0]
    psp(:,3) = [ 69.885301257628711d0,   27.694166675418500d0,   22.29490497385732d0,   60.165818171369821d0]
    psp(:,4) = [ 69.885301257628612d0,  -27.694166675418500d0,  -22.29490497385732d0,  -60.165818171369779d0]

    call compute_process_rcl(1,psp,'LO', r1)

    psp(:,1) = [  12.64150249837562d0, 0d0, 0d0,  12.64150249837562d0]
    psp(:,2) = [ 386.34294716922125d0, 0d0, 0d0, -386.34294716922125d0]
    psp(:,3) = [  38.62826083972509d0,  27.694166675418500d0,  22.29490497385732d0,  -15.10340623378218d0]
    psp(:,4) = [ 360.35618882787185d0, -27.694166675418500d0, -22.29490497385732d0, -358.59803843706345d0]

    call compute_process_rcl(1,psp,'LO', r2)

    call check_value(r1(1), r2(1), 5)

    call reset_recola_rcl()

    contains

    subroutine check_value(value1, value2, tolerance)
        real(kind=8), intent(in) :: value1, value2
        integer, intent(in) :: tolerance
        real(kind=8) :: diff, expected

        diff = abs(value1 - value2) / max(abs(value1), abs(value2))
        expected = 10d0**(-tolerance)
        if ( diff > expected) then
            write(*,*) 'Values do not match within tolerance: ', value1, value2
            stop 9
        else
            write(*,*) 'Values match within tolerance: ', value1, value2
        end if
    end subroutine check_value
end program momenta_correction

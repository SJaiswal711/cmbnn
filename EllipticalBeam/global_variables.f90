module global_variables
    use math_library
    implicit none
    ! Declare global variables
    character(len=100) :: filename
    real(kind=8) :: theta1, theta2, w1, w2, w3
    real(kind=8) :: sigma_x, sigma_y, sigma, fwhm_x, fwhm_y
    integer :: nside, npix
    real :: scan_time
    integer :: steps
    real(kind=8), allocatable :: time_periods(:)
    !integer :: nmaps
    !real(8) :: nullval
    !logical :: anynull

contains
    subroutine initialize_globals()
        real(kind=8), parameter :: pi = 3.141592653589793d0

        theta1 = 7.5 * pi / 180.0d0
        theta2 = 85.0 * pi / 180.0d0
        w1 = 2.0d0 * pi  ! rad/min
        w3 = 0.000011954d0  ! rad/min
		w2 = 2.0d0 * w3  ! rad/min
		
        nside = 128
        npix = 12 * nside**2

        fwhm_x = 10.0 / 60.0 * pi / 180.0
        fwhm_y = 15.0 / 60.0 * pi / 180.0

        sigma_x = fwhm_x / sqrt(8.0 * log(2.0))
        sigma_y = fwhm_y / sqrt(8.0 * log(2.0))
        sigma = max(sigma_x, sigma_y)
        
        scan_time = sqrt(4 * pi / npix) / w1
	filename = "map1.fits"
    end subroutine initialize_globals
    
    subroutine setup_time_periods(scan_time)
        real, intent(in) :: scan_time
        real(kind=8) :: start_time, duration
        integer :: i

        ! Define time step parameters
        start_time = 5.0d0
        duration = 1 ! Duration in minutes (one year)

        ! Compute number of steps
        steps = int(duration / scan_time)-1
        !steps = 1
        print*,"steps: ",steps
        allocate(time_periods(steps))
        call linspace(time_periods, start_time, (duration+start_time), steps)
        print*, "Allocated memory for time_period"
        
    end subroutine setup_time_periods
    
    !subroutine load_map()
    !   use fitstools
    !    implicit none
    !    nmaps = 1
    !    nullval =0
    !    anynull = .false.
    !    filename = "map1.fits"
    !    allocate(map(0:npix-1, 1:nmaps)) ! Allocate map as a 1D array with size npix
        ! Ensure read_bintab is correctly linked or available
        !call read_bintab(filename, map, npix, 1, nullval, anynull)
    !    call read_bintab(filename, map, npix, nmaps, nullval, anynull)
    !    print*, "Map is loaded"
    !    print*,npix, nmaps, nullval, anynull,"npix, nmaps, nullval, anynull"
    !end subroutine load_map

end module global_variables

subroutine get_vectors(t, result_r, result_s)
  use global_variables
  implicit none
  real(kind=8), intent(in) :: t
  real(kind=8), intent(out) :: result_r(3,1), result_s(3,1)
 
  !real, parameter :: pi = 3.141592653589793d0
  real(kind=8) :: cos_w1t, sin_w1t, cos_w2t, sin_w2t
  real(kind=8) :: cos_theta1, sin_theta1, cos_theta2, sin_theta2
  !real(kind=8) :: theta1, theta2, w1, w2, w3
  real(kind=8), dimension(3, 3) :: A, B, C,result1
  real(kind=8), dimension(3,1) :: D_R, D_S
   call initialize_globals()
  ! Calculate trigonometric terms
  cos_w1t = cos(w1 * t)
  sin_w1t = sin(w1 * t)
  
  cos_w2t = cos(w2 * t)
  sin_w2t = sin(w2 * t)
  
  cos_theta1 = cos(theta1)
  sin_theta1 = sin(theta1)
  
  cos_theta2 = cos(theta2)
  sin_theta2 = sin(theta2)

  ! Define the matrices
  A = reshape([cos_w3t, sin_w3t, 0.0d0, -sin_w3t, cos_w3t, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3, 3], order=[2,1])
  B = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, cos_w2t, sin_w2t, 0.0d0, -sin_w2t, cos_w2t], [3, 3], order=[2,1])
  C = reshape([cos_theta1, 0.0d0, sin_theta1, 0.0d0, 1.0d0, 0.0d0, -sin_theta1, 0.0d0, cos_theta1], [3, 3], order=[2,1])

  ! Define vectors
  D_R = reshape([cos_theta2, sin_theta2 * cosw1t, sin_theta2 * sinw1t], [3, 1], order=[2,1])
  D_S = reshape([1.0d0, 0.0d0, 0.0d0], [3, 1], order=[2,1])

  ! Perform matrix multiplications (optimized)
  result1 = matmul(matmul(A,B),C)
  result_r = matmul(result1, D_R)
  result_s = matmul(result1, D_S)
  !print*,"result1",result1
  !print*, "D_R = ", D_R
  !print*,"result_r",result_r
  !print*,"result_s",result_s
end subroutine get_vectors


subroutine anglev(vec1, vec2, angle)
    use math_library
    implicit none
    real(8), intent(in) :: vec1(3), vec2(3)
    real(8), intent(out) :: angle
    real(8) :: dot_prod, clipped_dp
    integer :: i=0

    dot_prod = dot_product(vec1, vec2)

    ! Clip dot_product to the valid range for arccos to avoid NaNs
    clipped_dp = max(-1.0d0, min(1.0d0, dot_prod))
    !print*,clipped_dp

    ! Calculate the angle
    angle = acos(clipped_dp)
end subroutine anglev


subroutine process_time_step(input_map, time_step, pix_ring, convolved_temperature)
    use math_library 
    use global_variables
    use healpix_modules
    implicit none
    !real(8), dimension(:,:), intent(in) :: input_map
    real(8), dimension(:,:), intent(in) :: input_map
    real(8), intent(in) :: time_step
    integer, intent(out) :: pix_ring
    real(8), intent(out) :: convolved_temperature
    
    integer, allocatable :: listpix(:)
    integer, parameter :: length = 20
    integer :: nlist, i

    real(8) :: x(length)
    real(8) :: y(length)

    real(8) :: t
    real(8) :: R(3), S(3), R_i(3)
    real(8) :: Z_t(3), I_t(3), N_t(3),nor_A_i(3), nor_N_t(3)
    real(8) :: Rc(3)
    real(8) :: theta_i, alpha_i, A_i(3)
    real(8), dimension(length) :: neighbor_temperatures
    real(8) :: exp_term
    allocate(listpix(length))

    call initialize_globals()
    t=time_step
    !print*,map
    ! 1. Calculate R(t) and S(t) vectors
    call get_vectors(t, R, S)
 
    ! 2. Calculate pixel number along R(t) vector (ring format)
    call vec2pix_ring(nside, R, pix_ring)
    !print*,"pix,R: ",pix_ring,R
    ! 3. Calculate Z_t, I_t, and N_t
    call cross_product(R, S, Z_t)
    call cross_product(R, Z_t, I_t)
    N_t = I_t

    ! 4. Find neighboring pixels in RING format
    call pix2vec_ring(nside, pix_ring, Rc)

    !print*,Rc,nside
    call query_disc(nside, Rc, 3.0d0 * sigma, listpix, nlist, nest=0, inclusive=1)
    !print*, "lp",listpix
    ! 5. Angular separation between central pixel and neighboring pixels
    do i = 1, nlist
        call pix2vec_ring(nside, listpix(i), R_i)
        
        !print*, "R_i = ",R_i
        call anglev(Rc, R_i, theta_i)

        ! 6. A_i = line joining central pixel and neighbor pixel
        A_i = Rc - R_i

        ! 7. Angle between N_t and A_i
        call normalize_vector(A_i, nor_A_i)
        call normalize_vector(N_t, nor_N_t)
        call anglev(nor_A_i, nor_N_t, alpha_i)

        ! 8. x_i and y_i
        x(i) = theta_i * cos(alpha_i)
        y(i) = theta_i * sin(alpha_i)
        !print*,"x(i),y(i) : ",x(i),y(i)
        
        !print*, "x_i = ",x(i)
        !print*, "y_i = ",y(i)
        ! 10. Retrieve the temperature of the neighboring pixel
        neighbor_temperatures(i) = input_map(listpix(i), 1)
        !print*,"neighbor_temperatures(i) : ",neighbor_temperatures(i)
        !print*,"neighbor_temperatures(i)=",listpix(i)
        
    end do
    print*,time_step, pix_ring,"neighbor_temperatures(i) : "
    ! 11. Apply elliptical convolution
    convolved_temperature = 0.0d0
    do i = 1, nlist
        ! Split the long line into multiple lines for readability
        exp_term = exp(-x(i)**2 / (2.0d0 * sigma_x**2) - y(i)**2 / (2.0d0 * sigma_y**2))
        convolved_temperature = convolved_temperature + neighbor_temperatures(i) * exp_term
    end do
end subroutine process_time_step

program main
    use global_variables
    use math_library
    use fitstools
    implicit none

    integer :: i, pix_ring,j
    real(kind=8) :: time_step, convolved_temperature
    real(kind=8) :: t
    real(8) :: nullval
    logical :: anynull
    integer :: nmaps
    real(8), dimension(:,:), allocatable :: map
    real(8) :: start_time, end_time, elapsed_time
    integer :: ios
    integer :: unit
    integer, parameter :: output_unit = 10  ! File unit number
        !print*,"1.pix_ring: ",pix_ring
    ! Explicit interface for process_time_step
    interface
        subroutine process_time_step(input_map, time_step, pix_ring, convolved_temperature)
            use math_library
            use global_variables
            use healpix_modules
            implicit none
            real(8), dimension(:,:), intent(in) :: input_map
            real(8), intent(in) :: time_step
            integer, intent(out) :: pix_ring
            real(8), intent(out) :: convolved_temperature
        end subroutine process_time_step
    end interface

    ! Initialize global variables
    call initialize_globals()
    
    nmaps = 1
    nullval = 0
    anynull = .false.

    allocate(map(0:npix-1, 1:nmaps)) ! Allocate map as a 1D array with size npix

    call input_map(filename, map, npix, nmaps)
    open(unit=unit, file='map.dat', status='replace', action='write')

    ! Write the data to the file
    do i = 0, npix-1
        ! Write index and values in a single line
        write(unit, '(I4, 1X, F10.5, A)', advance='no') i, map(i, 1), ', '
        do j = 2, nmaps
            write(unit, '(F10.5, A)', advance='no') map(i, j), ', '
        end do
        write(unit, *)  ! End the line
    end do

    ! Close the file
    close(unit)
    ! Set up the time periods based on scan time
    call setup_time_periods(scan_time)
    
    !call process_time_step(map, 501.0d-2, pix_ring, convolved_temperature)
    !print*,"pix_ring, convolved_temperature: ",pix_ring, convolved_temperature


! Open file for writing
    open(unit=output_unit, file='convolved_map.dat', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file for writing.'
        stop
    end if

    ! Record start time
    call cpu_time(start_time)
    
    ! Assume you want to process the map for each time step in `time_periods`
    do i = 1, steps
        t = time_periods(i)
        call process_time_step(map, t, pix_ring, convolved_temperature)
        
        ! Write the results to the file
        write(output_unit, '(F12.6, I10, F12.6)') t, pix_ring, convolved_temperature
    end do

    ! Record end time
    call cpu_time(end_time)
    
    ! Calculate and print elapsed time
    elapsed_time = end_time - start_time

    ! Write elapsed time to the file
    print*, 'Total execution time (seconds):', elapsed_time
    
    ! Close the file
    close(output_unit)

    ! Deallocate map
    deallocate(map)

end program main

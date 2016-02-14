! The Computer Language Benchmarks Game
! http://benchmarksgame.alioth.debian.org/
!
! Contributed by Jannis Teunissen. Adapted from Fortran versions by Jason
! Blevins, George R. Gonzalez and Simon Geard
!
! ifort -fast -openmp -o mandelbrot mandelbrot.f90
program mandelbrot
  implicit none

  integer, parameter         :: dp     = selected_real_kind(15, 307)
  integer, parameter         :: int8   = selected_int_kind(2)
  integer, parameter         :: iter   = 50
  real(dp), parameter        :: limit2 = 4.0_dp
  character(len=40)          :: tmp_string
  integer                    :: n_pixels, i, j, n
  real(dp)                   :: Cr(8), Ci, dx
  integer(int8), allocatable :: output_buf(:,:)
  character, allocatable :: test(:)


  ! Get number of pixels
  call get_command_argument(1, tmp_string)
  read(tmp_string, *) n_pixels
  if (modulo(n_pixels, 8) /= 0) stop "Argument should be multiple of 8"

  allocate(output_buf(n_pixels/8, n_pixels))
  dx = 2.0_dp / n_pixels        ! Precalculate constants

  !$OMP parallel do private(y, x, i, Cr, Ci)
  do j = 0, n_pixels - 1
     Ci = dx * j - 1.0_dp
     do i = 0, n_pixels/8 - 1
        Cr = dx * [(8 * i + n, n = 0, 7)] - 1.5_dp
        output_buf(i+1, j+1) = mandel_byte(Cr, Ci, iter)
     end do
  end do
  !$OMP end parallel do

  write(*, '("P4",/,i0," ",i0)') n_pixels, n_pixels ! PBM header
  allocate(test(n_pixels**2 / 8))
  write(*, *) transfer(output_buf, test)

contains

  function mandel_byte(Cr, Ci, iter) result(byte)
    real(dp), intent(in) :: Cr(8), Ci
    integer, intent(in)  :: iter
    integer(int8)        :: byte
    logical              :: lbyte(8)
    real(dp)             :: Zi(8), Zr(8), Ti(8), Tr(8)
    integer              :: i

    Zr = 0; Zi = 0; Tr = 0; Ti = 0

    do i = 1, iter
       Zi = 2 * Zr * Zi + Ci
       Zr = Tr - Ti + Cr
       Ti = Zi * Zi
       Tr = Zr * Zr
       if (all(Tr + Ti > limit2)) exit
    end do

    byte  = 0_int8
    lbyte = (Tr + Ti <= limit2)
    do i = 1, 8
       if (lbyte(i)) byte = ibset(byte, 8-i)
    end do
  end function mandel_byte

end program

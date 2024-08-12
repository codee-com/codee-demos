module time_module

   implicit none

   integer :: nb_ticks_sec, nb_ticks_max

contains

   subroutine time_inizialization()
      call SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
   end subroutine time_inizialization

   subroutine get_time(nb_ticks)
      integer, intent(inout) :: nb_ticks

      call SYSTEM_CLOCK(COUNT=nb_ticks)
   end subroutine get_time

   subroutine get_real_time(nb_ticks_initial, nb_ticks_final, elapsed_time)
      integer, intent(in) :: nb_ticks_initial, nb_ticks_final
      double precision, intent(out) :: elapsed_time
      integer :: nb_ticks

      nb_ticks = nb_ticks_final - nb_ticks_initial
      if (nb_ticks_final < nb_ticks_initial) &
         nb_ticks = nb_ticks + nb_ticks_max
      elapsed_time = REAL(nb_ticks)/(nb_ticks_sec)

   end subroutine get_real_time

end module time_module

program print_emoji
  implicit none
  character(len=:), allocatable :: emoji_string

  ! Assign a string containing emojis
  emoji_string = 'Fortran is 💪, 😎, 🔥!'

  ! Print the string to the console
  write(*,*) emoji_string

end program print_emoji

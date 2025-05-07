module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    real*8, intent(in) :: inTgas(cellsNumber)
    real*8 :: Tgas, T, invT
    integer :: icell
    real*8::n(cellsNumber,speciesNumber)

    call computeTotalDensityExcludingM(n)
  
  
    do icell = 1, cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0 / Tgas

      !O + O2 + M -> O3 + M
      krate(icell,1) = 3.11d-34*(T/298)**(-2.0)

      !O + O3 -> O2 + O2
      krate(icell,2) = 1.83d-11*exp(-2164/T)

      !H + HO2 -> O + H2O
      krate(icell,3) = 1.6d-12

      !O + CO + M -> CO2 + M
      krate(icell,4) = 1.7d-33*exp(-1510/T)

      !H + O2 + M -> HO2 + M
      krate(icell,5) = ((4.4d-32*(T/300)**(-1.3d0))*n(icell, patmo_idx_M)/(1d0+((4.4d-32*(T/300)**(-1.3d0))*n(icell, patmo_idx_M)/(7.5d-11*(T/300)**(2d-1))))*6d-1**((1d0+(log10((4.4d-32*(T/300)**(-1.3d0))*n(icell, patmo_idx_M)/(7.5d-11*(T/300)**(2d-1))))**(2d0))**(-1d0)))/n(icell, patmo_idx_M)
      
    end do

  end subroutine computeRates

  subroutine computeTotalDensityExcludingM(n)
    use patmo_commons
    use patmo_parameters
    implicit none

    real*8, intent(inout) :: n(:,:)  ! n(cellsNumber, speciesNumber)
    integer :: j

    n(:, patmo_idx_M) = 0.0d0
    do j = 1, chemSpeciesNumber
      if (j /= patmo_idx_M) then
        n(:, patmo_idx_M) = n(:, patmo_idx_M) + nAll(:, j)
      end if
    end do
  end subroutine computeTotalDensityExcludingM

end module patmo_rates

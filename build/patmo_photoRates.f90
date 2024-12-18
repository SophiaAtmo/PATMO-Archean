module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !O2 -> O + O
    krate(:,212) = integrateXsec(1, tau(:,:))

    !O3 -> O2 + O(1D)
    krate(:,213) = integrateXsec(2, tau(:,:))

    !O3 -> O2 + O
    krate(:,214) = integrateXsec(3, tau(:,:))

    !HO2 -> OH + O
    krate(:,215) = integrateXsec(4, tau(:,:))

    !H2O -> OH + H
    krate(:,216) = integrateXsec(5, tau(:,:))

    !H2O -> H2 + O
    krate(:,217) = integrateXsec(6, tau(:,:))

    !N2O -> N2 + O(1D)
    krate(:,218) = integrateXsec(7, tau(:,:))

    !NO2 -> NO + O
    krate(:,219) = integrateXsec(8, tau(:,:))

    !NO3 -> NO + O2
    krate(:,220) = integrateXsec(9, tau(:,:))

    !NO3 -> O + NO2
    krate(:,221) = integrateXsec(10, tau(:,:))

    !N2O5 -> NO2 + NO3
    krate(:,222) = integrateXsec(11, tau(:,:))

    !N2O5 -> O + NO + NO3
    krate(:,223) = integrateXsec(12, tau(:,:))

    !HNO3 -> OH + NO2
    krate(:,224) = integrateXsec(13, tau(:,:))

    !HNO3 -> H + NO3
    krate(:,225) = integrateXsec(14, tau(:,:))

    !CH4 -> CH3 + H
    krate(:,226) = integrateXsec(15, tau(:,:))

    !CH3OOH -> CH3O + OH
    krate(:,227) = integrateXsec(16, tau(:,:))

    !CH2O -> H + CHO
    krate(:,228) = integrateXsec(17, tau(:,:))

    !CH2O -> H2 + CO
    krate(:,229) = integrateXsec(18, tau(:,:))

    !CO2 -> CO + O
    krate(:,230) = integrateXsec(19, tau(:,:))

    !H2O2 -> OH + OH
    krate(:,231) = integrateXsec(20, tau(:,:))

    !H2O2 -> H + HO2
    krate(:,232) = integrateXsec(21, tau(:,:))

    !COS -> CO + S
    krate(:,233) = integrateXsec(22, tau(:,:))

    !SO -> S + O
    krate(:,234) = integrateXsec(23, tau(:,:))

    !CS2 -> CS + S
    krate(:,235) = integrateXsec(24, tau(:,:))

    !H2S -> SH + H
    krate(:,236) = integrateXsec(25, tau(:,:))

    !SO2 -> SO + O
    krate(:,237) = integrateXsec(26, tau(:,:))

    !SO3 -> SO2 + O
    krate(:,238) = integrateXsec(27, tau(:,:))

    !H2SO4 -> SO2 + OH + OH
    krate(:,239) = integrateXsec(28, tau(:,:))

    !S2O2 -> SO + SO
    krate(:,240) = integrateXsec(29, tau(:,:))

    !S2O -> SO + S
    krate(:,241) = integrateXsec(30, tau(:,:))

    !N2H4 -> H + N2H3
    krate(:,242) = integrateXsec(31, tau(:,:))

    !NH3 -> H + NH2
    krate(:,243) = integrateXsec(32, tau(:,:))

    !NH3 -> H2 + NH
    krate(:,244) = integrateXsec(33, tau(:,:))

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber), dE
    integer::j

    ! !loop on cells (stride photobins)
    ! do j=1,cellsNumber
    !    integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:) &
        !         /energyMid(:)*energySpan(:)*exp(-tau(:,j))) / planck_eV
    ! end do

    dE = 0.1!nm
    !loop on cells (stride photobins)
    do j=1,cellsNumber
      integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:)*exp(-2*tau(:,j))*dE)
    end do

  end function integrateXsec

end module patmo_photoRates

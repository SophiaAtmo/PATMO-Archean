module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !O2 -> O + O
    krate(:,215) = integrateXsec(1, tau(:,:))

    !O3 -> O2 + O(1D)
    krate(:,216) = integrateXsec(2, tau(:,:))

    !O3 -> O2 + O
    krate(:,217) = integrateXsec(3, tau(:,:))

    !OH -> O + H
    krate(:,218) = integrateXsec(4, tau(:,:))

    !OH -> O(1D) + H
    krate(:,219) = integrateXsec(5, tau(:,:))

    !HO2 -> OH + O
    krate(:,220) = integrateXsec(6, tau(:,:))

    !H2O -> OH + H
    krate(:,221) = integrateXsec(7, tau(:,:))

    !H2O -> H2 + O
    krate(:,222) = integrateXsec(8, tau(:,:))

    !H2 -> H + H
    krate(:,223) = integrateXsec(9, tau(:,:))

    !N2O -> N2 + O(1D)
    krate(:,224) = integrateXsec(10, tau(:,:))

    !NO2 -> NO + O
    krate(:,225) = integrateXsec(11, tau(:,:))

    !NO3 -> NO + O2
    krate(:,226) = integrateXsec(12, tau(:,:))

    !NO3 -> O + NO2
    krate(:,227) = integrateXsec(13, tau(:,:))

    !N2O5 -> NO2 + NO3
    krate(:,228) = integrateXsec(14, tau(:,:))

    !N2O5 -> O + NO + NO3
    krate(:,229) = integrateXsec(15, tau(:,:))

    !HNO3 -> OH + NO2
    krate(:,230) = integrateXsec(16, tau(:,:))

    !HNO3 -> H + NO3
    krate(:,231) = integrateXsec(17, tau(:,:))

    !CH4 -> CH3 + H
    krate(:,232) = integrateXsec(18, tau(:,:))

    !CH3OOH -> CH3O + OH
    krate(:,233) = integrateXsec(19, tau(:,:))

    !CH2O -> H + CHO
    krate(:,234) = integrateXsec(20, tau(:,:))

    !CH2O -> H2 + CO
    krate(:,235) = integrateXsec(21, tau(:,:))

    !CHO -> H + CO
    krate(:,236) = integrateXsec(22, tau(:,:))

    !CO2 -> CO + O
    krate(:,237) = integrateXsec(23, tau(:,:))

    !H2O2 -> OH + OH
    krate(:,238) = integrateXsec(24, tau(:,:))

    !H2O2 -> H + HO2
    krate(:,239) = integrateXsec(25, tau(:,:))

    !COS -> CO + S
    krate(:,240) = integrateXsec(26, tau(:,:))

    !SO -> S + O
    krate(:,241) = integrateXsec(27, tau(:,:))

    !CS2 -> CS + S
    krate(:,242) = integrateXsec(28, tau(:,:))

    !H2S -> SH + H
    krate(:,243) = integrateXsec(29, tau(:,:))

    !SO2 -> SO + O
    krate(:,244) = integrateXsec(30, tau(:,:))

    !SO3 -> SO2 + O
    krate(:,245) = integrateXsec(31, tau(:,:))

    !H2SO4 -> SO2 + OH + OH
    krate(:,246) = integrateXsec(32, tau(:,:))

    !CH3OH -> CH3 + OH
    krate(:,247) = integrateXsec(33, tau(:,:))

    !CH3OH -> CH3O + H
    krate(:,248) = integrateXsec(34, tau(:,:))

    !S2O2 -> SO + SO
    krate(:,249) = integrateXsec(35, tau(:,:))

    !S2O -> SO + S
    krate(:,250) = integrateXsec(36, tau(:,:))

    !N2H4 -> H + N2H3
    krate(:,251) = integrateXsec(37, tau(:,:))

    !NH3 -> H + NH2
    krate(:,252) = integrateXsec(38, tau(:,:))

    !NH3 -> H2 + NH
    krate(:,253) = integrateXsec(39, tau(:,:))

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

    !dE = (wavelengMax-wavelengMin)/photoBinsNumber (nm)
    dE = 0.1

    !loop on cells (stride photobins)
    do j=1,cellsNumber
      integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:)*exp(-2*tau(:,j))*dE)
    end do

  end function integrateXsec

end module patmo_photoRates

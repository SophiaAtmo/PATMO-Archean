module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !O2 -> O + O
    krate(:,228) = integrateXsec(1, tau(:,:))

    !O3 -> O2 + O(1D)
    krate(:,229) = integrateXsec(2, tau(:,:))

    !O3 -> O2 + O
    krate(:,230) = integrateXsec(3, tau(:,:))

    !OH -> O + H
    krate(:,231) = integrateXsec(4, tau(:,:))

    !OH -> O(1D) + H
    krate(:,232) = integrateXsec(5, tau(:,:))

    !HO2 -> OH + O
    krate(:,233) = integrateXsec(6, tau(:,:))

    !H2O -> OH + H
    krate(:,234) = integrateXsec(7, tau(:,:))

    !H2O -> H2 + O
    krate(:,235) = integrateXsec(8, tau(:,:))

    !H2 -> H + H
    krate(:,236) = integrateXsec(9, tau(:,:))

    !N2O -> N2 + O(1D)
    krate(:,237) = integrateXsec(10, tau(:,:))

    !NO2 -> NO + O
    krate(:,238) = integrateXsec(11, tau(:,:))

    !NO3 -> NO + O2
    krate(:,239) = integrateXsec(12, tau(:,:))

    !NO3 -> O + NO2
    krate(:,240) = integrateXsec(13, tau(:,:))

    !N2O5 -> NO2 + NO3
    krate(:,241) = integrateXsec(14, tau(:,:))

    !N2O5 -> O + NO + NO3
    krate(:,242) = integrateXsec(15, tau(:,:))

    !HNO3 -> OH + NO2
    krate(:,243) = integrateXsec(16, tau(:,:))

    !HNO3 -> H + NO3
    krate(:,244) = integrateXsec(17, tau(:,:))

    !CH4 -> CH3 + H
    krate(:,245) = integrateXsec(18, tau(:,:))

    !CH3OOH -> CH3O + OH
    krate(:,246) = integrateXsec(19, tau(:,:))

    !CH2O -> H + CHO
    krate(:,247) = integrateXsec(20, tau(:,:))

    !CH2O -> H2 + CO
    krate(:,248) = integrateXsec(21, tau(:,:))

    !CHO -> H + CO
    krate(:,249) = integrateXsec(22, tau(:,:))

    !CO2 -> CO + O
    krate(:,250) = integrateXsec(23, tau(:,:))

    !H2O2 -> OH + OH
    krate(:,251) = integrateXsec(24, tau(:,:))

    !H2O2 -> H + HO2
    krate(:,252) = integrateXsec(25, tau(:,:))

    !COS -> CO + S
    krate(:,253) = integrateXsec(26, tau(:,:))

    !SO -> S + O
    krate(:,254) = integrateXsec(27, tau(:,:))

    !CS2 -> CS + S
    krate(:,255) = integrateXsec(28, tau(:,:))

    !H2S -> SH + H
    krate(:,256) = integrateXsec(29, tau(:,:))

    !SO2 -> SO + O
    krate(:,257) = integrateXsec(30, tau(:,:))

    !SO3 -> SO2 + O
    krate(:,258) = integrateXsec(31, tau(:,:))

    !H2SO4 -> SO2 + OH + OH
    krate(:,259) = integrateXsec(32, tau(:,:))

    !CH3OH -> CH3 + OH
    krate(:,260) = integrateXsec(33, tau(:,:))

    !CH3OH -> CH3O + H
    krate(:,261) = integrateXsec(34, tau(:,:))

    !S2O2 -> SO + SO
    krate(:,262) = integrateXsec(35, tau(:,:))

    !S2O -> SO + S
    krate(:,263) = integrateXsec(36, tau(:,:))

    !N2H4 -> H + N2H3
    krate(:,264) = integrateXsec(37, tau(:,:))

    !NH3 -> H + NH2
    krate(:,265) = integrateXsec(38, tau(:,:))

    !NH3 -> H2 + NH
    krate(:,266) = integrateXsec(39, tau(:,:))

    !HCN -> CN + H
    krate(:,267) = integrateXsec(40, tau(:,:))

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber), dE, mu, coef
    integer::j

    ! !loop on cells (stride photobins)
    ! do j=1,cellsNumber
    !    integrateXsec(j) = sum(xsecAll(:,index)*photoFlux(:) &
        !         /energyMid(:)*energySpan(:)*exp(-tau(:,j))) / planck_eV
    ! end do

    !dE = (wavelengMax-wavelengMin)/photoBinsNumber (nm)
    dE = 0.1
    !mu =cosine(zenith_angle)
    mu = 0.500000
    !Parameter of incident solar flux = coef (Cronin, 2014- DOI:10.1175/JAS-D-13-0392.1)
    coef = 0.500000
    !loop on cells (stride photobins)
    do j=1,cellsNumber
      integrateXsec(j) = sum(xsecAll(:,index)*coef*photoFlux(:)*exp(-tau(:,j)/mu)*dE)
    end do

  end function integrateXsec

end module patmo_photoRates

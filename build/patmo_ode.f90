module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    use patmo_rates
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
      n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    call computeTotalDensityExcludingM(n)

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

    !forward grid points
    do j=1,cellsNumber-1
      dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j) = Tgas(j+1)
      ngas_p(j) = ngas(j+1)
      ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
      n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
      dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
      kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
      Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
      Tgas_m(j) = Tgas(j-1)
      ngas_m(j) = ngas(j-1)
      ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
      n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
      prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,i) =  dzzh_hp(:) &
          * (prem(:)*iTgas_hp(:) &
          - therm_hp(:))
      d_hm(:,i) = dzzh_hm(:) &
          * (prem(:)*iTgas_hm(:) &
          - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    dn(:,:) = 0d0
    dn(:,patmo_idx_O) = &
        - krate(:,1)*n(:,patmo_idx_O)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        - krate(:,2)*n(:,patmo_idx_O)*n(:,patmo_idx_O3) &
        + krate(:,4)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        + krate(:,4)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        + krate(:,5)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2) &
        + krate(:,6)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O2) &
        - krate(:,15)*n(:,patmo_idx_O)*n(:,patmo_idx_NO2) &
        - krate(:,23)*n(:,patmo_idx_O)*n(:,patmo_idx_OH) &
        - krate(:,25)*n(:,patmo_idx_O)*n(:,patmo_idx_HO2) &
        + krate(:,27)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,42)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,44)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,45)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,47)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,49)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,52)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,60)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,72)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,75)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,83)*n(:,patmo_idx_O)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        - krate(:,88)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        - krate(:,89)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        - krate(:,91)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        - krate(:,108)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,109)*n(:,patmo_idx_O)*n(:,patmo_idx_CH3) &
        - krate(:,110)*n(:,patmo_idx_O)*n(:,patmo_idx_CH3) &
        + krate(:,115)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        - krate(:,127)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        - krate(:,128)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        - krate(:,129)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        + krate(:,134)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,139)*n(:,patmo_idx_O)*n(:,patmo_idx_CH) &
        + krate(:,141)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        - krate(:,148)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        - krate(:,149)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        - krate(:,159)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_O) &
        + krate(:,169)*n(:,patmo_idx_N)*n(:,patmo_idx_O2) &
        + krate(:,170)*n(:,patmo_idx_N)*n(:,patmo_idx_NO) &
        - krate(:,172)*n(:,patmo_idx_O)*n(:,patmo_idx_NO3) &
        - krate(:,177)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        - krate(:,179)*n(:,patmo_idx_NH2)*n(:,patmo_idx_O) &
        - krate(:,184)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        - krate(:,190)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,191)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,195)*n(:,patmo_idx_O)*n(:,patmo_idx_N)*n(:,patmo_idx_M) &
        + krate(:,197)*n(:,patmo_idx_NO2)*n(:,patmo_idx_N) &
        - krate(:,198)*n(:,patmo_idx_O)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,198)*n(:,patmo_idx_O)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,212)*n(:,patmo_idx_O)*n(:,patmo_idx_H2) &
        + krate(:,215)*n(:,patmo_idx_O2) &
        + krate(:,215)*n(:,patmo_idx_O2) &
        + krate(:,217)*n(:,patmo_idx_O3) &
        + krate(:,218)*n(:,patmo_idx_OH) &
        + krate(:,220)*n(:,patmo_idx_HO2) &
        + krate(:,222)*n(:,patmo_idx_H2O) &
        + krate(:,225)*n(:,patmo_idx_NO2) &
        + krate(:,227)*n(:,patmo_idx_NO3) &
        + krate(:,229)*n(:,patmo_idx_N2O5) &
        + krate(:,237)*n(:,patmo_idx_CO2) &
        + krate(:,241)*n(:,patmo_idx_SO) &
        + krate(:,244)*n(:,patmo_idx_SO2) &
        + krate(:,245)*n(:,patmo_idx_SO3) &
        + krate(:,254)*n(:,patmo_idx_O3)*n(:,patmo_idx_M) &
        + krate(:,255)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        - krate(:,257)*n(:,patmo_idx_O2)*n(:,patmo_idx_O)*n(:,patmo_idx_O) &
        - krate(:,257)*n(:,patmo_idx_O2)*n(:,patmo_idx_O)*n(:,patmo_idx_O) &
        - krate(:,258)*n(:,patmo_idx_O)*n(:,patmo_idx_N2) &
        - krate(:,259)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        + krate(:,268)*n(:,patmo_idx_NO)*n(:,patmo_idx_O2) &
        + krate(:,276)*n(:,patmo_idx_H)*n(:,patmo_idx_O2) &
        + krate(:,278)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        - krate(:,280)*n(:,patmo_idx_O)*n(:,patmo_idx_H2O) &
        + krate(:,295)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        + krate(:,297)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,298)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,300)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,302)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,305)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,310)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,313)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,325)*n(:,patmo_idx_SO3)*n(:,patmo_idx_M) &
        + krate(:,328)*n(:,patmo_idx_SO2) &
        + krate(:,336)*n(:,patmo_idx_CO2)*n(:,patmo_idx_M) &
        + krate(:,341)*n(:,patmo_idx_H)*n(:,patmo_idx_CO2) &
        + krate(:,342)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        + krate(:,344)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        + krate(:,361)*n(:,patmo_idx_S)*n(:,patmo_idx_SO) &
        + krate(:,362)*n(:,patmo_idx_CH3O) &
        + krate(:,363)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H) &
        - krate(:,368)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O) &
        + krate(:,380)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H) &
        + krate(:,381)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        + krate(:,382)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        - krate(:,387)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        + krate(:,392)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        - krate(:,394)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        + krate(:,401)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        + krate(:,402)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        + krate(:,412)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        - krate(:,422)*n(:,patmo_idx_O)*n(:,patmo_idx_NO) &
        - krate(:,423)*n(:,patmo_idx_N2)*n(:,patmo_idx_O) &
        + krate(:,425)*n(:,patmo_idx_O2)*n(:,patmo_idx_NO2) &
        + krate(:,430)*n(:,patmo_idx_N)*n(:,patmo_idx_OH) &
        + krate(:,432)*n(:,patmo_idx_NH)*n(:,patmo_idx_OH) &
        + krate(:,437)*n(:,patmo_idx_NO)*n(:,patmo_idx_H) &
        + krate(:,443)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        + krate(:,444)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        + krate(:,448)*n(:,patmo_idx_NO)*n(:,patmo_idx_M) &
        - krate(:,450)*n(:,patmo_idx_N2O)*n(:,patmo_idx_O) &
        + krate(:,451)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,451)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,465)*n(:,patmo_idx_OH)*n(:,patmo_idx_H)

    dn(:,patmo_idx_O2) = &
        - krate(:,1)*n(:,patmo_idx_O)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,2)*n(:,patmo_idx_O)*n(:,patmo_idx_O3) &
        + krate(:,2)*n(:,patmo_idx_O)*n(:,patmo_idx_O3) &
        + krate(:,3)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        + krate(:,3)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        + krate(:,4)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        - krate(:,6)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O2) &
        + krate(:,6)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O2) &
        + krate(:,7)*n(:,patmo_idx_OH)*n(:,patmo_idx_O3) &
        + krate(:,8)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O3) &
        + krate(:,8)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O3) &
        + krate(:,9)*n(:,patmo_idx_OH)*n(:,patmo_idx_HO2) &
        + krate(:,13)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2O) &
        + krate(:,15)*n(:,patmo_idx_O)*n(:,patmo_idx_NO2) &
        + krate(:,16)*n(:,patmo_idx_NO)*n(:,patmo_idx_O3) &
        + krate(:,17)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O3) &
        + krate(:,22)*n(:,patmo_idx_H)*n(:,patmo_idx_O3) &
        + krate(:,23)*n(:,patmo_idx_O)*n(:,patmo_idx_OH) &
        - krate(:,24)*n(:,patmo_idx_H)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,25)*n(:,patmo_idx_O)*n(:,patmo_idx_HO2) &
        + krate(:,28)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,30)*n(:,patmo_idx_CH3)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,31)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_HO2) &
        - krate(:,35)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        - krate(:,37)*n(:,patmo_idx_CHO)*n(:,patmo_idx_O2) &
        + krate(:,39)*n(:,patmo_idx_HO2)*n(:,patmo_idx_HO2)*n(:,patmo_idx_M) &
        - krate(:,45)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,46)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,53)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        + krate(:,54)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,60)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,61)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,65)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        - krate(:,66)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,67)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,67)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,69)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,70)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,81)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        + krate(:,82)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        + krate(:,112)*n(:,patmo_idx_O3)*n(:,patmo_idx_CH3) &
        + krate(:,118)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        - krate(:,131)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,132)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,133)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,134)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,141)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        - krate(:,142)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        - krate(:,169)*n(:,patmo_idx_N)*n(:,patmo_idx_O2) &
        + krate(:,172)*n(:,patmo_idx_O)*n(:,patmo_idx_NO3) &
        + krate(:,198)*n(:,patmo_idx_O)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,214)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_O2) &
        - krate(:,215)*n(:,patmo_idx_O2) &
        + krate(:,216)*n(:,patmo_idx_O3) &
        + krate(:,217)*n(:,patmo_idx_O3) &
        + krate(:,226)*n(:,patmo_idx_NO3) &
        + krate(:,254)*n(:,patmo_idx_O3)*n(:,patmo_idx_M) &
        - krate(:,255)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        - krate(:,255)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        - krate(:,256)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        - krate(:,256)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        - krate(:,257)*n(:,patmo_idx_O2)*n(:,patmo_idx_O)*n(:,patmo_idx_O) &
        - krate(:,259)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        + krate(:,259)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        - krate(:,260)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O2) &
        - krate(:,261)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        - krate(:,261)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        - krate(:,262)*n(:,patmo_idx_H2O)*n(:,patmo_idx_O2) &
        - krate(:,266)*n(:,patmo_idx_N2)*n(:,patmo_idx_O2) &
        - krate(:,268)*n(:,patmo_idx_NO)*n(:,patmo_idx_O2) &
        - krate(:,269)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O2) &
        - krate(:,270)*n(:,patmo_idx_NO3)*n(:,patmo_idx_O2) &
        - krate(:,275)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        - krate(:,276)*n(:,patmo_idx_H)*n(:,patmo_idx_O2) &
        + krate(:,277)*n(:,patmo_idx_HO2)*n(:,patmo_idx_M) &
        - krate(:,278)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        - krate(:,281)*n(:,patmo_idx_H2)*n(:,patmo_idx_O2) &
        + krate(:,283)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_M) &
        - krate(:,284)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_O2) &
        + krate(:,288)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_HO2) &
        + krate(:,290)*n(:,patmo_idx_CO)*n(:,patmo_idx_HO2) &
        - krate(:,292)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,298)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,299)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,306)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        - krate(:,307)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,309)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,310)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,313)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,314)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        - krate(:,318)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        + krate(:,319)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,320)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        - krate(:,320)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        + krate(:,322)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        + krate(:,323)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,334)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        - krate(:,335)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_O2) &
        - krate(:,365)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        - krate(:,371)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O2) &
        + krate(:,384)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO2) &
        + krate(:,385)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO2) &
        + krate(:,386)*n(:,patmo_idx_CO)*n(:,patmo_idx_H2O) &
        + krate(:,387)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        + krate(:,394)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        + krate(:,395)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        + krate(:,422)*n(:,patmo_idx_O)*n(:,patmo_idx_NO) &
        - krate(:,425)*n(:,patmo_idx_O2)*n(:,patmo_idx_NO2) &
        - krate(:,451)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,467)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CO2)

    dn(:,patmo_idx_O3) = &
        + krate(:,1)*n(:,patmo_idx_O)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        - krate(:,2)*n(:,patmo_idx_O)*n(:,patmo_idx_O3) &
        - krate(:,3)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        - krate(:,4)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        - krate(:,7)*n(:,patmo_idx_OH)*n(:,patmo_idx_O3) &
        - krate(:,8)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O3) &
        - krate(:,16)*n(:,patmo_idx_NO)*n(:,patmo_idx_O3) &
        - krate(:,17)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O3) &
        - krate(:,22)*n(:,patmo_idx_H)*n(:,patmo_idx_O3) &
        - krate(:,46)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,54)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,61)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,65)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        - krate(:,67)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,112)*n(:,patmo_idx_O3)*n(:,patmo_idx_CH3) &
        - krate(:,216)*n(:,patmo_idx_O3) &
        - krate(:,217)*n(:,patmo_idx_O3) &
        - krate(:,254)*n(:,patmo_idx_O3)*n(:,patmo_idx_M) &
        + krate(:,255)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        + krate(:,256)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        + krate(:,257)*n(:,patmo_idx_O2)*n(:,patmo_idx_O)*n(:,patmo_idx_O) &
        + krate(:,260)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O2) &
        + krate(:,261)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        + krate(:,269)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O2) &
        + krate(:,270)*n(:,patmo_idx_NO3)*n(:,patmo_idx_O2) &
        + krate(:,275)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        + krate(:,299)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,307)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,309)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,314)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        + krate(:,318)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        + krate(:,320)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        + krate(:,365)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_O_1D) = &
        - krate(:,3)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        - krate(:,4)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O3) &
        - krate(:,5)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2) &
        - krate(:,6)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_O2) &
        - krate(:,10)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2O) &
        - krate(:,12)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2) &
        - krate(:,13)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2O) &
        - krate(:,14)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2O) &
        - krate(:,78)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,79)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,80)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,92)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2) &
        - krate(:,180)*n(:,patmo_idx_NH3)*n(:,patmo_idx_O_1D) &
        - krate(:,208)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CO2) &
        - krate(:,209)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2) &
        - krate(:,210)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_SO2) &
        + krate(:,216)*n(:,patmo_idx_O3) &
        + krate(:,219)*n(:,patmo_idx_OH) &
        + krate(:,224)*n(:,patmo_idx_N2O) &
        + krate(:,256)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        + krate(:,257)*n(:,patmo_idx_O2)*n(:,patmo_idx_O)*n(:,patmo_idx_O) &
        + krate(:,258)*n(:,patmo_idx_O)*n(:,patmo_idx_N2) &
        + krate(:,259)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        + krate(:,263)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,265)*n(:,patmo_idx_N2O) &
        + krate(:,266)*n(:,patmo_idx_N2)*n(:,patmo_idx_O2) &
        + krate(:,267)*n(:,patmo_idx_NO)*n(:,patmo_idx_NO) &
        + krate(:,331)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,332)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        + krate(:,333)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        + krate(:,345)*n(:,patmo_idx_H)*n(:,patmo_idx_OH) &
        + krate(:,433)*n(:,patmo_idx_NH2)*n(:,patmo_idx_OH) &
        + krate(:,461)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_CO2) &
        + krate(:,462)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_N2) &
        + krate(:,463)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_OH) = &
        - krate(:,7)*n(:,patmo_idx_OH)*n(:,patmo_idx_O3) &
        + krate(:,8)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O3) &
        - krate(:,9)*n(:,patmo_idx_OH)*n(:,patmo_idx_HO2) &
        + krate(:,10)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2O) &
        + krate(:,10)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2O) &
        + krate(:,11)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        - krate(:,19)*n(:,patmo_idx_NO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,20)*n(:,patmo_idx_HNO3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,21)*n(:,patmo_idx_HO2)*n(:,patmo_idx_NO) &
        + krate(:,22)*n(:,patmo_idx_H)*n(:,patmo_idx_O3) &
        - krate(:,23)*n(:,patmo_idx_O)*n(:,patmo_idx_OH) &
        + krate(:,25)*n(:,patmo_idx_O)*n(:,patmo_idx_HO2) &
        + krate(:,26)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        + krate(:,26)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,29)*n(:,patmo_idx_CH4)*n(:,patmo_idx_OH) &
        - krate(:,33)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        + krate(:,33)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        - krate(:,34)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        - krate(:,36)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        - krate(:,38)*n(:,patmo_idx_CO)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,40)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_OH) &
        - krate(:,41)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,43)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,48)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,49)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,53)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,58)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,62)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,66)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,73)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,74)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,76)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,78)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,87)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        + krate(:,89)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        + krate(:,91)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        + krate(:,92)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2) &
        - krate(:,93)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        + krate(:,94)*n(:,patmo_idx_SO)*n(:,patmo_idx_HO2) &
        - krate(:,114)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        - krate(:,115)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        - krate(:,116)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3)*n(:,patmo_idx_M) &
        + krate(:,117)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        - krate(:,135)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH2) &
        + krate(:,142)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        + krate(:,148)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        + krate(:,149)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        - krate(:,153)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        - krate(:,154)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        - krate(:,155)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        + krate(:,159)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_O) &
        + krate(:,160)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,164)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        + krate(:,171)*n(:,patmo_idx_H)*n(:,patmo_idx_NO2) &
        + krate(:,176)*n(:,patmo_idx_NH)*n(:,patmo_idx_NO) &
        + krate(:,177)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        + krate(:,179)*n(:,patmo_idx_NH2)*n(:,patmo_idx_O) &
        + krate(:,180)*n(:,patmo_idx_NH3)*n(:,patmo_idx_O_1D) &
        - krate(:,181)*n(:,patmo_idx_NH3)*n(:,patmo_idx_OH) &
        - krate(:,192)*n(:,patmo_idx_OH)*n(:,patmo_idx_NH2) &
        - krate(:,199)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        + krate(:,200)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_O_3P) &
        - krate(:,201)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_OH) &
        - krate(:,207)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,207)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,212)*n(:,patmo_idx_O)*n(:,patmo_idx_H2) &
        - krate(:,218)*n(:,patmo_idx_OH) &
        - krate(:,219)*n(:,patmo_idx_OH) &
        + krate(:,220)*n(:,patmo_idx_HO2) &
        + krate(:,221)*n(:,patmo_idx_H2O) &
        + krate(:,230)*n(:,patmo_idx_HNO3) &
        + krate(:,233)*n(:,patmo_idx_CH3OOH) &
        + krate(:,238)*n(:,patmo_idx_H2O2) &
        + krate(:,238)*n(:,patmo_idx_H2O2) &
        + krate(:,246)*n(:,patmo_idx_H2SO4) &
        + krate(:,246)*n(:,patmo_idx_H2SO4) &
        + krate(:,247)*n(:,patmo_idx_CH3OH) &
        + krate(:,260)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O2) &
        - krate(:,261)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        + krate(:,262)*n(:,patmo_idx_H2O)*n(:,patmo_idx_O2) &
        - krate(:,263)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        - krate(:,263)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        - krate(:,264)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        + krate(:,272)*n(:,patmo_idx_HNO3)*n(:,patmo_idx_M) &
        + krate(:,273)*n(:,patmo_idx_NO3)*n(:,patmo_idx_H2O)*n(:,patmo_idx_M) &
        - krate(:,274)*n(:,patmo_idx_OH)*n(:,patmo_idx_NO2) &
        - krate(:,275)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        + krate(:,276)*n(:,patmo_idx_H)*n(:,patmo_idx_O2) &
        - krate(:,278)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        - krate(:,279)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        - krate(:,279)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,282)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        - krate(:,286)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2O) &
        + krate(:,286)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2O) &
        + krate(:,287)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_H2O) &
        + krate(:,289)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H2O) &
        + krate(:,291)*n(:,patmo_idx_CO2)*n(:,patmo_idx_H)*n(:,patmo_idx_M) &
        + krate(:,293)*n(:,patmo_idx_HO2)*n(:,patmo_idx_H2O) &
        + krate(:,294)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,296)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,301)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,302)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,306)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,311)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,315)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,316)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,319)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,326)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_M) &
        + krate(:,327)*n(:,patmo_idx_SO2) &
        + krate(:,329)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)*n(:,patmo_idx_M) &
        - krate(:,331)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,340)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO) &
        - krate(:,342)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        - krate(:,344)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        - krate(:,345)*n(:,patmo_idx_H)*n(:,patmo_idx_OH) &
        + krate(:,346)*n(:,patmo_idx_H)*n(:,patmo_idx_H2O) &
        - krate(:,347)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,367)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        + krate(:,368)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O) &
        + krate(:,369)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_M) &
        - krate(:,370)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        + krate(:,388)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,395)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        - krate(:,401)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        - krate(:,402)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        + krate(:,406)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2O) &
        + krate(:,407)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O) &
        + krate(:,408)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        - krate(:,412)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        - krate(:,413)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,417)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O) &
        - krate(:,424)*n(:,patmo_idx_NO)*n(:,patmo_idx_OH) &
        - krate(:,429)*n(:,patmo_idx_N2)*n(:,patmo_idx_OH) &
        - krate(:,430)*n(:,patmo_idx_N)*n(:,patmo_idx_OH) &
        - krate(:,432)*n(:,patmo_idx_NH)*n(:,patmo_idx_OH) &
        - krate(:,433)*n(:,patmo_idx_NH2)*n(:,patmo_idx_OH) &
        + krate(:,434)*n(:,patmo_idx_NH2)*n(:,patmo_idx_H2O) &
        + krate(:,445)*n(:,patmo_idx_H2O)*n(:,patmo_idx_NH) &
        + krate(:,452)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_M) &
        - krate(:,453)*n(:,patmo_idx_CO2)*n(:,patmo_idx_OH) &
        + krate(:,454)*n(:,patmo_idx_CO2)*n(:,patmo_idx_H2O) &
        + krate(:,460)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_M) &
        + krate(:,460)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_M) &
        - krate(:,465)*n(:,patmo_idx_OH)*n(:,patmo_idx_H)

    dn(:,patmo_idx_HO2) = &
        + krate(:,7)*n(:,patmo_idx_OH)*n(:,patmo_idx_O3) &
        - krate(:,8)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O3) &
        - krate(:,9)*n(:,patmo_idx_OH)*n(:,patmo_idx_HO2) &
        - krate(:,21)*n(:,patmo_idx_HO2)*n(:,patmo_idx_NO) &
        + krate(:,24)*n(:,patmo_idx_H)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        - krate(:,25)*n(:,patmo_idx_O)*n(:,patmo_idx_HO2) &
        - krate(:,26)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,27)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,28)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,31)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_HO2) &
        + krate(:,35)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        + krate(:,37)*n(:,patmo_idx_CHO)*n(:,patmo_idx_O2) &
        - krate(:,39)*n(:,patmo_idx_HO2)*n(:,patmo_idx_HO2)*n(:,patmo_idx_M) &
        - krate(:,39)*n(:,patmo_idx_HO2)*n(:,patmo_idx_HO2)*n(:,patmo_idx_M) &
        + krate(:,40)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_OH) &
        - krate(:,51)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,69)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        + krate(:,70)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        - krate(:,94)*n(:,patmo_idx_SO)*n(:,patmo_idx_HO2) &
        + krate(:,113)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_CH3) &
        - krate(:,117)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        - krate(:,118)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        + krate(:,163)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O2) &
        - krate(:,165)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_HO2) &
        + krate(:,214)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_O2) &
        - krate(:,220)*n(:,patmo_idx_HO2) &
        + krate(:,239)*n(:,patmo_idx_H2O2) &
        - krate(:,260)*n(:,patmo_idx_HO2)*n(:,patmo_idx_O2) &
        + krate(:,261)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2) &
        + krate(:,262)*n(:,patmo_idx_H2O)*n(:,patmo_idx_O2) &
        + krate(:,274)*n(:,patmo_idx_OH)*n(:,patmo_idx_NO2) &
        - krate(:,277)*n(:,patmo_idx_HO2)*n(:,patmo_idx_M) &
        + krate(:,278)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        + krate(:,279)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,280)*n(:,patmo_idx_O)*n(:,patmo_idx_H2O) &
        + krate(:,281)*n(:,patmo_idx_H2)*n(:,patmo_idx_O2) &
        + krate(:,284)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_O2) &
        - krate(:,288)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_HO2) &
        - krate(:,290)*n(:,patmo_idx_CO)*n(:,patmo_idx_HO2) &
        + krate(:,292)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,292)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        - krate(:,293)*n(:,patmo_idx_HO2)*n(:,patmo_idx_H2O) &
        + krate(:,304)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        + krate(:,316)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,322)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        - krate(:,323)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        + krate(:,347)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,366)*n(:,patmo_idx_CH4)*n(:,patmo_idx_HO2) &
        + krate(:,370)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        + krate(:,371)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O2) &
        - krate(:,416)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_HO2) &
        + krate(:,418)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O2) &
        - krate(:,467)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CO2)

    dn(:,patmo_idx_H2O) = &
        + krate(:,9)*n(:,patmo_idx_OH)*n(:,patmo_idx_HO2) &
        - krate(:,10)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2O) &
        - krate(:,11)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        + krate(:,20)*n(:,patmo_idx_HNO3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,27)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        + krate(:,29)*n(:,patmo_idx_CH4)*n(:,patmo_idx_OH) &
        + krate(:,33)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        + krate(:,34)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        + krate(:,36)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        + krate(:,40)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_OH) &
        + krate(:,48)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,51)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,71)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        + krate(:,87)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        + krate(:,93)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        + krate(:,133)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,143)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH) &
        + krate(:,150)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        + krate(:,153)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        + krate(:,154)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        + krate(:,155)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        + krate(:,164)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        + krate(:,178)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NO) &
        + krate(:,181)*n(:,patmo_idx_NH3)*n(:,patmo_idx_OH) &
        + krate(:,192)*n(:,patmo_idx_OH)*n(:,patmo_idx_NH2) &
        + krate(:,201)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_OH) &
        + krate(:,202)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CH3) &
        + krate(:,204)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_H) &
        - krate(:,221)*n(:,patmo_idx_H2O) &
        - krate(:,222)*n(:,patmo_idx_H2O) &
        - krate(:,262)*n(:,patmo_idx_H2O)*n(:,patmo_idx_O2) &
        + krate(:,263)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,264)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        - krate(:,273)*n(:,patmo_idx_NO3)*n(:,patmo_idx_H2O)*n(:,patmo_idx_M) &
        - krate(:,280)*n(:,patmo_idx_O)*n(:,patmo_idx_H2O) &
        - krate(:,282)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        - krate(:,286)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2O) &
        - krate(:,287)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_H2O) &
        - krate(:,289)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H2O) &
        - krate(:,293)*n(:,patmo_idx_HO2)*n(:,patmo_idx_H2O) &
        - krate(:,301)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,304)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        + krate(:,324)*n(:,patmo_idx_H2SO4) &
        - krate(:,340)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO) &
        - krate(:,346)*n(:,patmo_idx_H)*n(:,patmo_idx_H2O) &
        - krate(:,386)*n(:,patmo_idx_CO)*n(:,patmo_idx_H2O) &
        + krate(:,396)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,403)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        - krate(:,406)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2O) &
        - krate(:,407)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O) &
        - krate(:,408)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        - krate(:,417)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O) &
        - krate(:,431)*n(:,patmo_idx_N2)*n(:,patmo_idx_H2O) &
        - krate(:,434)*n(:,patmo_idx_NH2)*n(:,patmo_idx_H2O) &
        - krate(:,445)*n(:,patmo_idx_H2O)*n(:,patmo_idx_NH) &
        - krate(:,454)*n(:,patmo_idx_CO2)*n(:,patmo_idx_H2O) &
        - krate(:,455)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH2CO) &
        - krate(:,457)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO)

    dn(:,patmo_idx_H) = &
        - krate(:,11)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        - krate(:,22)*n(:,patmo_idx_H)*n(:,patmo_idx_O3) &
        + krate(:,23)*n(:,patmo_idx_O)*n(:,patmo_idx_OH) &
        - krate(:,24)*n(:,patmo_idx_H)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        - krate(:,26)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,27)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        - krate(:,28)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        + krate(:,38)*n(:,patmo_idx_CO)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,50)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        + krate(:,52)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,58)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,62)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,79)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,84)*n(:,patmo_idx_H)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        - krate(:,85)*n(:,patmo_idx_H)*n(:,patmo_idx_CHO) &
        + krate(:,88)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        - krate(:,90)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        + krate(:,92)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2) &
        + krate(:,93)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        - krate(:,100)*n(:,patmo_idx_SH)*n(:,patmo_idx_H) &
        + krate(:,110)*n(:,patmo_idx_O)*n(:,patmo_idx_CH3) &
        - krate(:,111)*n(:,patmo_idx_H)*n(:,patmo_idx_CH3)*n(:,patmo_idx_M) &
        + krate(:,114)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        + krate(:,121)*n(:,patmo_idx_CH3) &
        + krate(:,126)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH3) &
        + krate(:,127)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        + krate(:,128)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        + krate(:,128)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        - krate(:,130)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,131)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        + krate(:,131)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        + krate(:,135)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH2) &
        + krate(:,139)*n(:,patmo_idx_O)*n(:,patmo_idx_CH) &
        + krate(:,143)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH) &
        + krate(:,144)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,150)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        - krate(:,151)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        - krate(:,152)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        + krate(:,155)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        - krate(:,160)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,161)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,162)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,171)*n(:,patmo_idx_H)*n(:,patmo_idx_NO2) &
        - krate(:,174)*n(:,patmo_idx_N2H4)*n(:,patmo_idx_H) &
        - krate(:,175)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H) &
        - krate(:,182)*n(:,patmo_idx_NH2)*n(:,patmo_idx_H)*n(:,patmo_idx_M) &
        + krate(:,183)*n(:,patmo_idx_NH)*n(:,patmo_idx_NO) &
        + krate(:,184)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        - krate(:,186)*n(:,patmo_idx_COS)*n(:,patmo_idx_H) &
        + krate(:,189)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        - krate(:,196)*n(:,patmo_idx_H)*n(:,patmo_idx_N)*n(:,patmo_idx_M) &
        - krate(:,204)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_H) &
        - krate(:,205)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_H) &
        + krate(:,212)*n(:,patmo_idx_O)*n(:,patmo_idx_H2) &
        - krate(:,213)*n(:,patmo_idx_H)*n(:,patmo_idx_H) &
        - krate(:,213)*n(:,patmo_idx_H)*n(:,patmo_idx_H) &
        + krate(:,218)*n(:,patmo_idx_OH) &
        + krate(:,219)*n(:,patmo_idx_OH) &
        + krate(:,221)*n(:,patmo_idx_H2O) &
        + krate(:,223)*n(:,patmo_idx_H2) &
        + krate(:,223)*n(:,patmo_idx_H2) &
        + krate(:,231)*n(:,patmo_idx_HNO3) &
        + krate(:,232)*n(:,patmo_idx_CH4) &
        + krate(:,234)*n(:,patmo_idx_CH2O) &
        + krate(:,236)*n(:,patmo_idx_CHO) &
        + krate(:,239)*n(:,patmo_idx_H2O2) &
        + krate(:,243)*n(:,patmo_idx_H2S) &
        + krate(:,248)*n(:,patmo_idx_CH3OH) &
        + krate(:,251)*n(:,patmo_idx_N2H4) &
        + krate(:,252)*n(:,patmo_idx_NH3) &
        + krate(:,264)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        + krate(:,275)*n(:,patmo_idx_OH)*n(:,patmo_idx_O2) &
        - krate(:,276)*n(:,patmo_idx_H)*n(:,patmo_idx_O2) &
        + krate(:,277)*n(:,patmo_idx_HO2)*n(:,patmo_idx_M) &
        + krate(:,279)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,280)*n(:,patmo_idx_O)*n(:,patmo_idx_H2O) &
        + krate(:,281)*n(:,patmo_idx_H2)*n(:,patmo_idx_O2) &
        - krate(:,291)*n(:,patmo_idx_CO2)*n(:,patmo_idx_H)*n(:,patmo_idx_M) &
        + krate(:,303)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        - krate(:,305)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,311)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,315)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,332)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        + krate(:,337)*n(:,patmo_idx_CHO)*n(:,patmo_idx_M) &
        + krate(:,338)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        - krate(:,341)*n(:,patmo_idx_H)*n(:,patmo_idx_CO2) &
        + krate(:,343)*n(:,patmo_idx_H2)*n(:,patmo_idx_CHO) &
        - krate(:,345)*n(:,patmo_idx_H)*n(:,patmo_idx_OH) &
        - krate(:,346)*n(:,patmo_idx_H)*n(:,patmo_idx_H2O) &
        + krate(:,353)*n(:,patmo_idx_H2)*n(:,patmo_idx_S) &
        - krate(:,363)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H) &
        + krate(:,364)*n(:,patmo_idx_CH4)*n(:,patmo_idx_M) &
        - krate(:,367)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        - krate(:,374)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        - krate(:,379)*n(:,patmo_idx_CH4)*n(:,patmo_idx_H) &
        - krate(:,380)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H) &
        - krate(:,381)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        - krate(:,381)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        + krate(:,383)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,384)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO2) &
        - krate(:,384)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO2) &
        - krate(:,388)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,392)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        - krate(:,396)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,397)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,403)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        + krate(:,404)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2) &
        + krate(:,405)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2) &
        - krate(:,408)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        + krate(:,413)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,414)*n(:,patmo_idx_CH3OH) &
        + krate(:,415)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        + krate(:,424)*n(:,patmo_idx_NO)*n(:,patmo_idx_OH) &
        + krate(:,427)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H2) &
        + krate(:,428)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH2) &
        + krate(:,435)*n(:,patmo_idx_NH3)*n(:,patmo_idx_M) &
        - krate(:,436)*n(:,patmo_idx_N2O)*n(:,patmo_idx_H) &
        - krate(:,437)*n(:,patmo_idx_NO)*n(:,patmo_idx_H) &
        + krate(:,439)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        - krate(:,442)*n(:,patmo_idx_COS)*n(:,patmo_idx_H) &
        + krate(:,449)*n(:,patmo_idx_NH)*n(:,patmo_idx_M) &
        + krate(:,457)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO) &
        + krate(:,458)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO2) &
        - krate(:,465)*n(:,patmo_idx_OH)*n(:,patmo_idx_H) &
        + krate(:,466)*n(:,patmo_idx_H2) &
        + krate(:,466)*n(:,patmo_idx_H2)

    dn(:,patmo_idx_H2) = &
        + krate(:,11)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        + krate(:,28)*n(:,patmo_idx_H)*n(:,patmo_idx_HO2) &
        + krate(:,50)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        + krate(:,80)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        + krate(:,85)*n(:,patmo_idx_H)*n(:,patmo_idx_CHO) &
        + krate(:,90)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,92)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_H2) &
        - krate(:,93)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        + krate(:,100)*n(:,patmo_idx_SH)*n(:,patmo_idx_H) &
        + krate(:,122)*n(:,patmo_idx_CH3) &
        - krate(:,126)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH3) &
        + krate(:,129)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        + krate(:,130)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,132)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,144)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,145)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        + krate(:,151)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        + krate(:,152)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        + krate(:,162)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        + krate(:,174)*n(:,patmo_idx_N2H4)*n(:,patmo_idx_H) &
        + krate(:,205)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_H) &
        - krate(:,212)*n(:,patmo_idx_O)*n(:,patmo_idx_H2) &
        + krate(:,213)*n(:,patmo_idx_H)*n(:,patmo_idx_H) &
        + krate(:,222)*n(:,patmo_idx_H2O) &
        - krate(:,223)*n(:,patmo_idx_H2) &
        + krate(:,235)*n(:,patmo_idx_CH2O) &
        + krate(:,253)*n(:,patmo_idx_NH3) &
        - krate(:,264)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2) &
        - krate(:,281)*n(:,patmo_idx_H2)*n(:,patmo_idx_O2) &
        - krate(:,303)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        - krate(:,333)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        - krate(:,338)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        - krate(:,343)*n(:,patmo_idx_H2)*n(:,patmo_idx_CHO) &
        + krate(:,345)*n(:,patmo_idx_H)*n(:,patmo_idx_OH) &
        + krate(:,346)*n(:,patmo_idx_H)*n(:,patmo_idx_H2O) &
        - krate(:,353)*n(:,patmo_idx_H2)*n(:,patmo_idx_S) &
        - krate(:,375)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        + krate(:,379)*n(:,patmo_idx_CH4)*n(:,patmo_idx_H) &
        - krate(:,382)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        - krate(:,383)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,385)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO2) &
        + krate(:,397)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,398)*n(:,patmo_idx_CH3) &
        - krate(:,404)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2) &
        - krate(:,405)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2) &
        - krate(:,415)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        - krate(:,427)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H2) &
        - krate(:,458)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO2) &
        + krate(:,465)*n(:,patmo_idx_OH)*n(:,patmo_idx_H) &
        - krate(:,466)*n(:,patmo_idx_H2)

    dn(:,patmo_idx_N2O) = &
        + krate(:,12)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2) &
        - krate(:,13)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2O) &
        - krate(:,14)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2O) &
        + krate(:,183)*n(:,patmo_idx_NH)*n(:,patmo_idx_NO) &
        + krate(:,197)*n(:,patmo_idx_NO2)*n(:,patmo_idx_N) &
        - krate(:,224)*n(:,patmo_idx_N2O) &
        - krate(:,265)*n(:,patmo_idx_N2O) &
        + krate(:,266)*n(:,patmo_idx_N2)*n(:,patmo_idx_O2) &
        + krate(:,267)*n(:,patmo_idx_NO)*n(:,patmo_idx_NO) &
        - krate(:,436)*n(:,patmo_idx_N2O)*n(:,patmo_idx_H) &
        - krate(:,450)*n(:,patmo_idx_N2O)*n(:,patmo_idx_O)

    dn(:,patmo_idx_NO) = &
        + krate(:,14)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2O) &
        + krate(:,14)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2O) &
        + krate(:,15)*n(:,patmo_idx_O)*n(:,patmo_idx_NO2) &
        - krate(:,16)*n(:,patmo_idx_NO)*n(:,patmo_idx_O3) &
        - krate(:,21)*n(:,patmo_idx_HO2)*n(:,patmo_idx_NO) &
        - krate(:,32)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_NO) &
        + krate(:,55)*n(:,patmo_idx_SH)*n(:,patmo_idx_NO2) &
        + krate(:,59)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        + krate(:,64)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO2) &
        + krate(:,68)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO2) &
        + krate(:,140)*n(:,patmo_idx_CH)*n(:,patmo_idx_NO2) &
        + krate(:,169)*n(:,patmo_idx_N)*n(:,patmo_idx_O2) &
        - krate(:,170)*n(:,patmo_idx_N)*n(:,patmo_idx_NO) &
        + krate(:,171)*n(:,patmo_idx_H)*n(:,patmo_idx_NO2) &
        - krate(:,176)*n(:,patmo_idx_NH)*n(:,patmo_idx_NO) &
        - krate(:,178)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NO) &
        - krate(:,183)*n(:,patmo_idx_NH)*n(:,patmo_idx_NO) &
        + krate(:,184)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        + krate(:,188)*n(:,patmo_idx_CS)*n(:,patmo_idx_NO2) &
        + krate(:,195)*n(:,patmo_idx_O)*n(:,patmo_idx_N)*n(:,patmo_idx_M) &
        + krate(:,225)*n(:,patmo_idx_NO2) &
        + krate(:,226)*n(:,patmo_idx_NO3) &
        + krate(:,229)*n(:,patmo_idx_N2O5) &
        - krate(:,267)*n(:,patmo_idx_NO)*n(:,patmo_idx_NO) &
        - krate(:,267)*n(:,patmo_idx_NO)*n(:,patmo_idx_NO) &
        - krate(:,268)*n(:,patmo_idx_NO)*n(:,patmo_idx_O2) &
        + krate(:,269)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O2) &
        + krate(:,274)*n(:,patmo_idx_OH)*n(:,patmo_idx_NO2) &
        + krate(:,285)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_NO2) &
        - krate(:,308)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO) &
        - krate(:,312)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        - krate(:,317)*n(:,patmo_idx_SO3)*n(:,patmo_idx_NO) &
        - krate(:,321)*n(:,patmo_idx_NO)*n(:,patmo_idx_HSO2) &
        - krate(:,393)*n(:,patmo_idx_CHO)*n(:,patmo_idx_NO) &
        - krate(:,422)*n(:,patmo_idx_O)*n(:,patmo_idx_NO) &
        + krate(:,423)*n(:,patmo_idx_N2)*n(:,patmo_idx_O) &
        - krate(:,424)*n(:,patmo_idx_NO)*n(:,patmo_idx_OH) &
        + krate(:,429)*n(:,patmo_idx_N2)*n(:,patmo_idx_OH) &
        + krate(:,431)*n(:,patmo_idx_N2)*n(:,patmo_idx_H2O) &
        + krate(:,436)*n(:,patmo_idx_N2O)*n(:,patmo_idx_H) &
        - krate(:,437)*n(:,patmo_idx_NO)*n(:,patmo_idx_H) &
        - krate(:,441)*n(:,patmo_idx_COS)*n(:,patmo_idx_NO) &
        - krate(:,448)*n(:,patmo_idx_NO)*n(:,patmo_idx_M)

    dn(:,patmo_idx_NO2) = &
        - krate(:,15)*n(:,patmo_idx_O)*n(:,patmo_idx_NO2) &
        + krate(:,16)*n(:,patmo_idx_NO)*n(:,patmo_idx_O3) &
        - krate(:,17)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O3) &
        - krate(:,18)*n(:,patmo_idx_NO2)*n(:,patmo_idx_NO3)*n(:,patmo_idx_M) &
        - krate(:,19)*n(:,patmo_idx_NO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,21)*n(:,patmo_idx_HO2)*n(:,patmo_idx_NO) &
        + krate(:,32)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_NO) &
        - krate(:,55)*n(:,patmo_idx_SH)*n(:,patmo_idx_NO2) &
        - krate(:,59)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        - krate(:,64)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO2) &
        - krate(:,68)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO2) &
        - krate(:,140)*n(:,patmo_idx_CH)*n(:,patmo_idx_NO2) &
        - krate(:,171)*n(:,patmo_idx_H)*n(:,patmo_idx_NO2) &
        + krate(:,172)*n(:,patmo_idx_O)*n(:,patmo_idx_NO3) &
        - krate(:,188)*n(:,patmo_idx_CS)*n(:,patmo_idx_NO2) &
        - krate(:,197)*n(:,patmo_idx_NO2)*n(:,patmo_idx_N) &
        - krate(:,225)*n(:,patmo_idx_NO2) &
        + krate(:,227)*n(:,patmo_idx_NO3) &
        + krate(:,228)*n(:,patmo_idx_N2O5) &
        + krate(:,230)*n(:,patmo_idx_HNO3) &
        + krate(:,268)*n(:,patmo_idx_NO)*n(:,patmo_idx_O2) &
        - krate(:,269)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O2) &
        + krate(:,270)*n(:,patmo_idx_NO3)*n(:,patmo_idx_O2) &
        + krate(:,271)*n(:,patmo_idx_N2O5)*n(:,patmo_idx_M) &
        + krate(:,272)*n(:,patmo_idx_HNO3)*n(:,patmo_idx_M) &
        - krate(:,274)*n(:,patmo_idx_OH)*n(:,patmo_idx_NO2) &
        - krate(:,285)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_NO2) &
        + krate(:,308)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO) &
        + krate(:,312)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        + krate(:,317)*n(:,patmo_idx_SO3)*n(:,patmo_idx_NO) &
        + krate(:,321)*n(:,patmo_idx_NO)*n(:,patmo_idx_HSO2) &
        + krate(:,393)*n(:,patmo_idx_CHO)*n(:,patmo_idx_NO) &
        + krate(:,424)*n(:,patmo_idx_NO)*n(:,patmo_idx_OH) &
        - krate(:,425)*n(:,patmo_idx_O2)*n(:,patmo_idx_NO2) &
        + krate(:,441)*n(:,patmo_idx_COS)*n(:,patmo_idx_NO) &
        + krate(:,450)*n(:,patmo_idx_N2O)*n(:,patmo_idx_O)

    dn(:,patmo_idx_NO3) = &
        + krate(:,17)*n(:,patmo_idx_NO2)*n(:,patmo_idx_O3) &
        - krate(:,18)*n(:,patmo_idx_NO2)*n(:,patmo_idx_NO3)*n(:,patmo_idx_M) &
        + krate(:,20)*n(:,patmo_idx_HNO3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,172)*n(:,patmo_idx_O)*n(:,patmo_idx_NO3) &
        - krate(:,226)*n(:,patmo_idx_NO3) &
        - krate(:,227)*n(:,patmo_idx_NO3) &
        + krate(:,228)*n(:,patmo_idx_N2O5) &
        + krate(:,229)*n(:,patmo_idx_N2O5) &
        + krate(:,231)*n(:,patmo_idx_HNO3) &
        - krate(:,270)*n(:,patmo_idx_NO3)*n(:,patmo_idx_O2) &
        + krate(:,271)*n(:,patmo_idx_N2O5)*n(:,patmo_idx_M) &
        - krate(:,273)*n(:,patmo_idx_NO3)*n(:,patmo_idx_H2O)*n(:,patmo_idx_M) &
        + krate(:,425)*n(:,patmo_idx_O2)*n(:,patmo_idx_NO2)

    dn(:,patmo_idx_N2O5) = &
        + krate(:,18)*n(:,patmo_idx_NO2)*n(:,patmo_idx_NO3)*n(:,patmo_idx_M) &
        - krate(:,228)*n(:,patmo_idx_N2O5) &
        - krate(:,229)*n(:,patmo_idx_N2O5) &
        - krate(:,271)*n(:,patmo_idx_N2O5)*n(:,patmo_idx_M)

    dn(:,patmo_idx_HNO3) = &
        + krate(:,19)*n(:,patmo_idx_NO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,20)*n(:,patmo_idx_HNO3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,230)*n(:,patmo_idx_HNO3) &
        - krate(:,231)*n(:,patmo_idx_HNO3) &
        - krate(:,272)*n(:,patmo_idx_HNO3)*n(:,patmo_idx_M) &
        + krate(:,273)*n(:,patmo_idx_NO3)*n(:,patmo_idx_H2O)*n(:,patmo_idx_M)

    dn(:,patmo_idx_CH4) = &
        - krate(:,29)*n(:,patmo_idx_CH4)*n(:,patmo_idx_OH) &
        - krate(:,78)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,79)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,80)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        + krate(:,111)*n(:,patmo_idx_H)*n(:,patmo_idx_CH3)*n(:,patmo_idx_M) &
        + krate(:,113)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_CH3) &
        + krate(:,115)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        + krate(:,118)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        + krate(:,119)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH3) &
        + krate(:,120)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3) &
        + krate(:,123)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        + krate(:,124)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        + krate(:,126)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH3) &
        + krate(:,156)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        + krate(:,157)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        + krate(:,185)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2S) &
        + krate(:,203)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CH3) &
        - krate(:,211)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        - krate(:,232)*n(:,patmo_idx_CH4) &
        + krate(:,282)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        + krate(:,331)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,332)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        + krate(:,333)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        - krate(:,364)*n(:,patmo_idx_CH4)*n(:,patmo_idx_M) &
        - krate(:,366)*n(:,patmo_idx_CH4)*n(:,patmo_idx_HO2) &
        - krate(:,368)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O) &
        - krate(:,371)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O2) &
        - krate(:,372)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CO) &
        - krate(:,373)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        - krate(:,376)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        - krate(:,377)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        - krate(:,379)*n(:,patmo_idx_CH4)*n(:,patmo_idx_H) &
        - krate(:,409)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH3O) &
        - krate(:,410)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2OH) &
        - krate(:,438)*n(:,patmo_idx_CH4)*n(:,patmo_idx_SH) &
        - krate(:,456)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CO2) &
        + krate(:,464)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3)

    dn(:,patmo_idx_CH3) = &
        + krate(:,29)*n(:,patmo_idx_CH4)*n(:,patmo_idx_OH) &
        - krate(:,30)*n(:,patmo_idx_CH3)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,78)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        - krate(:,109)*n(:,patmo_idx_O)*n(:,patmo_idx_CH3) &
        - krate(:,110)*n(:,patmo_idx_O)*n(:,patmo_idx_CH3) &
        - krate(:,111)*n(:,patmo_idx_H)*n(:,patmo_idx_CH3)*n(:,patmo_idx_M) &
        - krate(:,112)*n(:,patmo_idx_O3)*n(:,patmo_idx_CH3) &
        - krate(:,113)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_CH3) &
        - krate(:,114)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        - krate(:,115)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        - krate(:,116)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3)*n(:,patmo_idx_M) &
        - krate(:,117)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        - krate(:,118)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        - krate(:,119)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH3) &
        - krate(:,120)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3) &
        - krate(:,120)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3) &
        - krate(:,121)*n(:,patmo_idx_CH3) &
        - krate(:,122)*n(:,patmo_idx_CH3) &
        - krate(:,123)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        - krate(:,124)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        - krate(:,125)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3) &
        - krate(:,126)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH3) &
        + krate(:,136)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH2) &
        + krate(:,145)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        + krate(:,146)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        + krate(:,147)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        + krate(:,150)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        - krate(:,156)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        - krate(:,157)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        + krate(:,158)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2) &
        + krate(:,160)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,185)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2S) &
        - krate(:,202)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CH3) &
        - krate(:,203)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CH3) &
        + krate(:,211)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        + krate(:,211)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        + krate(:,232)*n(:,patmo_idx_CH4) &
        + krate(:,247)*n(:,patmo_idx_CH3OH) &
        - krate(:,282)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        + krate(:,283)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_M) &
        - krate(:,331)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,362)*n(:,patmo_idx_CH3O) &
        + krate(:,363)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H) &
        + krate(:,364)*n(:,patmo_idx_CH4)*n(:,patmo_idx_M) &
        + krate(:,365)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        + krate(:,366)*n(:,patmo_idx_CH4)*n(:,patmo_idx_HO2) &
        + krate(:,367)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        + krate(:,368)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O) &
        + krate(:,369)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_M) &
        + krate(:,370)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        + krate(:,371)*n(:,patmo_idx_CH4)*n(:,patmo_idx_O2) &
        + krate(:,372)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CO) &
        + krate(:,373)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        + krate(:,373)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        + krate(:,374)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,375)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        + krate(:,376)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        + krate(:,377)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        + krate(:,378)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O) &
        + krate(:,379)*n(:,patmo_idx_CH4)*n(:,patmo_idx_H) &
        - krate(:,389)*n(:,patmo_idx_CO)*n(:,patmo_idx_CH3) &
        - krate(:,398)*n(:,patmo_idx_CH3) &
        - krate(:,399)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        - krate(:,400)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        - krate(:,403)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        + krate(:,409)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH3O) &
        + krate(:,410)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2OH) &
        - krate(:,411)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3) &
        - krate(:,413)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,438)*n(:,patmo_idx_CH4)*n(:,patmo_idx_SH) &
        + krate(:,455)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH2CO) &
        + krate(:,456)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CO2) &
        - krate(:,464)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3) &
        - krate(:,464)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3)

    dn(:,patmo_idx_CH3O2) = &
        + krate(:,30)*n(:,patmo_idx_CH3)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        - krate(:,31)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_HO2) &
        - krate(:,32)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_NO) &
        + krate(:,34)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        - krate(:,81)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        - krate(:,81)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        - krate(:,82)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        - krate(:,82)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        - krate(:,125)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3) &
        - krate(:,137)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH2) &
        - krate(:,283)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_M) &
        + krate(:,284)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_O2) &
        + krate(:,285)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_NO2) &
        - krate(:,287)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_H2O) &
        + krate(:,334)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        + krate(:,334)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        + krate(:,335)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_O2) &
        + krate(:,335)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_O2) &
        + krate(:,378)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O) &
        + krate(:,390)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3O)

    dn(:,patmo_idx_CH3OOH) = &
        + krate(:,31)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_HO2) &
        - krate(:,33)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        - krate(:,34)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        - krate(:,233)*n(:,patmo_idx_CH3OOH) &
        - krate(:,284)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_O2) &
        + krate(:,286)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2O) &
        + krate(:,287)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_H2O)

    dn(:,patmo_idx_CH3O) = &
        + krate(:,32)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_NO) &
        - krate(:,35)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        + krate(:,79)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        + krate(:,81)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        + krate(:,81)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        + krate(:,109)*n(:,patmo_idx_O)*n(:,patmo_idx_CH3) &
        + krate(:,112)*n(:,patmo_idx_O3)*n(:,patmo_idx_CH3) &
        + krate(:,114)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3) &
        + krate(:,117)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CH3) &
        - krate(:,123)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        + krate(:,125)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3) &
        + krate(:,125)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3) &
        + krate(:,137)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH2) &
        + krate(:,146)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        + krate(:,148)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        + krate(:,151)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        + krate(:,153)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        + krate(:,156)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        + krate(:,233)*n(:,patmo_idx_CH3OOH) &
        + krate(:,248)*n(:,patmo_idx_CH3OH) &
        - krate(:,285)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_NO2) &
        + krate(:,288)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_HO2) &
        - krate(:,332)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        - krate(:,334)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        - krate(:,334)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        - krate(:,362)*n(:,patmo_idx_CH3O) &
        - krate(:,365)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        - krate(:,367)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H) &
        - krate(:,370)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        + krate(:,376)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        - krate(:,378)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O) &
        - krate(:,378)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3O) &
        - krate(:,390)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3O) &
        - krate(:,399)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        - krate(:,401)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        - krate(:,404)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2) &
        - krate(:,406)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2O) &
        - krate(:,409)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH3O)

    dn(:,patmo_idx_CH2O) = &
        + krate(:,33)*n(:,patmo_idx_CH3OOH)*n(:,patmo_idx_OH) &
        + krate(:,35)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_O2) &
        - krate(:,36)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        + krate(:,80)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CH4) &
        + krate(:,82)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        + krate(:,86)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CHO) &
        - krate(:,90)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,91)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        - krate(:,101)*n(:,patmo_idx_SH)*n(:,patmo_idx_CH2O) &
        + krate(:,110)*n(:,patmo_idx_O)*n(:,patmo_idx_CH3) &
        + krate(:,123)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        + krate(:,124)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        + krate(:,134)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        + krate(:,135)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH2) &
        + krate(:,137)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH2) &
        + krate(:,138)*n(:,patmo_idx_CO2)*n(:,patmo_idx_CH2) &
        + krate(:,143)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH) &
        + krate(:,155)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        + krate(:,158)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2) &
        + krate(:,159)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_O) &
        + krate(:,162)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        + krate(:,164)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        + krate(:,165)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_HO2) &
        + krate(:,167)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        + krate(:,167)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        + krate(:,168)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2OH) &
        - krate(:,234)*n(:,patmo_idx_CH2O) &
        - krate(:,235)*n(:,patmo_idx_CH2O) &
        - krate(:,286)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH)*n(:,patmo_idx_H2O) &
        - krate(:,288)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_HO2) &
        + krate(:,289)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H2O) &
        - krate(:,333)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        - krate(:,335)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_O2) &
        - krate(:,339)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CO) &
        + krate(:,343)*n(:,patmo_idx_H2)*n(:,patmo_idx_CHO) &
        + krate(:,344)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        + krate(:,354)*n(:,patmo_idx_H2S)*n(:,patmo_idx_CHO) &
        - krate(:,363)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H) &
        - krate(:,376)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        - krate(:,377)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        - krate(:,387)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        - krate(:,388)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,390)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3O) &
        - krate(:,391)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CO) &
        - krate(:,396)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        - krate(:,408)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        - krate(:,411)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3) &
        - krate(:,412)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        - krate(:,415)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        - krate(:,417)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O) &
        - krate(:,418)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O2) &
        - krate(:,420)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH2O) &
        - krate(:,420)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH2O) &
        - krate(:,421)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3OH)

    dn(:,patmo_idx_CHO) = &
        + krate(:,36)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        - krate(:,37)*n(:,patmo_idx_CHO)*n(:,patmo_idx_O2) &
        + krate(:,84)*n(:,patmo_idx_H)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        - krate(:,85)*n(:,patmo_idx_H)*n(:,patmo_idx_CHO) &
        - krate(:,86)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CHO) &
        - krate(:,86)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CHO) &
        - krate(:,87)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        - krate(:,88)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        - krate(:,89)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        + krate(:,90)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        + krate(:,91)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        + krate(:,101)*n(:,patmo_idx_SH)*n(:,patmo_idx_CH2O) &
        - krate(:,119)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH3) &
        + krate(:,127)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        - krate(:,136)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH2) &
        + krate(:,140)*n(:,patmo_idx_CH)*n(:,patmo_idx_NO2) &
        + krate(:,141)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        - krate(:,166)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        - krate(:,167)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        + krate(:,234)*n(:,patmo_idx_CH2O) &
        - krate(:,236)*n(:,patmo_idx_CHO) &
        - krate(:,289)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H2O) &
        + krate(:,290)*n(:,patmo_idx_CO)*n(:,patmo_idx_HO2) &
        - krate(:,337)*n(:,patmo_idx_CHO)*n(:,patmo_idx_M) &
        + krate(:,338)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        + krate(:,339)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CO) &
        + krate(:,339)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CO) &
        + krate(:,340)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO) &
        + krate(:,341)*n(:,patmo_idx_H)*n(:,patmo_idx_CO2) &
        + krate(:,342)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        - krate(:,343)*n(:,patmo_idx_H2)*n(:,patmo_idx_CHO) &
        - krate(:,344)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        - krate(:,354)*n(:,patmo_idx_H2S)*n(:,patmo_idx_CHO) &
        + krate(:,372)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CO) &
        - krate(:,380)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H) &
        + krate(:,389)*n(:,patmo_idx_CO)*n(:,patmo_idx_CH3) &
        - krate(:,393)*n(:,patmo_idx_CHO)*n(:,patmo_idx_NO) &
        - krate(:,394)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        + krate(:,419)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CO) &
        + krate(:,420)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH2O)

    dn(:,patmo_idx_CO) = &
        + krate(:,37)*n(:,patmo_idx_CHO)*n(:,patmo_idx_O2) &
        - krate(:,38)*n(:,patmo_idx_CO)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,42)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,47)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,83)*n(:,patmo_idx_O)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        - krate(:,84)*n(:,patmo_idx_H)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        + krate(:,85)*n(:,patmo_idx_H)*n(:,patmo_idx_CHO) &
        + krate(:,86)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CHO) &
        + krate(:,87)*n(:,patmo_idx_OH)*n(:,patmo_idx_CHO) &
        + krate(:,89)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        + krate(:,119)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH3) &
        + krate(:,128)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        + krate(:,129)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        + krate(:,133)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        + krate(:,136)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH2) &
        + krate(:,138)*n(:,patmo_idx_CO2)*n(:,patmo_idx_CH2) &
        + krate(:,139)*n(:,patmo_idx_O)*n(:,patmo_idx_CH) &
        + krate(:,142)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        + krate(:,166)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        + krate(:,186)*n(:,patmo_idx_COS)*n(:,patmo_idx_H) &
        + krate(:,187)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        - krate(:,189)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        + krate(:,190)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,199)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        + krate(:,204)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_H) &
        - krate(:,206)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CO) &
        + krate(:,235)*n(:,patmo_idx_CH2O) &
        + krate(:,236)*n(:,patmo_idx_CHO) &
        + krate(:,237)*n(:,patmo_idx_CO2) &
        + krate(:,240)*n(:,patmo_idx_COS) &
        - krate(:,290)*n(:,patmo_idx_CO)*n(:,patmo_idx_HO2) &
        + krate(:,291)*n(:,patmo_idx_CO2)*n(:,patmo_idx_H)*n(:,patmo_idx_M) &
        - krate(:,295)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,300)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,336)*n(:,patmo_idx_CO2)*n(:,patmo_idx_M) &
        + krate(:,337)*n(:,patmo_idx_CHO)*n(:,patmo_idx_M) &
        - krate(:,338)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        - krate(:,339)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CO) &
        - krate(:,340)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO) &
        - krate(:,342)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        - krate(:,372)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CO) &
        - krate(:,381)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        - krate(:,382)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        - krate(:,386)*n(:,patmo_idx_CO)*n(:,patmo_idx_H2O) &
        - krate(:,389)*n(:,patmo_idx_CO)*n(:,patmo_idx_CH3) &
        - krate(:,391)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CO) &
        - krate(:,392)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        - krate(:,395)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        - krate(:,419)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CO) &
        - krate(:,439)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        - krate(:,440)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        + krate(:,442)*n(:,patmo_idx_COS)*n(:,patmo_idx_H) &
        - krate(:,443)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        + krate(:,452)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_M) &
        - krate(:,457)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO) &
        + krate(:,459)*n(:,patmo_idx_COCOOH)

    dn(:,patmo_idx_H2O2) = &
        + krate(:,39)*n(:,patmo_idx_HO2)*n(:,patmo_idx_HO2)*n(:,patmo_idx_M) &
        - krate(:,40)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_OH) &
        - krate(:,113)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_CH3) &
        - krate(:,163)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O2) &
        + krate(:,165)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_HO2) &
        + krate(:,207)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,238)*n(:,patmo_idx_H2O2) &
        - krate(:,239)*n(:,patmo_idx_H2O2) &
        - krate(:,292)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_M) &
        + krate(:,293)*n(:,patmo_idx_HO2)*n(:,patmo_idx_H2O) &
        + krate(:,366)*n(:,patmo_idx_CH4)*n(:,patmo_idx_HO2) &
        + krate(:,416)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_HO2) &
        - krate(:,418)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O2) &
        - krate(:,460)*n(:,patmo_idx_H2O2)*n(:,patmo_idx_M)

    dn(:,patmo_idx_COS) = &
        - krate(:,41)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,42)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,43)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,45)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,46)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,186)*n(:,patmo_idx_COS)*n(:,patmo_idx_H) &
        - krate(:,187)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        + krate(:,188)*n(:,patmo_idx_CS)*n(:,patmo_idx_NO2) &
        + krate(:,189)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        + krate(:,191)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,240)*n(:,patmo_idx_COS) &
        + krate(:,294)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,295)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,296)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,298)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,299)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,439)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        + krate(:,440)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        - krate(:,441)*n(:,patmo_idx_COS)*n(:,patmo_idx_NO) &
        - krate(:,442)*n(:,patmo_idx_COS)*n(:,patmo_idx_H) &
        - krate(:,444)*n(:,patmo_idx_COS)*n(:,patmo_idx_S)

    dn(:,patmo_idx_SH) = &
        + krate(:,41)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        + krate(:,43)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,48)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,49)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,50)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,52)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        - krate(:,53)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,54)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,55)*n(:,patmo_idx_SH)*n(:,patmo_idx_NO2) &
        + krate(:,67)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,99)*n(:,patmo_idx_SH)*n(:,patmo_idx_SH) &
        - krate(:,99)*n(:,patmo_idx_SH)*n(:,patmo_idx_SH) &
        - krate(:,100)*n(:,patmo_idx_SH)*n(:,patmo_idx_H) &
        - krate(:,101)*n(:,patmo_idx_SH)*n(:,patmo_idx_CH2O) &
        + krate(:,185)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2S) &
        + krate(:,186)*n(:,patmo_idx_COS)*n(:,patmo_idx_H) &
        - krate(:,189)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        + krate(:,243)*n(:,patmo_idx_H2S) &
        - krate(:,294)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        - krate(:,296)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,301)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,302)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,303)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,305)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        + krate(:,306)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,307)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,308)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO) &
        - krate(:,320)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        + krate(:,352)*n(:,patmo_idx_S)*n(:,patmo_idx_H2S) &
        + krate(:,352)*n(:,patmo_idx_S)*n(:,patmo_idx_H2S) &
        + krate(:,353)*n(:,patmo_idx_H2)*n(:,patmo_idx_S) &
        + krate(:,354)*n(:,patmo_idx_H2S)*n(:,patmo_idx_CHO) &
        - krate(:,438)*n(:,patmo_idx_CH4)*n(:,patmo_idx_SH) &
        - krate(:,439)*n(:,patmo_idx_CO)*n(:,patmo_idx_SH) &
        + krate(:,442)*n(:,patmo_idx_COS)*n(:,patmo_idx_H)

    dn(:,patmo_idx_SO) = &
        + krate(:,42)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,44)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,52)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,53)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,58)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,59)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        + krate(:,60)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,61)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,62)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        - krate(:,94)*n(:,patmo_idx_SO)*n(:,patmo_idx_HO2) &
        - krate(:,95)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO)*n(:,patmo_idx_M) &
        - krate(:,95)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO)*n(:,patmo_idx_M) &
        - krate(:,96)*n(:,patmo_idx_SO)*n(:,patmo_idx_S2O2) &
        - krate(:,97)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO) &
        - krate(:,97)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO) &
        - krate(:,98)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO3) &
        + krate(:,108)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,241)*n(:,patmo_idx_SO) &
        + krate(:,244)*n(:,patmo_idx_SO2) &
        + krate(:,249)*n(:,patmo_idx_S2O2) &
        + krate(:,249)*n(:,patmo_idx_S2O2) &
        + krate(:,250)*n(:,patmo_idx_S2O) &
        - krate(:,295)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,297)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,305)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,306)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,309)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,310)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,311)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,312)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        - krate(:,313)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,314)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        - krate(:,315)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        + krate(:,347)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,348)*n(:,patmo_idx_S2O2)*n(:,patmo_idx_M) &
        + krate(:,348)*n(:,patmo_idx_S2O2)*n(:,patmo_idx_M) &
        + krate(:,349)*n(:,patmo_idx_SO2)*n(:,patmo_idx_S2O) &
        + krate(:,350)*n(:,patmo_idx_S)*n(:,patmo_idx_SO2) &
        + krate(:,350)*n(:,patmo_idx_S)*n(:,patmo_idx_SO2) &
        + krate(:,351)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SO2) &
        - krate(:,361)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_CS2) = &
        - krate(:,43)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,44)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,190)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,191)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,242)*n(:,patmo_idx_CS2) &
        + krate(:,296)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,297)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,443)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        + krate(:,444)*n(:,patmo_idx_COS)*n(:,patmo_idx_S)

    dn(:,patmo_idx_CS) = &
        + krate(:,44)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,45)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,46)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,47)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,188)*n(:,patmo_idx_CS)*n(:,patmo_idx_NO2) &
        + krate(:,242)*n(:,patmo_idx_CS2) &
        - krate(:,297)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,298)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,299)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,300)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,441)*n(:,patmo_idx_COS)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_S) = &
        + krate(:,47)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,60)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,61)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,62)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,97)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO) &
        + krate(:,99)*n(:,patmo_idx_SH)*n(:,patmo_idx_SH) &
        + krate(:,100)*n(:,patmo_idx_SH)*n(:,patmo_idx_H) &
        - krate(:,102)*n(:,patmo_idx_S)*n(:,patmo_idx_S)*n(:,patmo_idx_M) &
        - krate(:,102)*n(:,patmo_idx_S)*n(:,patmo_idx_S)*n(:,patmo_idx_M) &
        - krate(:,103)*n(:,patmo_idx_S)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        - krate(:,104)*n(:,patmo_idx_S)*n(:,patmo_idx_S3)*n(:,patmo_idx_M) &
        + krate(:,107)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        + krate(:,107)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        + krate(:,108)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,187)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        + krate(:,191)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,240)*n(:,patmo_idx_COS) &
        + krate(:,241)*n(:,patmo_idx_SO) &
        + krate(:,242)*n(:,patmo_idx_CS2) &
        + krate(:,250)*n(:,patmo_idx_S2O) &
        - krate(:,300)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,313)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,314)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        + krate(:,315)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,350)*n(:,patmo_idx_S)*n(:,patmo_idx_SO2) &
        - krate(:,352)*n(:,patmo_idx_S)*n(:,patmo_idx_H2S) &
        - krate(:,353)*n(:,patmo_idx_H2)*n(:,patmo_idx_S) &
        + krate(:,355)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        + krate(:,355)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        + krate(:,356)*n(:,patmo_idx_S3)*n(:,patmo_idx_M) &
        + krate(:,357)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        - krate(:,360)*n(:,patmo_idx_S)*n(:,patmo_idx_S)*n(:,patmo_idx_M) &
        - krate(:,360)*n(:,patmo_idx_S)*n(:,patmo_idx_S)*n(:,patmo_idx_M) &
        - krate(:,361)*n(:,patmo_idx_S)*n(:,patmo_idx_SO) &
        + krate(:,440)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        - krate(:,444)*n(:,patmo_idx_COS)*n(:,patmo_idx_S)

    dn(:,patmo_idx_H2S) = &
        - krate(:,48)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        - krate(:,49)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,50)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,51)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        + krate(:,99)*n(:,patmo_idx_SH)*n(:,patmo_idx_SH) &
        + krate(:,101)*n(:,patmo_idx_SH)*n(:,patmo_idx_CH2O) &
        - krate(:,185)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2S) &
        - krate(:,243)*n(:,patmo_idx_H2S) &
        + krate(:,301)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        + krate(:,302)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,303)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,304)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        - krate(:,352)*n(:,patmo_idx_S)*n(:,patmo_idx_H2S) &
        - krate(:,354)*n(:,patmo_idx_H2S)*n(:,patmo_idx_CHO) &
        + krate(:,438)*n(:,patmo_idx_CH4)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_HSO) = &
        + krate(:,51)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        + krate(:,54)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,55)*n(:,patmo_idx_SH)*n(:,patmo_idx_NO2) &
        - krate(:,66)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,67)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,68)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO2) &
        - krate(:,304)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        - krate(:,307)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,308)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO) &
        + krate(:,319)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,320)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
        + krate(:,321)*n(:,patmo_idx_NO)*n(:,patmo_idx_HSO2)

    dn(:,patmo_idx_SO2) = &
        + krate(:,56)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,57)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,58)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,59)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        - krate(:,64)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO2) &
        - krate(:,65)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,66)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,69)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,72)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,73)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,74)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,75)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        + krate(:,76)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,94)*n(:,patmo_idx_SO)*n(:,patmo_idx_HO2) &
        + krate(:,96)*n(:,patmo_idx_SO)*n(:,patmo_idx_S2O2) &
        + krate(:,97)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO) &
        + krate(:,98)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO3) &
        + krate(:,98)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO3) &
        - krate(:,210)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_SO2) &
        + krate(:,210)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_SO2) &
        - krate(:,244)*n(:,patmo_idx_SO2) &
        + krate(:,245)*n(:,patmo_idx_SO3) &
        + krate(:,246)*n(:,patmo_idx_H2SO4) &
        - krate(:,309)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,310)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,311)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,312)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        + krate(:,316)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        + krate(:,317)*n(:,patmo_idx_SO3)*n(:,patmo_idx_NO) &
        + krate(:,318)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,319)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,322)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        + krate(:,325)*n(:,patmo_idx_SO3)*n(:,patmo_idx_M) &
        + krate(:,326)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_M) &
        - krate(:,327)*n(:,patmo_idx_SO2) &
        - krate(:,328)*n(:,patmo_idx_SO2) &
        - krate(:,329)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)*n(:,patmo_idx_M) &
        - krate(:,347)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,349)*n(:,patmo_idx_SO2)*n(:,patmo_idx_S2O) &
        - krate(:,350)*n(:,patmo_idx_S)*n(:,patmo_idx_SO2) &
        - krate(:,351)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SO2) &
        - krate(:,351)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SO2) &
        - krate(:,463)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_SO2) &
        + krate(:,463)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_SO3) = &
        + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,64)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO2) &
        + krate(:,65)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,70)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        - krate(:,71)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        + krate(:,72)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,98)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO3) &
        - krate(:,245)*n(:,patmo_idx_SO3) &
        - krate(:,316)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,317)*n(:,patmo_idx_SO3)*n(:,patmo_idx_NO) &
        - krate(:,318)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,323)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        + krate(:,324)*n(:,patmo_idx_H2SO4) &
        - krate(:,325)*n(:,patmo_idx_SO3)*n(:,patmo_idx_M) &
        + krate(:,351)*n(:,patmo_idx_SO2)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_HSO2) = &
        + krate(:,68)*n(:,patmo_idx_HSO)*n(:,patmo_idx_NO2) &
        - krate(:,69)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,321)*n(:,patmo_idx_NO)*n(:,patmo_idx_HSO2) &
        + krate(:,322)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_HSO3) = &
        - krate(:,70)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,73)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,323)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,326)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_M)

    dn(:,patmo_idx_H2SO4) = &
        + krate(:,71)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,77)*n(:,patmo_idx_H2SO4) &
        - krate(:,246)*n(:,patmo_idx_H2SO4) &
        - krate(:,324)*n(:,patmo_idx_H2SO4) &
        + krate(:,330)*n(:,patmo_idx_SO4)

    dn(:,patmo_idx_CH3SCH3) = &
        - krate(:,74)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,75)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,76)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,327)*n(:,patmo_idx_SO2) &
        + krate(:,328)*n(:,patmo_idx_SO2) &
        + krate(:,329)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)*n(:,patmo_idx_M)

    dn(:,patmo_idx_CH4O3S) = &
        + krate(:,76)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,329)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)*n(:,patmo_idx_M)

    dn(:,patmo_idx_SO4) = &
        + krate(:,77)*n(:,patmo_idx_H2SO4) &
        - krate(:,330)*n(:,patmo_idx_SO4)

    dn(:,patmo_idx_CH3OH) = &
        + krate(:,82)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH3O2) &
        + krate(:,116)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH3)*n(:,patmo_idx_M) &
        - krate(:,146)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        - krate(:,147)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        - krate(:,148)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        - krate(:,149)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        - krate(:,150)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        - krate(:,151)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        - krate(:,152)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        - krate(:,153)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        - krate(:,154)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        - krate(:,155)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        - krate(:,156)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        - krate(:,157)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        + krate(:,161)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        + krate(:,163)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O2) &
        + krate(:,166)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        + krate(:,168)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2OH) &
        - krate(:,247)*n(:,patmo_idx_CH3OH) &
        - krate(:,248)*n(:,patmo_idx_CH3OH) &
        - krate(:,335)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_O2) &
        - krate(:,369)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_M) &
        + krate(:,399)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        + krate(:,400)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        + krate(:,401)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_OH) &
        + krate(:,402)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        + krate(:,403)*n(:,patmo_idx_CH3)*n(:,patmo_idx_H2O) &
        + krate(:,404)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2) &
        + krate(:,405)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2) &
        + krate(:,406)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_H2O) &
        + krate(:,407)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O) &
        + krate(:,408)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O)*n(:,patmo_idx_H) &
        + krate(:,409)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH3O) &
        + krate(:,410)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2OH) &
        - krate(:,414)*n(:,patmo_idx_CH3OH) &
        - krate(:,416)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_HO2) &
        - krate(:,419)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CO) &
        - krate(:,421)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3OH)

    dn(:,patmo_idx_S2O2) = &
        + krate(:,95)*n(:,patmo_idx_SO)*n(:,patmo_idx_SO)*n(:,patmo_idx_M) &
        - krate(:,96)*n(:,patmo_idx_SO)*n(:,patmo_idx_S2O2) &
        - krate(:,249)*n(:,patmo_idx_S2O2) &
        - krate(:,348)*n(:,patmo_idx_S2O2)*n(:,patmo_idx_M) &
        + krate(:,349)*n(:,patmo_idx_SO2)*n(:,patmo_idx_S2O)

    dn(:,patmo_idx_S2O) = &
        + krate(:,96)*n(:,patmo_idx_SO)*n(:,patmo_idx_S2O2) &
        - krate(:,250)*n(:,patmo_idx_S2O) &
        - krate(:,349)*n(:,patmo_idx_SO2)*n(:,patmo_idx_S2O)

    dn(:,patmo_idx_S2) = &
        + krate(:,102)*n(:,patmo_idx_S)*n(:,patmo_idx_S)*n(:,patmo_idx_M) &
        - krate(:,103)*n(:,patmo_idx_S)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        - krate(:,105)*n(:,patmo_idx_S2)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        - krate(:,105)*n(:,patmo_idx_S2)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        - krate(:,107)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        - krate(:,108)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        + krate(:,187)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        + krate(:,190)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,355)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        + krate(:,356)*n(:,patmo_idx_S3)*n(:,patmo_idx_M) &
        + krate(:,358)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        + krate(:,358)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        + krate(:,360)*n(:,patmo_idx_S)*n(:,patmo_idx_S)*n(:,patmo_idx_M) &
        + krate(:,361)*n(:,patmo_idx_S)*n(:,patmo_idx_SO) &
        - krate(:,440)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2) &
        - krate(:,443)*n(:,patmo_idx_CO)*n(:,patmo_idx_S2)

    dn(:,patmo_idx_S3) = &
        + krate(:,103)*n(:,patmo_idx_S)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        - krate(:,104)*n(:,patmo_idx_S)*n(:,patmo_idx_S3)*n(:,patmo_idx_M) &
        - krate(:,356)*n(:,patmo_idx_S3)*n(:,patmo_idx_M) &
        + krate(:,357)*n(:,patmo_idx_S4)*n(:,patmo_idx_M)

    dn(:,patmo_idx_S4) = &
        + krate(:,104)*n(:,patmo_idx_S)*n(:,patmo_idx_S3)*n(:,patmo_idx_M) &
        + krate(:,105)*n(:,patmo_idx_S2)*n(:,patmo_idx_S2)*n(:,patmo_idx_M) &
        - krate(:,106)*n(:,patmo_idx_S4)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        - krate(:,106)*n(:,patmo_idx_S4)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        - krate(:,357)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        - krate(:,358)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        + krate(:,359)*n(:,patmo_idx_S8)*n(:,patmo_idx_M) &
        + krate(:,359)*n(:,patmo_idx_S8)*n(:,patmo_idx_M)

    dn(:,patmo_idx_S8) = &
        + krate(:,106)*n(:,patmo_idx_S4)*n(:,patmo_idx_S4)*n(:,patmo_idx_M) &
        - krate(:,359)*n(:,patmo_idx_S8)*n(:,patmo_idx_M)

    dn(:,patmo_idx_CH2) = &
        + krate(:,120)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3) &
        + krate(:,121)*n(:,patmo_idx_CH3) &
        - krate(:,127)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        - krate(:,128)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        - krate(:,129)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2) &
        - krate(:,130)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        - krate(:,131)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,132)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,133)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,134)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH2) &
        - krate(:,135)*n(:,patmo_idx_OH)*n(:,patmo_idx_CH2) &
        - krate(:,136)*n(:,patmo_idx_CHO)*n(:,patmo_idx_CH2) &
        - krate(:,137)*n(:,patmo_idx_CH3O2)*n(:,patmo_idx_CH2) &
        - krate(:,138)*n(:,patmo_idx_CO2)*n(:,patmo_idx_CH2) &
        + krate(:,144)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,146)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        - krate(:,147)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        - krate(:,158)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2) &
        - krate(:,211)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        - krate(:,373)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2) &
        - krate(:,374)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,380)*n(:,patmo_idx_CHO)*n(:,patmo_idx_H) &
        + krate(:,381)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        + krate(:,382)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO) &
        + krate(:,383)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        + krate(:,384)*n(:,patmo_idx_H)*n(:,patmo_idx_H)*n(:,patmo_idx_CO2) &
        + krate(:,385)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO2) &
        + krate(:,386)*n(:,patmo_idx_CO)*n(:,patmo_idx_H2O) &
        + krate(:,387)*n(:,patmo_idx_O)*n(:,patmo_idx_CH2O) &
        + krate(:,388)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        + krate(:,389)*n(:,patmo_idx_CO)*n(:,patmo_idx_CH3) &
        + krate(:,390)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3O) &
        + krate(:,391)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CO) &
        - krate(:,397)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,399)*n(:,patmo_idx_CH3O)*n(:,patmo_idx_CH3) &
        + krate(:,400)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        + krate(:,411)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3) &
        + krate(:,464)*n(:,patmo_idx_CH3)*n(:,patmo_idx_CH3)

    dn(:,patmo_idx_CH) = &
        + krate(:,122)*n(:,patmo_idx_CH3) &
        + krate(:,130)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        - krate(:,139)*n(:,patmo_idx_O)*n(:,patmo_idx_CH) &
        - krate(:,140)*n(:,patmo_idx_CH)*n(:,patmo_idx_NO2) &
        - krate(:,141)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        - krate(:,142)*n(:,patmo_idx_O2)*n(:,patmo_idx_CH) &
        - krate(:,143)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH) &
        - krate(:,144)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,145)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,375)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        - krate(:,383)*n(:,patmo_idx_H2)*n(:,patmo_idx_CH) &
        + krate(:,392)*n(:,patmo_idx_H)*n(:,patmo_idx_CO) &
        + krate(:,393)*n(:,patmo_idx_CHO)*n(:,patmo_idx_NO) &
        + krate(:,394)*n(:,patmo_idx_O)*n(:,patmo_idx_CHO) &
        + krate(:,395)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO) &
        + krate(:,396)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2O) &
        + krate(:,397)*n(:,patmo_idx_H)*n(:,patmo_idx_CH2) &
        + krate(:,398)*n(:,patmo_idx_CH3)

    dn(:,patmo_idx_CH2OH) = &
        - krate(:,124)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        + krate(:,147)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH2) &
        + krate(:,149)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_O) &
        + krate(:,152)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_H) &
        + krate(:,154)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_OH) &
        + krate(:,157)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CH3) &
        - krate(:,158)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2) &
        - krate(:,159)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_O) &
        - krate(:,160)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,161)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,162)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H) &
        - krate(:,163)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O2) &
        - krate(:,164)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        - krate(:,165)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_HO2) &
        - krate(:,166)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        - krate(:,167)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CHO) &
        - krate(:,168)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2OH) &
        - krate(:,168)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH2OH) &
        + krate(:,377)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH4) &
        - krate(:,400)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_CH3) &
        - krate(:,402)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_OH) &
        - krate(:,405)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2) &
        - krate(:,407)*n(:,patmo_idx_CH2OH)*n(:,patmo_idx_H2O) &
        - krate(:,410)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CH2OH) &
        + krate(:,411)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3) &
        + krate(:,412)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_OH) &
        + krate(:,413)*n(:,patmo_idx_CH3)*n(:,patmo_idx_OH) &
        + krate(:,414)*n(:,patmo_idx_CH3OH) &
        + krate(:,415)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2) &
        + krate(:,416)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_HO2) &
        + krate(:,417)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O) &
        + krate(:,418)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_H2O2) &
        + krate(:,419)*n(:,patmo_idx_CH3OH)*n(:,patmo_idx_CO) &
        + krate(:,420)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH2O) &
        + krate(:,421)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3OH) &
        + krate(:,421)*n(:,patmo_idx_CH2O)*n(:,patmo_idx_CH3OH)

    dn(:,patmo_idx_N) = &
        - krate(:,169)*n(:,patmo_idx_N)*n(:,patmo_idx_O2) &
        - krate(:,170)*n(:,patmo_idx_N)*n(:,patmo_idx_NO) &
        + krate(:,177)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        + krate(:,193)*n(:,patmo_idx_NH)*n(:,patmo_idx_NH) &
        + krate(:,194)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH) &
        - krate(:,195)*n(:,patmo_idx_O)*n(:,patmo_idx_N)*n(:,patmo_idx_M) &
        - krate(:,196)*n(:,patmo_idx_H)*n(:,patmo_idx_N)*n(:,patmo_idx_M) &
        - krate(:,197)*n(:,patmo_idx_NO2)*n(:,patmo_idx_N) &
        + krate(:,422)*n(:,patmo_idx_O)*n(:,patmo_idx_NO) &
        + krate(:,423)*n(:,patmo_idx_N2)*n(:,patmo_idx_O) &
        - krate(:,430)*n(:,patmo_idx_N)*n(:,patmo_idx_OH) &
        - krate(:,446)*n(:,patmo_idx_NH2)*n(:,patmo_idx_N) &
        - krate(:,447)*n(:,patmo_idx_NH3)*n(:,patmo_idx_N) &
        + krate(:,448)*n(:,patmo_idx_NO)*n(:,patmo_idx_M) &
        + krate(:,449)*n(:,patmo_idx_NH)*n(:,patmo_idx_M) &
        + krate(:,450)*n(:,patmo_idx_N2O)*n(:,patmo_idx_O)

    dn(:,patmo_idx_NH2) = &
        - krate(:,173)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH2)*n(:,patmo_idx_M) &
        - krate(:,173)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH2)*n(:,patmo_idx_M) &
        + krate(:,175)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H) &
        + krate(:,175)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H) &
        - krate(:,178)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NO) &
        - krate(:,179)*n(:,patmo_idx_NH2)*n(:,patmo_idx_O) &
        + krate(:,180)*n(:,patmo_idx_NH3)*n(:,patmo_idx_O_1D) &
        + krate(:,181)*n(:,patmo_idx_NH3)*n(:,patmo_idx_OH) &
        - krate(:,182)*n(:,patmo_idx_NH2)*n(:,patmo_idx_H)*n(:,patmo_idx_M) &
        - krate(:,192)*n(:,patmo_idx_OH)*n(:,patmo_idx_NH2) &
        + krate(:,193)*n(:,patmo_idx_NH)*n(:,patmo_idx_NH) &
        - krate(:,194)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH) &
        + krate(:,252)*n(:,patmo_idx_NH3) &
        + krate(:,426)*n(:,patmo_idx_N2H4)*n(:,patmo_idx_M) &
        + krate(:,426)*n(:,patmo_idx_N2H4)*n(:,patmo_idx_M) &
        - krate(:,428)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH2) &
        - krate(:,428)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH2) &
        + krate(:,431)*n(:,patmo_idx_N2)*n(:,patmo_idx_H2O) &
        + krate(:,432)*n(:,patmo_idx_NH)*n(:,patmo_idx_OH) &
        - krate(:,433)*n(:,patmo_idx_NH2)*n(:,patmo_idx_OH) &
        - krate(:,434)*n(:,patmo_idx_NH2)*n(:,patmo_idx_H2O) &
        + krate(:,435)*n(:,patmo_idx_NH3)*n(:,patmo_idx_M) &
        + krate(:,445)*n(:,patmo_idx_H2O)*n(:,patmo_idx_NH) &
        - krate(:,446)*n(:,patmo_idx_NH2)*n(:,patmo_idx_N) &
        + krate(:,447)*n(:,patmo_idx_NH3)*n(:,patmo_idx_N)

    dn(:,patmo_idx_N2H4) = &
        + krate(:,173)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH2)*n(:,patmo_idx_M) &
        - krate(:,174)*n(:,patmo_idx_N2H4)*n(:,patmo_idx_H) &
        - krate(:,251)*n(:,patmo_idx_N2H4) &
        - krate(:,426)*n(:,patmo_idx_N2H4)*n(:,patmo_idx_M) &
        + krate(:,427)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H2)

    dn(:,patmo_idx_N2H3) = &
        + krate(:,174)*n(:,patmo_idx_N2H4)*n(:,patmo_idx_H) &
        - krate(:,175)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H) &
        + krate(:,251)*n(:,patmo_idx_N2H4) &
        - krate(:,427)*n(:,patmo_idx_N2H3)*n(:,patmo_idx_H2) &
        + krate(:,428)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH2)

    dn(:,patmo_idx_NH) = &
        - krate(:,176)*n(:,patmo_idx_NH)*n(:,patmo_idx_NO) &
        - krate(:,177)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        + krate(:,179)*n(:,patmo_idx_NH2)*n(:,patmo_idx_O) &
        - krate(:,183)*n(:,patmo_idx_NH)*n(:,patmo_idx_NO) &
        - krate(:,184)*n(:,patmo_idx_NH)*n(:,patmo_idx_O) &
        + krate(:,192)*n(:,patmo_idx_OH)*n(:,patmo_idx_NH2) &
        - krate(:,193)*n(:,patmo_idx_NH)*n(:,patmo_idx_NH) &
        - krate(:,193)*n(:,patmo_idx_NH)*n(:,patmo_idx_NH) &
        - krate(:,194)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH) &
        + krate(:,196)*n(:,patmo_idx_H)*n(:,patmo_idx_N)*n(:,patmo_idx_M) &
        + krate(:,253)*n(:,patmo_idx_NH3) &
        + krate(:,429)*n(:,patmo_idx_N2)*n(:,patmo_idx_OH) &
        + krate(:,430)*n(:,patmo_idx_N)*n(:,patmo_idx_OH) &
        - krate(:,432)*n(:,patmo_idx_NH)*n(:,patmo_idx_OH) &
        + krate(:,436)*n(:,patmo_idx_N2O)*n(:,patmo_idx_H) &
        + krate(:,437)*n(:,patmo_idx_NO)*n(:,patmo_idx_H) &
        - krate(:,445)*n(:,patmo_idx_H2O)*n(:,patmo_idx_NH) &
        + krate(:,446)*n(:,patmo_idx_NH2)*n(:,patmo_idx_N) &
        + krate(:,446)*n(:,patmo_idx_NH2)*n(:,patmo_idx_N) &
        + krate(:,447)*n(:,patmo_idx_NH3)*n(:,patmo_idx_N) &
        - krate(:,449)*n(:,patmo_idx_NH)*n(:,patmo_idx_M)

    dn(:,patmo_idx_NH3) = &
        - krate(:,180)*n(:,patmo_idx_NH3)*n(:,patmo_idx_O_1D) &
        - krate(:,181)*n(:,patmo_idx_NH3)*n(:,patmo_idx_OH) &
        + krate(:,182)*n(:,patmo_idx_NH2)*n(:,patmo_idx_H)*n(:,patmo_idx_M) &
        + krate(:,194)*n(:,patmo_idx_NH2)*n(:,patmo_idx_NH) &
        - krate(:,252)*n(:,patmo_idx_NH3) &
        - krate(:,253)*n(:,patmo_idx_NH3) &
        + krate(:,433)*n(:,patmo_idx_NH2)*n(:,patmo_idx_OH) &
        + krate(:,434)*n(:,patmo_idx_NH2)*n(:,patmo_idx_H2O) &
        - krate(:,435)*n(:,patmo_idx_NH3)*n(:,patmo_idx_M) &
        - krate(:,447)*n(:,patmo_idx_NH3)*n(:,patmo_idx_N)

    dn(:,patmo_idx_HOCO) = &
        + krate(:,199)*n(:,patmo_idx_OH)*n(:,patmo_idx_CO)*n(:,patmo_idx_M) &
        - krate(:,200)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_O_3P) &
        - krate(:,201)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_OH) &
        - krate(:,202)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CH3) &
        - krate(:,203)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CH3) &
        - krate(:,204)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_H) &
        - krate(:,205)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_H) &
        - krate(:,206)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CO) &
        - krate(:,214)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_O2) &
        - krate(:,452)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_M) &
        + krate(:,453)*n(:,patmo_idx_CO2)*n(:,patmo_idx_OH) &
        + krate(:,454)*n(:,patmo_idx_CO2)*n(:,patmo_idx_H2O) &
        + krate(:,455)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH2CO) &
        + krate(:,456)*n(:,patmo_idx_CH4)*n(:,patmo_idx_CO2) &
        + krate(:,457)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CO) &
        + krate(:,458)*n(:,patmo_idx_H2)*n(:,patmo_idx_CO2) &
        + krate(:,459)*n(:,patmo_idx_COCOOH) &
        + krate(:,467)*n(:,patmo_idx_HO2)*n(:,patmo_idx_CO2)

    dn(:,patmo_idx_O_3P) = &
        - krate(:,200)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_O_3P) &
        + krate(:,208)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_CO2) &
        + krate(:,209)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_N2) &
        + krate(:,210)*n(:,patmo_idx_O_1D)*n(:,patmo_idx_SO2) &
        + krate(:,453)*n(:,patmo_idx_CO2)*n(:,patmo_idx_OH) &
        - krate(:,461)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_CO2) &
        - krate(:,462)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_N2) &
        - krate(:,463)*n(:,patmo_idx_O_3P)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_CH2CO) = &
        + krate(:,202)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CH3) &
        - krate(:,455)*n(:,patmo_idx_H2O)*n(:,patmo_idx_CH2CO)

    dn(:,patmo_idx_COCOOH) = &
        + krate(:,206)*n(:,patmo_idx_HOCO)*n(:,patmo_idx_CO) &
        - krate(:,459)*n(:,patmo_idx_COCOOH)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

    do i=1,chemSpeciesNumber
      dn(:,i) = dn(:,i) &
          + (k_hp(:)-d_hp(:,i)) * ngas_hpp(:) * n_p(:,i) &
          - ((k_hp(:)+d_hp(:,i)) * ngas_hpz(:) &
          + (k_hm(:)-d_hm(:,i)) * ngas_hmz(:)) * n(:,i) &
          + (k_hm(:)+d_hm(:,i)) * ngas_hmm(:) * n_m(:,i)
    end do

    !Chemical Species with constant concentration
    dn(:,patmo_idx_CO2) = 0d0
    dn(:,patmo_idx_N2) = 0d0
    dn(:,patmo_idx_M) = 0d0

    ! Dry Deposition: assumed a deposition rate of 0.1 cm/s
    !dn(1,patmo_idx_A)=dn(1,patmo_idx_A) - 0.1/1.0d5*n(1,patmo_idx_A)
    ! Fix the mixing ratio of CH4 and O2 at the bottom layer as a constant (Claire et al., 2014; Zahnle et al., 2006)
    if (n(1,patmo_idx_CH4) > 0.1/1d5) then
      dn(1,patmo_idx_CH4) = dn(1,patmo_idx_CH4) - (0.1/1d5) * n(1,patmo_idx_CH4)
    end if
    CH4Flux = -1d5 * dn(1, patmo_idx_CH4)
    dn(1,patmo_idx_CH4) = 0d0
    if (n(1,patmo_idx_O2) > 0.1/1d5) then
      dn(1,patmo_idx_O2) = dn(1,patmo_idx_O2) - (0.1/1d5) * n(1,patmo_idx_O2)
    end if
    O2Flux = -1d5 * dn(1, patmo_idx_O2)
    dn(1,patmo_idx_O2) = 0d0
    if (n(1,patmo_idx_S8) > 0.1/1d5) then
      dn(1,patmo_idx_S8) = dn(1,patmo_idx_S8) - (0.1/1d5) * n(1,patmo_idx_S8)
    end if
    if (n(1,patmo_idx_H2O2) > 0.1/1d5) then
      dn(1,patmo_idx_H2O2) = dn(1,patmo_idx_H2O2) - (0.1/1d5) * n(1,patmo_idx_H2O2)
    end if
    if (n(1,patmo_idx_HO2) > 0.1/1d5) then
      dn(1,patmo_idx_HO2) = dn(1,patmo_idx_HO2) - (0.1/1d5) * n(1,patmo_idx_HO2)
    end if
    if (n(1,patmo_idx_CH2O) > 0.1/1d5) then
      dn(1,patmo_idx_CH2O) = dn(1,patmo_idx_CH2O) - (0.1/1d5) * n(1,patmo_idx_CH2O)
    end if
    if (n(1,patmo_idx_CHO) > 0.1/1d5) then
      dn(1,patmo_idx_CHO) = dn(1,patmo_idx_CHO) - (0.1/1d5) * n(1,patmo_idx_CHO)
    end if
    if (n(1,patmo_idx_OH) > 0.1/1d5) then
      dn(1,patmo_idx_OH) = dn(1,patmo_idx_OH) - (0.1/1d5) * n(1,patmo_idx_OH)
    end if
    if (n(1,patmo_idx_O) > 0.1/1d5) then
      dn(1,patmo_idx_O) = dn(1,patmo_idx_O) - (0.1/1d5) * n(1,patmo_idx_O)
    end if
    if (n(1,patmo_idx_H) > 0.1/1d5) then
      dn(1,patmo_idx_H) = dn(1,patmo_idx_H) - (0.1/1d5) * n(1,patmo_idx_H)
    end if
    if (n(1,patmo_idx_CO) > 0.1/1d5) then
      dn(1,patmo_idx_CO) = dn(1,patmo_idx_CO) - (0.1/1d5) * n(1,patmo_idx_CO)
    end if
    if (n(1,patmo_idx_HOCO) > 0.1/1d5) then
      dn(1,patmo_idx_HOCO) = dn(1,patmo_idx_HOCO) - (0.1/1d5) * n(1,patmo_idx_HOCO)
    end if
    if (n(1,patmo_idx_CH3) > 0.1/1d5) then
      dn(1,patmo_idx_CH3) = dn(1,patmo_idx_CH3) - (0.1/1d5) * n(1,patmo_idx_CH3)
    end if

    !Volcanic emission
    !The release of 1 Tmol/year from Claire et al., 2014, with an H2S:SO2 ratio of 1:10
    !The release of molecular hydrogen 3 Tmol/year from Claire et al., 2014
    dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) + 3.85d9/1d5/11d0
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) + 10d0*3.85d9/1d5/11d0
    dn(1,patmo_idx_CO) = dn(1,patmo_idx_CO) + 3d9/1d5
    dn(1,patmo_idx_H2) = dn(1,patmo_idx_H2) + 1d10/1d5
    dn(1,patmo_idx_NH3) = dn(1,patmo_idx_NH3) + 5.3d9/1d5

    ! Water Removal
    dn(:,patmo_idx_H2O) = dn(:,patmo_idx_H2O) - n(:,patmo_idx_H2O) * condenseH2O(:)

    ! Wet Deposition
    do j=12, 2, -1
      do i = 1, chemSpeciesNumber
        dn(j,     i) = dn(j,     i) - wetdep(j, i) * n(j, i)
        dn(j - 1, i) = dn(j - 1, i) + wetdep(j, i) * n(j, i)
      end do
    end do
    do i = 1, chemSpeciesNumber
      dn(1, i) = dn(1, i) - wetdep(1, i) * n(1, i)
    end do

    !aerosol formation
    do i=13,34
      if (va(i) <= n(i, patmo_idx_H2SO4) .AND. pa(i) >= n(i, patmo_idx_H2SO4)) then
        dn(i, patmo_idx_H2SO4) = dn(i, patmo_idx_H2SO4) - (n(i, patmo_idx_H2SO4) - va(i))
        dn(i, patmo_idx_SO4)   = dn(i, patmo_idx_SO4)   + (n(i, patmo_idx_H2SO4) - va(i))
      end if
    end do

    ! Gravity Settling
    do j = cellsNumber, 2, -1
      dn(j    , patmo_idx_SO4) = dn(j    , patmo_idx_SO4) - gd(j) * n(j, patmo_idx_SO4)
      dn(j - 1, patmo_idx_SO4) = dn(j - 1, patmo_idx_SO4) + gd(j) * n(j, patmo_idx_SO4)
    end do
    SO4SurFall = gd(1) * n(1, patmo_idx_SO4)
    dn(1, patmo_idx_SO4) = dn(1, patmo_idx_SO4) - SO4SurFall
    do j = cellsNumber, 2, -1
      dn(j    , patmo_idx_S8) = dn(j    , patmo_idx_S8) - 2.62d-4 * n(j, patmo_idx_S8)
      dn(j - 1, patmo_idx_S8) = dn(j - 1, patmo_idx_S8) + 2.62d-4 * n(j, patmo_idx_S8)
    end do
    S8SurFall = 2.62d-4 * n(1, patmo_idx_S8)
    dn(1, patmo_idx_S8) = dn(1, patmo_idx_S8) - S8SurFall

    ! Hydrogen Escape
    if (n(cellsNumber, patmo_idx_H) > Hesc) then
      dn(cellsNumber, patmo_idx_H) = dn(cellsNumber, patmo_idx_H) - Hesc
      !print *, "triggered H escape"
    else
      n(cellsNumber, patmo_idx_H) = 0d0
      dn(cellsNumber, patmo_idx_H) = 0d0
    end if

    if (n(cellsNumber, patmo_idx_H2) > H2esc) then
      dn(cellsNumber, patmo_idx_H2) = dn(cellsNumber, patmo_idx_H2) - H2esc
      !print *, "triggered H2 escape"
    else
      n(cellsNumber, patmo_idx_H2) = 0d0
      dn(cellsNumber, patmo_idx_H2) = 0d0
    end if

    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode

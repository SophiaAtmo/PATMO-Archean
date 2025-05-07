module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,M(cellsNumber)
    integer::icell
    real*8::n(cellsNumber,speciesNumber)

    call computeTotalDensityExcludingM(n)
    !total density per layer
    M(:) = n(:, patmo_idx_M)

    !loop on cells
    do icell=1,cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0/Tgas
      !O + O2 + M -> O3 + M
      krate(icell,1) = 3.11d-34*(T/298)**(-2.0)

      !O + O3 -> O2 + O2
      krate(icell,2) = 1.83d-11*exp(-2164/T)

      !O(1D) + O3 -> O2 + O2
      krate(icell,3) = 1.2d-10

      !O(1D) + O3 -> O2 + O + O
      krate(icell,4) = 1.2d-10

      !O(1D) + N2 -> O + N2
      krate(icell,5) = 2.15d-11*exp(110/T)

      !O(1D) + O2 -> O + O2
      krate(icell,6) = 3.3d-11*exp(55/T)

      !OH + O3 -> HO2 + O2
      krate(icell,7) = 1.7d-12*exp(-940/T)

      !HO2 + O3 -> OH + O2 + O2
      krate(icell,8) = 1.0d-14*exp(-490/T)

      !OH + HO2 -> H2O + O2
      krate(icell,9) = 4.8d-11*exp(250/T)

      !O(1D) + H2O -> OH + OH
      krate(icell,10) = 1.63d-10*exp(60/T)

      !H2O + H -> OH + H2
      krate(icell,11) = 6.82d-12*(T/298)**(1.6d0)*exp(-9720/T)

      !O(1D) + N2 -> N2O
      krate(icell,12) = 2.8d-36*(T/300)**(-0.9)

      !O(1D) + N2O -> N2 + O2
      krate(icell,13) = 0.39*1.19d-10*exp(20/T)

      !O(1D) + N2O -> NO + NO
      krate(icell,14) = 0.61*1.19d-10*exp(20/T)

      !O + NO2 -> NO + O2
      krate(icell,15) = 5.1d-12*exp(210/T)

      !NO + O3 -> NO2 + O2
      krate(icell,16) = 3.0d-12*exp(-1500/T)

      !NO2 + O3 -> NO3 + O2
      krate(icell,17) = 1.2d-13*exp(-2450/T)

      !NO2 + NO3 + M -> N2O5 + M
      krate(icell,18) = ((2.4d-30*(T/300)**(-3d0))*M(icell)/(1d0+((2.4d-30*(T/300)**(-3d0))*M(icell)/(4.0d-12*(T/300)**(1d-1))))*6d-1**((1d0+(log10((2.4d-30*(T/300)**(-3d0))*M(icell)/(4.0d-12*(T/300)**(1d-1))))**(2d0))**(-1d0)))/M(icell)

      !NO2 + OH + M -> HNO3 + M
      krate(icell,19) = ((1.8d-30*(T/300)**(-3d0))*M(icell)/(1d0+((1.8d-30*(T/300)**(-3d0))*M(icell)/(2.8d-11)))*6d-1**((1d0+(log10((1.8d-30*(T/300)**(-3d0))*M(icell)/(2.8d-11)))**(2d0))**(-1d0)))/M(icell)

      !HNO3 + OH + M -> NO3 + H2O + M
      krate(icell,20) = (2.4d-14*exp(460/T)+(6.5d-34*exp(1335/T)*M(icell))/(1d0+(6.5d-34*exp(1335/T)*M(icell))/2.7d-17*exp(2119/T)))/M(icell)

      !HO2 + NO -> OH + NO2
      krate(icell,21) = 3.3d-12*exp(270/T)

      !H + O3 -> OH + O2
      krate(icell,22) = 1.4d-10*exp(-470/T)

      !O + OH -> H + O2
      krate(icell,23) = 1.8d-11*exp(180/T)

      !H + O2 + M -> HO2 + M
      krate(icell,24) = ((4.4d-32*(T/300)**(-1.3d0))*M(icell)/(1d0+((4.4d-32*(T/300)**(-1.3d0))*M(icell)/(7.5d-11*(T/300)**(2d-1))))*6d-1**((1d0+(log10((4.4d-32*(T/300)**(-1.3d0))*M(icell)/(7.5d-11*(T/300)**(2d-1))))**(2d0))**(-1d0)))/M(icell)

      !O + HO2 -> OH + O2
      krate(icell,25) = 3.0d-11*exp(200/T)

      !H + HO2 -> OH + OH
      krate(icell,26) = 7.2d-11

      !H + HO2 -> O + H2O
      krate(icell,27) = 1.6d-12

      !H + HO2 -> H2 + O2
      krate(icell,28) = 5.0d-12

      !CH4 + OH -> CH3 + H2O
      krate(icell,29) = 2.45d-12*exp(-1775/T)

      !CH3 + O2 + M -> CH3O2 + M
      krate(icell,30) = ((4.1d-31*(T/300)**(-3.6d0))*M(icell)/(1d0+((4.1d-31*(T/300)**(-3.6d0))*M(icell)/(1.2d-12*(T/300)**(1.1d0))))*6d-1**((1d0+(log10((4.1d-31*(T/300)**(-3.6d0))*M(icell)/(1.2d-12*(T/300)**(1.1d0))))**(2d0))**(-1d0)))/M(icell)

      !CH3O2 + HO2 -> CH3OOH + O2
      krate(icell,31) = 4.1d-13*exp(750/T)

      !CH3O2 + NO -> CH3O + NO2
      krate(icell,32) = 2.8d-12*exp(300/T)

      !CH3OOH + OH -> CH2O + OH + H2O
      krate(icell,33) = 3d-1*3.8d-12*exp(200/T)

      !CH3OOH + OH -> CH3O2 + H2O
      krate(icell,34) = 7d-1*3.8d-12*exp(200/T)

      !CH3O + O2 -> CH2O + HO2
      krate(icell,35) = 3.9d-14*exp(-900/T)

      !CH2O + OH -> CHO + H2O
      krate(icell,36) = 5.5d-12*exp(125/T)

      !CHO + O2 -> CO + HO2
      krate(icell,37) = 5.6d-12*(T/298)**(-0.4)

      !CO + OH + M -> CO2 + H + M
      krate(icell,38) = ((5.9d-33*(T/300)**(-1d0)*M(icell)/(1d0+(5.9d-33*(T/300)**(-1d0)*M(icell)/1.1d-12*(T/300)**(1.3d0))))*6d-1**((1d0+(log10(5.9d-33*(T/300)**(-1d0)*M(icell)/1.1d-12*(T/300)**(1.3d0)))**(2d0))**(-1d0))+1.5d-13)/M(icell)

      !HO2 + HO2 + M -> H2O2 + O2 + M
      krate(icell,39) = 2.1d-33*exp(920/T)

      !H2O2 + OH -> HO2 + H2O
      krate(icell,40) = 1.8d-12

      !COS + OH -> CO2 + SH
      krate(icell,41) = 1.1d-13*exp(-1200/T)

      !COS + O -> CO + SO
      krate(icell,42) = 2.1d-11*exp(-2200/T)

      !CS2 + OH -> SH + COS
      krate(icell,43) = 2d-15

      !CS2 + O -> CS + SO
      krate(icell,44) = 3.20d-11*exp(-650/T)

      !CS + O2 -> COS + O
      krate(icell,45) = 2.9d-19

      !CS + O3 -> COS + O2
      krate(icell,46) = 3.0d-16

      !CS + O -> CO + S
      krate(icell,47) = 2.70d-10*exp(-761/T)

      !H2S + OH -> H2O + SH
      krate(icell,48) = 6.10d-12*exp(-75/T)

      !H2S + O -> OH + SH
      krate(icell,49) = 9.22d-12*exp(-1803/T)

      !H2S + H -> H2 + SH
      krate(icell,50) = 8.00d-13

      !H2S + HO2 -> H2O + HSO
      krate(icell,51) = 3.00d-15

      !SH + O -> H + SO
      krate(icell,52) = 1.60d-10

      !SH + O2 -> OH + SO
      krate(icell,53) = 4.00d-19

      !SH + O3 -> HSO + O2
      krate(icell,54) = 9.00d-12*exp(-280/T)

      !SH + NO2 -> HSO + NO
      krate(icell,55) = 3.00d-11*exp(250/T)

      !SO + O3 -> SO2 + O2
      krate(icell,56) = 4.50d-12*exp(-1170/T)

      !SO + O2 -> SO2 + O
      krate(icell,57) = 1.60d-13*exp(-2282/T)

      !SO + OH -> SO2 + H
      krate(icell,58) = 2.70d-11*exp(335/T)

      !SO + NO2 -> SO2 + NO
      krate(icell,59) = 1.40d-11

      !S + O2 -> SO + O
      krate(icell,60) = 2.31d-12

      !S + O3 -> O2 + SO
      krate(icell,61) = 1.20d-11

      !S + OH -> H + SO
      krate(icell,62) = 6.59d-11

      !SO2 + HO2 -> OH + SO3
      krate(icell,63) = 1.00d-18

      !SO2 + NO2 -> SO3 + NO
      krate(icell,64) = 2.32d-26

      !SO2 + O3 -> SO3 + O2
      krate(icell,65) = 3.00d-12*exp(-7000/T)

      !HSO + O2 -> SO2 + OH
      krate(icell,66) = 1.69d-15

      !HSO + O3 -> O2 + O2 + SH
      krate(icell,67) = 2.54d-13*exp(-392.4/T)

      !HSO + NO2 -> NO + HSO2
      krate(icell,68) = 9.60d-12

      !HSO2 + O2 -> HO2 + SO2
      krate(icell,69) = 3.01d-13

      !HSO3 + O2 -> HO2 + SO3
      krate(icell,70) = 1.30d-12*exp(-330/T)

      !SO3 + H2O -> H2SO4
      krate(icell,71) = 1.20d-15

      !SO2 + O + M -> SO3 + M
      krate(icell,72) = ((1.8d-33*(T/298)**(2d0))*M(icell)/(1d0+((1.8d-33*(T/298)**(2d0))*M(icell)/(4.1d-14*(T/298)**(1.8d0))))*6d-1**((1d0+(log10((1.8d-33*(T/298)**(2d0))*M(icell)/(4.1d-14*(T/298)**(1.8d0))))**(2d0))**(-1d0)))/M(icell)

      !SO2 + OH + M -> HSO3 + M
      krate(icell,73) = ((2.9d-31*(T/298)**(-4.1d0))*M(icell)/(1d0+((2.9d-31*(T/298)**(-4.1d0))*M(icell)/(1.7d-12*(T/298)**(0.2d0))))*6d-1**((1d0+(log10((2.9d-31*(T/298)**(-4.1d0))*M(icell)/(1.7d-12*(T/298)**(0.2d0))))**(2d0))**(-1d0)))/M(icell)

      !CH3SCH3 + OH -> SO2
      krate(icell,74) = 1.20d-11*exp(-260/T)

      !CH3SCH3 + O -> SO2
      krate(icell,75) = 1.00d-11*exp(410/T)

      !CH3SCH3 + OH + M -> SO2 + CH4O3S + M
      krate(icell,76) = (3.04d-12*exp(350/T)*5.53d-31*exp(7460/T)*M(icell)/(1d0+5.53d-31*exp(7460/T)*M(icell)))/M(icell)

      !H2SO4 -> SO4
      krate(icell,77) = 0d0

      !O(1D) + CH4 -> CH3 + OH
      krate(icell,78) = 0.75d0*1.75d-10

      !O(1D) + CH4 -> CH3O + H
      krate(icell,79) = 0.2d0*1.75d-10

      !O(1D) + CH4 -> CH2O + H2
      krate(icell,80) = 0.05d0*1.75d-10

      !CH3O2 + CH3O2 -> CH3O + CH3O + O2
      krate(icell,81) = (26.2d0*exp(-1130/T))/(1+26.2d0*exp(-1130/T))*9.5d-14*exp(390/T)

      !CH3O2 + CH3O2 -> CH3OH + CH2O + O2
      krate(icell,82) = (1d0/(1+26.2*exp(-1130/T)))*9.5d-14*exp(390/T)

      !O + CO + M -> CO2 + M
      krate(icell,83) = 1.7d-33*exp(-1510/T)

      !H + CO + M -> CHO + M
      krate(icell,84) = 5.29d-34*exp(-370/T)

      !H + CHO -> H2 + CO
      krate(icell,85) = 1.83d-10

      !CHO + CHO -> CH2O + CO
      krate(icell,86) = 4.48d-11

      !OH + CHO -> H2O + CO
      krate(icell,87) = 1.69d-10

      !O + CHO -> H + CO2
      krate(icell,88) = 5d-11

      !O + CHO -> OH + CO
      krate(icell,89) = 5d-11

      !H + CH2O -> H2 + CHO
      krate(icell,90) = 1.44d-11*exp(-1744/T)

      !O + CH2O -> OH + CHO
      krate(icell,91) = 1.78d-11*(T/298)**(0.57d0)*exp(-1390/T)

      !O(1D) + H2 -> H + OH
      krate(icell,92) = 1.2d-10

      !OH + H2 -> H + H2O
      krate(icell,93) = 7.7d-12*exp(-2100/T)

      !SO + HO2 -> SO2 + OH
      krate(icell,94) = 2.8d-11

      !SO + SO + M -> S2O2 + M
      krate(icell,95) = 4.41d-31

      !SO + S2O2 -> SO2 + S2O
      krate(icell,96) = 3.3d-14

      !SO + SO -> S + SO2
      krate(icell,97) = 3.5d-15

      !SO + SO3 -> SO2 + SO2
      krate(icell,98) = 1.99d-15

      !SH + SH -> S + H2S
      krate(icell,99) = 1.5d-11

      !SH + H -> H2 + S
      krate(icell,100) = 3d-11

      !SH + CH2O -> H2S + CHO
      krate(icell,101) = 1.7d-11*exp(-800/T)

      !S + S + M -> S2 + M
      krate(icell,102) = ((1.98d-33*exp(205/T))*M(icell)/(1d0+((1.98d-33*exp(205/T))*M(icell)/(2.26d-14*exp(414/T))))*6d-1**((1d0+(log10((1.98d-33*exp(205/T))*M(icell)/(2.26d-14*exp(414/T))))**(2d0))**(-1d0)))/M(icell)

      !S + S2 + M -> S3 + M
      krate(icell,103) = 1.0d-25*T**(-2d0)

      !S + S3 + M -> S4 + M
      krate(icell,104) = 1.0d-25*T**(-2d0)

      !S2 + S2 + M -> S4 + M
      krate(icell,105) = 1d-29

      !S4 + S4 + M -> S8 + M
      krate(icell,106) = 4.0d-31*exp(900/T)

      !S2 + M -> S + S + M
      krate(icell,107) = 8.97d-9*exp(-50637/T)

      !S2 + O -> S + SO
      krate(icell,108) = 1.66d-11

      !O + CH3 -> CH3O
      krate(icell,109) = 7.51d-14*(T/298)**(-2.12d0)*exp(-314/T)

      !O + CH3 -> CH2O + H
      krate(icell,110) = 1.4d-10*exp(-470/T)

      !H + CH3 + M -> CH4 + M
      krate(icell,111) = ((6.01d-29*(T/298)**(-1.8d0))*M(icell)/(1d0+((6.01d-29*(T/298)**(-1.8d0))*M(icell)/(3.5d-10)))*6d-1**((1d0+(log10((6.01d-29*(T/298)**(-1.8d0))*M(icell)/(3.5d-10)))**(2d0))**(-1d0)))/M(icell)

      !O3 + CH3 -> CH3O + O2
      krate(icell,112) = 9.68d-14

      !H2O2 + CH3 -> CH4 + HO2
      krate(icell,113) = 2.01d-14*exp(299/T)

      !OH + CH3 -> CH3O + H
      krate(icell,114) = 2.57d-12*(T/298)**(-2.3d-1)*exp(-7010/T)

      !OH + CH3 -> CH4 + O
      krate(icell,115) = 3.22d-14*(T/298)**(2.2d0)*exp(-2240/T)

      !OH + CH3 + M -> CH3OH + M
      krate(icell,116) = ((2.3d-27)*M(icell)/(1d0+((2.3d-27)*M(icell)/(1d-10)))*6d-1**((1d0+(log10((2.3d-27)*M(icell)/(1d-10)))**(2d0))**(-1d0)))/M(icell)

      !HO2 + CH3 -> CH3O + OH
      krate(icell,117) = 3.01d-11

      !HO2 + CH3 -> CH4 + O2
      krate(icell,118) = 5.99d-12

      !CHO + CH3 -> CH4 + CO
      krate(icell,119) = 2.01d-10

      !CH3 + CH3 -> CH4 + CH2
      krate(icell,120) = 7.14d-12*exp(-5051/T)

      !CH3 -> H + CH2
      krate(icell,121) = 1.69d-8*exp(-45586/T)

      !CH3 -> H2 + CH
      krate(icell,122) = 8.3d-9*exp(-42819/T)

      !CH3O + CH3 -> CH2O + CH4
      krate(icell,123) = 4d-11

      !CH2OH + CH3 -> CH2O + CH4
      krate(icell,124) = 4d-12

      !CH3O2 + CH3 -> CH3O + CH3O
      krate(icell,125) = 4d-11

      !H2 + CH3 -> CH4 + H
      krate(icell,126) = 6.86d-14*(T/298)**(2.74d0)*exp(-4740/T)

      !O + CH2 -> CHO + H
      krate(icell,127) = 5.01d-11

      !O + CH2 -> H + H + CO
      krate(icell,128) = 1.33d-10

      !O + CH2 -> H2 + CO
      krate(icell,129) = 6.64d-11

      !H + CH2 -> H2 + CH
      krate(icell,130) = 1d-11*exp(900/T)

      !O2 + CH2 -> H + H + CO2
      krate(icell,131) = 8.3d-14

      !O2 + CH2 -> H2 + CO2
      krate(icell,132) = 2.99d-11*(T/298)**(-3.3d0)*exp(-1440/T)

      !O2 + CH2 -> CO + H2O
      krate(icell,133) = 4d-13

      !O2 + CH2 -> O + CH2O
      krate(icell,134) = 6.64d-14

      !OH + CH2 -> H + CH2O
      krate(icell,135) = 3.01d-11

      !CHO + CH2 -> CO + CH3
      krate(icell,136) = 3.01d-11

      !CH3O2 + CH2 -> CH2O + CH3O
      krate(icell,137) = 3.01d-11

      !CO2 + CH2 -> CH2O + CO
      krate(icell,138) = 3.90d-14

      !O + CH -> H + CO
      krate(icell,139) = 6.59d-11

      !CH + NO2 -> CHO + NO
      krate(icell,140) = 1.45d-10

      !O2 + CH -> O + CHO
      krate(icell,141) = 1.66d-11

      !O2 + CH -> OH + CO
      krate(icell,142) = 8.3d-11

      !H2O + CH -> H + CH2O
      krate(icell,143) = 2.92d-11*(T/298)**(-1.22d0)*exp(-12/T)

      !H2 + CH -> H + CH2
      krate(icell,144) = 3.75d-10*exp(-1660/T)

      !H2 + CH -> CH3
      krate(icell,145) = 2.01d-10*(T/298)**(1.5d-1)

      !CH3OH + CH2 -> CH3O + CH3
      krate(icell,146) = 1.12d-15*(T/298)**(3.1d0)*exp(-3490/T)

      !CH3OH + CH2 -> CH2OH + CH3
      krate(icell,147) = 4.38d-15*(T/298)**(3.2d0)*exp(-3611/T)

      !CH3OH + O -> CH3O + OH
      krate(icell,148) = 1.66d-11*exp(-2360/T)

      !CH3OH + O -> CH2OH + OH
      krate(icell,149) = 1.63d-11*exp(-2270/T)

      !CH3OH + H -> CH3 + H2O
      krate(icell,150) = 1.91d-28

      !CH3OH + H -> CH3O + H2
      krate(icell,151) = 6.64d-11*exp(-3071/T)

      !CH3OH + H -> CH2OH + H2
      krate(icell,152) = 2.42d-12*(T/298)**(2d0)*exp(-2270/T)

      !CH3OH + OH -> CH3O + H2O
      krate(icell,153) = 1.4d-13

      !CH3OH + OH -> CH2OH + H2O
      krate(icell,154) = 3.1d-12*exp(-360/T)

      !CH3OH + OH -> CH2O + H2O + H
      krate(icell,155) = 1.1d-12*(T/298)**(1.44d0)*exp(-57/T)

      !CH3OH + CH3 -> CH4 + CH3O
      krate(icell,156) = 1.12d-15*(T/298)**(3.1d0)*exp(-3490/T)

      !CH3OH + CH3 -> CH4 + CH2OH
      krate(icell,157) = 4.38d-15*(T/298)**(3.2d0)*exp(-3611/T)

      !CH2OH + CH2 -> CH2O + CH3
      krate(icell,158) = 2.01d-12

      !CH2OH + O -> CH2O + OH
      krate(icell,159) = 7.01d-11

      !CH2OH + H -> CH3 + OH
      krate(icell,160) = 1.6d-10

      !CH2OH + H -> CH3OH
      krate(icell,161) = 2.89d-10*(T/298)**(4d-2)

      !CH2OH + H -> CH2O + H2
      krate(icell,162) = 1d-11

      !CH2OH + H2O2 -> CH3OH + HO2
      krate(icell,163) = 5d-15*exp(-1300/T)

      !CH2OH + OH -> CH2O + H2O
      krate(icell,164) = 4d-11

      !CH2OH + HO2 -> CH2O + H2O2
      krate(icell,165) = 2.01d-11

      !CH2OH + CHO -> CH3OH + CO
      krate(icell,166) = 2.01d-10

      !CH2OH + CHO -> CH2O + CH2O
      krate(icell,167) = 3.01d-10

      !CH2OH + CH2OH -> CH2O + CH3OH
      krate(icell,168) = 8d-12

      !N + O2 -> O + NO
      krate(icell,169) = 1.5d-11*exp(-3600/T)

      !N + NO -> N2 + O
      krate(icell,170) = 2.1d-11*exp(100/T)

      !H + NO2 -> NO + OH
      krate(icell,171) = 4d-10*exp(-340/T)

      !O + NO3 -> O2 + NO2
      krate(icell,172) = 1d-11

      !NH2 + NH2 + M -> N2H4 + M
      krate(icell,173) = ((1.96d-29*(T/298)**(-3.9))*M(icell)/(1d0+((1.96d-29*(T/298)**(-3.9))*M(icell)/(1.18d-10*(T/298)**(0.27))))*6d-1**((1d0+(log10((1.96d-29*(T/298)**(-3.9))*M(icell)/(1.18d-10*(T/298)**(0.27))))**(2d0))**(-1d0)))/M(icell)

      !N2H4 + H -> N2H3 + H2
      krate(icell,174) = 1.17d-11*exp(-1261/T)

      !N2H3 + H -> NH2 + NH2
      krate(icell,175) = 2.66d-12

      !NH + NO -> N2 + OH
      krate(icell,176) = 1.44d-11

      !NH + O -> N + OH
      krate(icell,177) = 1.16d-11

      !NH2 + NO -> N2 + H2O
      krate(icell,178) = 1.19d-11*(T/298)**(-1.85)

      !NH2 + O -> NH + OH
      krate(icell,179) = 1.16d-11

      !NH3 + O(1D) -> NH2 + OH
      krate(icell,180) = 2.51d-10

      !NH3 + OH -> NH2 + H2O
      krate(icell,181) = 3.5d-12*exp(-925/T)

      !NH2 + H + M -> NH3 + M
      krate(icell,182) = 3.01d-30

      !NH + NO -> N2O + H
      krate(icell,183) = 3.12d-11

      !NH + O -> NO + H
      krate(icell,184) = 1.16d-11

      !CH3 + H2S -> CH4 + SH
      krate(icell,185) = 2.29d-13*exp(-1101/T)

      !COS + H -> CO + SH
      krate(icell,186) = 9.07d-12*exp(-1940/T)

      !COS + S -> CO + S2
      krate(icell,187) = 1.5d-12*exp(-1261/T)

      !CS + NO2 -> COS + NO
      krate(icell,188) = 7.6d-17

      !CO + SH -> COS + H
      krate(icell,189) = 4.15d-14*exp(-7661/T)

      !CS2 + O -> CO + S2
      krate(icell,190) = 5.81d-14

      !CS2 + O -> COS + S
      krate(icell,191) = 3.65d-12*exp(-701/T)

      !OH + NH2 -> H2O + NH
      krate(icell,192) = 7.69d-13*(T/298)**(1.5d0)*exp(230/T)

      !NH + NH -> NH2 + N
      krate(icell,193) = 3.74d-15*(T/298)**(3.88d0)*exp(-172/T)

      !NH2 + NH -> NH3 + N
      krate(icell,194) = 1.94d-14*(T/298)**(2.46d0)*exp(-54/T)

      !O + N + M -> NO + M
      krate(icell,195) = 5.46d-33*exp(155/T)

      !H + N + M -> NH + M
      krate(icell,196) = 5.02d-32

      !NO2 + N -> N2O + O
      krate(icell,197) = 5.8d-12*exp(220/T)

      !O + O + M -> O2 + M
      krate(icell,198) = 5.21d-35*exp(900/T)

      !OH + CO + M -> HOCO + M
      krate(icell,199) = ((6.90d-33*(298/T)**(2.1d0))*M(icell)/(1d0+((6.90d-33*(298/T)**(2.1d0))*M(icell)/(1.10d-12*(298/T)**(-1.3))))*6d-1**((1d0+(log10((6.90d-33*(298/T)**(2.1d0))*M(icell)/(1.10d-12*(298/T)**(-1.3))))**(2d0))**(-1d0)))/M(icell)

      !HOCO + O(3P) -> CO2 + OH
      krate(icell,200) = 1.44d-11

      !HOCO + OH -> CO2 + H2O
      krate(icell,201) = 1.03d-11

      !HOCO + CH3 -> H2O + CH2CO
      krate(icell,202) = (1.52+T*1.95d-4)*3.24d-11*T**(0.1024)/(2.52+T*1.95d-4)

      !HOCO + CH3 -> CH4 + CO2
      krate(icell,203) = 3.24d-11*T**(0.1024)/(2.52+T*1.95d-4)

      !HOCO + H -> H2O + CO
      krate(icell,204) = 1.07d-10*0.13

      !HOCO + H -> H2 + CO2
      krate(icell,205) = 1.07d-10*0.87

      !HOCO + CO -> COCOOH
      krate(icell,206) = 2.99d-13*exp(-2586/T)

      !OH + OH + M -> H2O2 + M
      krate(icell,207) = ((6.9d-31*(T/298)**(-1d0))*M(icell)/(1d0+((6.9d-31*(T/298)**(-1d0))*M(icell)/(2.6d-11)))*6d-1**((1d0+(log10((6.9d-31*(T/298)**(-1d0))*M(icell)/(2.6d-11)))**(2d0))**(-1d0)))/M(icell)

      !O(1D) + CO2 -> O(3P) + CO2
      krate(icell,208) = 7.5d-11*exp(115/T)

      !O(1D) + N2 -> O(3P) + N2
      krate(icell,209) = 2.15d-11*exp(110/T)

      !O(1D) + SO2 -> O(3P) + SO2
      krate(icell,210) = 0.24*2.2d-10

      !CH4 + CH2 -> CH3 + CH3
      krate(icell,211) = 5.44d-11*(T/298)**(1.63d0)*exp(-3841/T)

      !O + H2 -> OH + H
      krate(icell,212) = 6.34d-12*exp(-4000/T)

      !H + H -> H2
      krate(icell,213) = 1.6d-32* (298/T)**(2.27d0)

      !HOCO + O2 -> HO2 + CO2
      krate(icell,214) = 2.00d-12

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

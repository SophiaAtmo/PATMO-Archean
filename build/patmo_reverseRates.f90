module patmo_reverseRates
contains

  !compute reverse rates using thermochemistry polynomials
  subroutine computeReverseRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(:)
    real*8::Tgas(cellsNumber)
    real*8::lnTgas(cellsNumber)
    real*8::Tgas2(cellsNumber)
    real*8::Tgas3(cellsNumber)
    real*8::Tgas4(cellsNumber)
    real*8::invTgas(cellsNumber)
    real*8::ntot(cellsNumber)
    integer::i

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)

    !extrapolate lower and upper limits
    do i=1,cellsNumber
      Tgas(i) = max(inTgas(i),2d2)
      Tgas(i) = min(Tgas(i),3d3)
    end do

    lnTgas(:) = log(Tgas(:))
    Tgas2(:) = Tgas(:)**2
    Tgas3(:) = Tgas(:)**3
    Tgas4(:) = Tgas(:)**4
    invTgas(:) = 1d0/Tgas(:)

    !O3 + M -> O + O2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,254) = krate(i,1)*exp(3.543341d0*(lnTgas(i)-1d0) &
            - 4.164922d-3*Tgas(i) &
            + 4.402935d-7*Tgas2(i) &
            + 5.434827d-10*Tgas3(i) &
            - 2.202172d-13*Tgas4(i) &
            - 1.219382d4*invTgas(i) &
            - 2.572867d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,254) = krate(i,1)*exp(-6.125694d0*(lnTgas(i)-1d0) &
            + 6.280764d-3*Tgas(i) &
            - 1.355459d-6*Tgas2(i) &
            + 1.4979d-10*Tgas3(i) &
            - 6.392726d-15*Tgas4(i) &
            - 1.533445d4*invTgas(i) &
            + 4.921999d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,254) = 0d0
      end if
    end do

    !O2 + O2 -> O + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,255) = krate(i,2)*exp(-9.892634d-1*(lnTgas(i)-1d0) &
            + 2.38397d-3*Tgas(i) &
            + 1.328442d-7*Tgas2(i) &
            - 7.580525d-10*Tgas3(i) &
            + 2.692968d-13*Tgas4(i) &
            - 4.711464d4*invTgas(i) &
            + 3.019058d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,255) = krate(i,2)*exp(7.552007d0*(lnTgas(i)-1d0) &
            - 6.636263d-3*Tgas(i) &
            + 1.377587d-6*Tgas2(i) &
            - 1.506792d-10*Tgas3(i) &
            + 6.409727d-15*Tgas4(i) &
            - 4.433355d4*invTgas(i) &
            - 4.279077d1)
      else
        krate(i,255) = 0d0
      end if
    end do

    !O2 + O2 -> O(1D) + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,256) = krate(i,3)*exp(-1.657593d0*(lnTgas(i)-1d0) &
            + 4.023715d-3*Tgas(i) &
            - 9.743907d-7*Tgas2(i) &
            - 2.473494d-10*Tgas3(i) &
            + 1.636552d-13*Tgas4(i) &
            - 6.998892d4*invTgas(i) &
            + 5.58397d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,256) = krate(i,3)*exp(7.502055d0*(lnTgas(i)-1d0) &
            - 6.615724d-3*Tgas(i) &
            + 1.376612d-6*Tgas2(i) &
            - 1.508621d-10*Tgas3(i) &
            + 6.423625d-15*Tgas4(i) &
            - 6.710617d4*invTgas(i) &
            - 4.306255d1)
      else
        krate(i,256) = 0d0
      end if
    end do

    !O2 + O + O -> O(1D) + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,257) = krate(i,4)*exp(-4.21167d0*(lnTgas(i)-1d0) &
            + 5.804667d-3*Tgas(i) &
            - 1.547528d-6*Tgas2(i) &
            - 3.277961d-11*Tgas3(i) &
            + 1.145757d-13*Tgas4(i) &
            - 1.068045d4*invTgas(i) &
            + 5.137779d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,257) = krate(i,4)*exp(6.075742d0*(lnTgas(i)-1d0) &
            - 6.260225d-3*Tgas(i) &
            + 1.354484d-6*Tgas2(i) &
            - 1.499729d-10*Tgas3(i) &
            + 6.406623d-15*Tgas4(i) &
            - 7.438167d3*invTgas(i) &
            - 4.949178d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,257) = 0d0
      end if
    end do

    !O + N2 -> O(1D) + N2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,258) = krate(i,5)*exp(-6.683292d-1*(lnTgas(i)-1d0) &
            + 1.639745d-3*Tgas(i) &
            - 1.107235d-6*Tgas2(i) &
            + 5.107031d-10*Tgas3(i) &
            - 1.056415d-13*Tgas4(i) &
            - 2.287427d4*invTgas(i) &
            + 2.564912d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,258) = krate(i,5)*exp(-4.995222d-2*(lnTgas(i)-1d0) &
            + 2.053902d-5*Tgas(i) &
            - 9.749684d-10*Tgas2(i) &
            - 1.828914d-13*Tgas3(i) &
            + 1.389781d-17*Tgas4(i) &
            - 2.277262d4*invTgas(i) &
            - 2.717851d-1)
      else
        krate(i,258) = 0d0
      end if
    end do

    !O + O2 -> O(1D) + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,259) = krate(i,6)*exp(-6.683292d-1*(lnTgas(i)-1d0) &
            + 1.639745d-3*Tgas(i) &
            - 1.107235d-6*Tgas2(i) &
            + 5.107031d-10*Tgas3(i) &
            - 1.056415d-13*Tgas4(i) &
            - 2.287427d4*invTgas(i) &
            + 2.564912d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,259) = krate(i,6)*exp(-4.995222d-2*(lnTgas(i)-1d0) &
            + 2.053902d-5*Tgas(i) &
            - 9.749684d-10*Tgas2(i) &
            - 1.828914d-13*Tgas3(i) &
            + 1.389781d-17*Tgas4(i) &
            - 2.277262d4*invTgas(i) &
            - 2.717851d-1)
      else
        krate(i,259) = 0d0
      end if
    end do

    !HO2 + O2 -> OH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,260) = krate(i,7)*exp(-6.84888d-1*(lnTgas(i)-1d0) &
            + 3.69929d-3*Tgas(i) &
            - 2.090058d-6*Tgas2(i) &
            + 6.456141d-10*Tgas3(i) &
            - 7.060267d-14*Tgas4(i) &
            - 2.003332d4*invTgas(i) &
            + 8.041394d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,260) = krate(i,7)*exp(7.335574d0*(lnTgas(i)-1d0) &
            - 6.681304d-3*Tgas(i) &
            + 1.363473d-6*Tgas2(i) &
            - 1.474934d-10*Tgas3(i) &
            + 6.23879d-15*Tgas4(i) &
            - 1.755835d4*invTgas(i) &
            - 4.141043d1)
      else
        krate(i,260) = 0d0
      end if
    end do

    !OH + O2 + O2 -> HO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,261) = krate(i,8)*exp(-3.847717d0*(lnTgas(i)-1d0) &
            + 2.849602d-3*Tgas(i) &
            + 1.782609d-6*Tgas2(i) &
            - 1.947149d-9*Tgas3(i) &
            + 5.601166d-13*Tgas4(i) &
            - 1.488751d4*invTgas(i) &
            + 4.787785d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,261) = krate(i,8)*exp(6.342127d0*(lnTgas(i)-1d0) &
            - 6.235723d-3*Tgas(i) &
            + 1.369572d-6*Tgas2(i) &
            - 1.529758d-10*Tgas3(i) &
            + 6.563662d-15*Tgas4(i) &
            - 1.144075d4*invTgas(i) &
            - 5.060033d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,261) = 0d0
      end if
    end do

    !H2O + O2 -> OH + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,262) = krate(i,9)*exp(3.126908d-1*(lnTgas(i)-1d0) &
            - 1.058526d-3*Tgas(i) &
            + 1.567881d-6*Tgas2(i) &
            - 1.082194d-9*Tgas3(i) &
            + 2.819875d-13*Tgas4(i) &
            - 3.499059d4*invTgas(i) &
            + 8.03997d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,262) = krate(i,9)*exp(6.728182d-1*(lnTgas(i)-1d0) &
            - 3.204791d-4*Tgas(i) &
            + 4.57735d-8*Tgas2(i) &
            - 4.456473d-12*Tgas3(i) &
            + 1.660746d-16*Tgas4(i) &
            - 3.48307d4*invTgas(i) &
            - 1.49529d0)
      else
        krate(i,262) = 0d0
      end if
    end do

    !OH + OH -> O(1D) + H2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,263) = krate(i,10)*exp(-1.285395d0*(lnTgas(i)-1d0) &
            + 1.382952d-3*Tgas(i) &
            - 4.52214d-7*Tgas2(i) &
            + 1.892309d-10*Tgas3(i) &
            - 4.772965d-14*Tgas4(i) &
            - 1.496501d4*invTgas(i) &
            + 3.975833d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,263) = krate(i,10)*exp(-5.06337d-1*(lnTgas(i)-1d0) &
            + 3.860588d-4*Tgas(i) &
            - 3.263476d-8*Tgas2(i) &
            + 1.087824d-12*Tgas3(i) &
            + 1.876006d-17*Tgas4(i) &
            - 1.471712d4*invTgas(i) &
            - 1.568335d-1)
      else
        krate(i,263) = 0d0
      end if
    end do

    !OH + H2 -> H2O + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,264) = krate(i,11)*exp(3.623198d-1*(lnTgas(i)-1d0) &
            - 3.807928d-3*Tgas(i) &
            + 3.563642d-6*Tgas2(i) &
            - 1.813831d-9*Tgas3(i) &
            + 3.892445d-13*Tgas4(i) &
            + 7.271029d3*invTgas(i) &
            - 1.874704d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,264) = krate(i,11)*exp(-5.943572d-1*(lnTgas(i)-1d0) &
            + 5.195803d-4*Tgas(i) &
            - 5.556105d-8*Tgas2(i) &
            + 3.079467d-12*Tgas3(i) &
            - 5.786477d-17*Tgas4(i) &
            + 7.296976d3*invTgas(i) &
            + 1.615249d0)
      else
        krate(i,264) = 0d0
      end if
    end do

    !N2O -> O(1D) + N2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,265) = krate(i,12)*exp(3.773793d0*(lnTgas(i)-1d0) &
            - 5.714109d-3*Tgas(i) &
            + 2.194662d-6*Tgas2(i) &
            - 6.038586d-10*Tgas3(i) &
            + 7.608674d-14*Tgas4(i) &
            - 4.220778d4*invTgas(i) &
            - 3.173676d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,265) = krate(i,12)*exp(6.231882d-1*(lnTgas(i)-1d0) &
            - 6.081815d-4*Tgas(i) &
            + 7.597284d-8*Tgas2(i) &
            - 6.553832d-12*Tgas3(i) &
            + 2.48304d-16*Tgas4(i) &
            - 4.300128d4*invTgas(i) &
            + 1.272412d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,265) = 0d0
      end if
    end do

    !N2 + O2 -> O(1D) + N2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,266) = krate(i,13)*exp(-2.556374d0*(lnTgas(i)-1d0) &
            + 7.212648d-3*Tgas(i) &
            - 3.835994d-6*Tgas2(i) &
            + 1.410695d-9*Tgas3(i) &
            - 2.382903d-13*Tgas4(i) &
            - 6.284923d4*invTgas(i) &
            + 8.749691d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,266) = krate(i,13)*exp(7.032206d-1*(lnTgas(i)-1d0) &
            + 2.937603d-4*Tgas(i) &
            - 5.57946d-8*Tgas2(i) &
            + 5.298869d-12*Tgas3(i) &
            - 2.035071d-16*Tgas4(i) &
            - 6.221196d4*invTgas(i) &
            - 6.838462d0)
      else
        krate(i,266) = 0d0
      end if
    end do

    !NO + NO -> O(1D) + N2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,267) = krate(i,14)*exp(-3.68011d0*(lnTgas(i)-1d0) &
            + 1.029233d-2*Tgas(i) &
            - 5.960046d-6*Tgas2(i) &
            + 2.363622d-9*Tgas3(i) &
            - 4.270993d-13*Tgas4(i) &
            - 4.104811d4*invTgas(i) &
            + 1.081362d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,267) = krate(i,14)*exp(7.95333d-1*(lnTgas(i)-1d0) &
            + 1.293821d-4*Tgas(i) &
            - 1.838392d-8*Tgas2(i) &
            + 1.989246d-12*Tgas3(i) &
            - 9.55457d-17*Tgas4(i) &
            - 4.022917d4*invTgas(i) &
            - 1.028922d1)
      else
        krate(i,267) = 0d0
      end if
    end do

    !NO + O2 -> O + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,268) = krate(i,15)*exp(-8.88757d-1*(lnTgas(i)-1d0) &
            + 1.385934d-3*Tgas(i) &
            + 4.015448d-7*Tgas2(i) &
            - 6.318035d-10*Tgas3(i) &
            + 1.949219d-13*Tgas4(i) &
            - 2.323772d4*invTgas(i) &
            + 2.42564d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,268) = krate(i,15)*exp(5.06718d-1*(lnTgas(i)-1d0) &
            + 1.48851d-4*Tgas(i) &
            - 4.366452d-8*Tgas2(i) &
            + 6.033498d-12*Tgas3(i) &
            - 2.829179d-16*Tgas4(i) &
            - 2.283706d4*invTgas(i) &
            - 4.97949d0)
      else
        krate(i,268) = 0d0
      end if
    end do

    !NO2 + O2 -> NO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,269) = krate(i,16)*exp(-1.005064d-1*(lnTgas(i)-1d0) &
            + 9.980363d-4*Tgas(i) &
            - 2.687007d-7*Tgas2(i) &
            - 1.26249d-10*Tgas3(i) &
            + 7.437481d-14*Tgas4(i) &
            - 2.387692d4*invTgas(i) &
            + 5.934182d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,269) = krate(i,16)*exp(7.045289d0*(lnTgas(i)-1d0) &
            - 6.785114d-3*Tgas(i) &
            + 1.421251d-6*Tgas2(i) &
            - 1.567127d-10*Tgas3(i) &
            + 6.692645d-15*Tgas4(i) &
            - 2.149649d4*invTgas(i) &
            - 3.781128d1)
      else
        krate(i,269) = 0d0
      end if
    end do

    !NO3 + O2 -> NO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,270) = krate(i,17)*exp(1.395364d0*(lnTgas(i)-1d0) &
            - 3.512586d-3*Tgas(i) &
            + 1.601975d-6*Tgas2(i) &
            - 4.140915d-10*Tgas3(i) &
            + 3.468321d-14*Tgas4(i) &
            - 1.201215d4*invTgas(i) &
            - 3.665417d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,270) = krate(i,17)*exp(6.070608d0*(lnTgas(i)-1d0) &
            - 6.497085d-3*Tgas(i) &
            + 1.384992d-6*Tgas2(i) &
            - 1.506137d-10*Tgas3(i) &
            + 6.378973d-15*Tgas4(i) &
            - 1.007815d4*invTgas(i) &
            - 3.02533d1)
      else
        krate(i,270) = 0d0
      end if
    end do

    !N2O5 + M -> NO2 + NO3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,271) = krate(i,18)*exp(2.42995d0*(lnTgas(i)-1d0) &
            - 1.515362d-2*Tgas(i) &
            + 1.384701d-5*Tgas2(i) &
            - 7.55345d-9*Tgas3(i) &
            + 1.730976d-12*Tgas4(i) &
            - 1.12828d4*invTgas(i) &
            + 8.717415d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,271) = krate(i,18)*exp(-7.425733d-1*(lnTgas(i)-1d0) &
            - 6.212383d-5*Tgas(i) &
            + 6.326582d-9*Tgas2(i) &
            + 1.117869d-12*Tgas3(i) &
            - 8.167839d-17*Tgas4(i) &
            - 1.130573d4*invTgas(i) &
            + 2.040852d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,271) = 0d0
      end if
    end do

    !HNO3 + M -> NO2 + OH + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,272) = krate(i,19)*exp(6.242724d0*(lnTgas(i)-1d0) &
            - 1.150163d-2*Tgas(i) &
            + 4.921037d-6*Tgas2(i) &
            - 1.524454d-9*Tgas3(i) &
            + 2.272941d-13*Tgas4(i) &
            - 2.368541d4*invTgas(i) &
            - 1.097599d1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,272) = krate(i,19)*exp(-3.077051d-1*(lnTgas(i)-1d0) &
            - 5.948887d-4*Tgas(i) &
            + 1.004209d-7*Tgas2(i) &
            - 7.667598d-12*Tgas3(i) &
            + 2.538238d-16*Tgas4(i) &
            - 2.532812d4*invTgas(i) &
            + 2.198918d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,272) = 0d0
      end if
    end do

    !NO3 + H2O + M -> HNO3 + OH + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,273) = krate(i,20)*exp(-6.869527d-1*(lnTgas(i)-1d0) &
            + 4.080918d-3*Tgas(i) &
            - 3.533789d-6*Tgas2(i) &
            + 1.975318d-9*Tgas3(i) &
            - 4.7074d-13*Tgas4(i) &
            - 8.429824d3*invTgas(i) &
            + 3.326785d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,273) = krate(i,20)*exp(7.090038d-1*(lnTgas(i)-1d0) &
            + 1.304827d-5*Tgas(i) &
            - 3.922799d-8*Tgas2(i) &
            + 5.573204d-12*Tgas3(i) &
            - 2.724385d-16*Tgas4(i) &
            - 8.139979d3*invTgas(i) &
            - 3.137444d0)
      else
        krate(i,273) = 0d0
      end if
    end do

    !OH + NO2 -> HO2 + NO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,274) = krate(i,21)*exp(5.843816d-1*(lnTgas(i)-1d0) &
            - 2.701253d-3*Tgas(i) &
            + 1.821357d-6*Tgas2(i) &
            - 7.718631d-10*Tgas3(i) &
            + 1.449775d-13*Tgas4(i) &
            - 3.843602d3*invTgas(i) &
            - 2.107212d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,274) = krate(i,21)*exp(-2.902846d-1*(lnTgas(i)-1d0) &
            - 1.038104d-4*Tgas(i) &
            + 5.777823d-8*Tgas2(i) &
            - 9.219256d-12*Tgas3(i) &
            + 4.538547d-16*Tgas4(i) &
            - 3.938146d3*invTgas(i) &
            + 3.599152d0)
      else
        krate(i,274) = 0d0
      end if
    end do

    !OH + O2 -> H + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,275) = krate(i,22)*exp(-1.867058d0*(lnTgas(i)-1d0) &
            + 3.725796d-3*Tgas(i) &
            - 1.025562d-7*Tgas2(i) &
            - 7.308913d-10*Tgas3(i) &
            + 2.576904d-13*Tgas4(i) &
            - 3.90332d4*invTgas(i) &
            + 4.282116d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,275) = krate(i,22)*exp(8.3308d0*(lnTgas(i)-1d0) &
            - 6.848129d-3*Tgas(i) &
            + 1.40376d-6*Tgas2(i) &
            - 1.528829d-10*Tgas3(i) &
            + 6.489893d-15*Tgas4(i) &
            - 3.566741d4*invTgas(i) &
            - 5.058933d1)
      else
        krate(i,275) = 0d0
      end if
    end do

    !H + O2 -> O + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,276) = krate(i,23)*exp(8.77795d-1*(lnTgas(i)-1d0) &
            - 1.341826d-3*Tgas(i) &
            + 2.354004d-7*Tgas2(i) &
            - 2.716118d-11*Tgas3(i) &
            + 1.160632d-14*Tgas4(i) &
            - 8.081441d3*invTgas(i) &
            - 1.263058d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,276) = krate(i,23)*exp(-7.787933d-1*(lnTgas(i)-1d0) &
            + 2.118654d-4*Tgas(i) &
            - 2.617348d-8*Tgas2(i) &
            + 2.203741d-12*Tgas3(i) &
            - 8.016591d-17*Tgas4(i) &
            - 8.666137d3*invTgas(i) &
            + 7.798561d0)
      else
        krate(i,276) = 0d0
      end if
    end do

    !HO2 + M -> H + O2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,277) = krate(i,24)*exp(1.980658d0*(lnTgas(i)-1d0) &
            + 8.761934d-4*Tgas(i) &
            - 1.885165d-6*Tgas2(i) &
            + 1.216258d-9*Tgas3(i) &
            - 3.024262d-13*Tgas4(i) &
            - 2.41457d4*invTgas(i) &
            - 5.056693d-1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,277) = krate(i,24)*exp(1.988673d0*(lnTgas(i)-1d0) &
            - 6.124052d-4*Tgas(i) &
            + 3.418794d-8*Tgas2(i) &
            + 9.283655d-14*Tgas3(i) &
            - 7.376956d-17*Tgas4(i) &
            - 2.422666d4*invTgas(i) &
            + 1.100322d-2)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,277) = 0d0
      end if
    end do

    !OH + O2 -> O + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,278) = krate(i,25)*exp(-3.043754d-1*(lnTgas(i)-1d0) &
            - 1.31532d-3*Tgas(i) &
            + 2.222902d-6*Tgas2(i) &
            - 1.403667d-9*Tgas3(i) &
            + 3.398994d-13*Tgas4(i) &
            - 2.708132d4*invTgas(i) &
            + 2.214918d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,278) = krate(i,25)*exp(2.164334d-1*(lnTgas(i)-1d0) &
            + 4.504066d-5*Tgas(i) &
            + 1.411371d-8*Tgas2(i) &
            - 3.185758d-12*Tgas3(i) &
            + 1.709368d-16*Tgas4(i) &
            - 2.67752d4*invTgas(i) &
            - 1.380338d0)
      else
        krate(i,278) = 0d0
      end if
    end do

    !OH + OH -> H + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,279) = krate(i,26)*exp(-1.18217d0*(lnTgas(i)-1d0) &
            + 2.650606d-5*Tgas(i) &
            + 1.987502d-6*Tgas2(i) &
            - 1.376505d-9*Tgas3(i) &
            + 3.282931d-13*Tgas4(i) &
            - 1.899988d4*invTgas(i) &
            + 3.477976d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,279) = krate(i,26)*exp(9.952267d-1*(lnTgas(i)-1d0) &
            - 1.668248d-4*Tgas(i) &
            + 4.028719d-8*Tgas2(i) &
            - 5.389499d-12*Tgas3(i) &
            + 2.511027d-16*Tgas4(i) &
            - 1.810906d4*invTgas(i) &
            - 9.178899d0)
      else
        krate(i,279) = 0d0
      end if
    end do

    !O + H2O -> H + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,280) = krate(i,27)*exp(-5.651042d-1*(lnTgas(i)-1d0) &
            + 2.832998d-4*Tgas(i) &
            + 1.332481d-6*Tgas2(i) &
            - 1.055033d-9*Tgas3(i) &
            + 2.703812d-13*Tgas4(i) &
            - 2.690915d4*invTgas(i) &
            + 2.067055d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,280) = krate(i,27)*exp(1.451612d0*(lnTgas(i)-1d0) &
            - 5.323445d-4*Tgas(i) &
            + 7.194698d-8*Tgas2(i) &
            - 6.660215d-12*Tgas3(i) &
            + 2.462405d-16*Tgas4(i) &
            - 2.616456d4*invTgas(i) &
            - 9.293851d0)
      else
        krate(i,280) = 0d0
      end if
    end do

    !H2 + O2 -> H + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,281) = krate(i,28)*exp(6.750106d-1*(lnTgas(i)-1d0) &
            - 4.866454d-3*Tgas(i) &
            + 5.131523d-6*Tgas2(i) &
            - 2.896025d-9*Tgas3(i) &
            + 6.712321d-13*Tgas4(i) &
            - 2.771956d4*invTgas(i) &
            - 1.070707d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,281) = krate(i,28)*exp(7.846101d-2*(lnTgas(i)-1d0) &
            + 1.991012d-4*Tgas(i) &
            - 9.787549d-9*Tgas2(i) &
            - 1.377007d-12*Tgas3(i) &
            + 1.082098d-16*Tgas4(i) &
            - 2.753372d4*invTgas(i) &
            + 1.199597d-1)
      else
        krate(i,281) = 0d0
      end if
    end do

    !CH3 + H2O -> CH4 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,282) = krate(i,29)*exp(1.284427d0*(lnTgas(i)-1d0) &
            - 8.095752d-3*Tgas(i) &
            + 7.002142d-6*Tgas2(i) &
            - 3.414031d-9*Tgas3(i) &
            + 7.067625d-13*Tgas4(i) &
            - 6.994586d3*invTgas(i) &
            - 5.561752d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,282) = krate(i,29)*exp(-9.048432d-1*(lnTgas(i)-1d0) &
            + 9.695294d-4*Tgas(i) &
            - 1.54755d-7*Tgas2(i) &
            + 1.492797d-11*Tgas3(i) &
            - 6.083582d-16*Tgas4(i) &
            - 6.974975d3*invTgas(i) &
            + 2.722335d0)
      else
        krate(i,282) = 0d0
      end if
    end do

    !CH3O2 + M -> CH3 + O2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,283) = krate(i,30)*exp(2.467941d0*(lnTgas(i)-1d0) &
            + 2.211715d-3*Tgas(i) &
            - 5.404621d-6*Tgas2(i) &
            + 3.450597d-9*Tgas3(i) &
            - 8.256281d-13*Tgas4(i) &
            - 1.548779d4*invTgas(i) &
            + 2.516199d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,283) = krate(i,30)*exp(1.083776d0*(lnTgas(i)-1d0) &
            - 1.334072d-3*Tgas(i) &
            + 1.869645d-7*Tgas2(i) &
            - 1.590301d-11*Tgas3(i) &
            + 5.833432d-16*Tgas4(i) &
            - 1.632923d4*invTgas(i) &
            + 1.212943d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,283) = 0d0
      end if
    end do

    !CH3OOH + O2 -> CH3O2 + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,284) = krate(i,31)*exp(2.585628d0*(lnTgas(i)-1d0) &
            - 1.227271d-2*Tgas(i) &
            + 8.960328d-6*Tgas2(i) &
            - 3.918244d-9*Tgas3(i) &
            + 7.41685d-13*Tgas4(i) &
            - 1.80884d4*invTgas(i) &
            - 8.5002d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,284) = krate(i,31)*exp(-1.698749d0*(lnTgas(i)-1d0) &
            + 8.660874d-4*Tgas(i) &
            - 7.726249d-8*Tgas2(i) &
            + 4.080141d-12*Tgas3(i) &
            - 9.370492d-17*Tgas4(i) &
            - 1.85093d4*invTgas(i) &
            + 9.949995d0)
      else
        krate(i,284) = 0d0
      end if
    end do

    !CH3O + NO2 -> CH3O2 + NO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,285) = krate(i,32)*exp(1.534458d0*(lnTgas(i)-1d0) &
            - 2.771692d-3*Tgas(i) &
            + 7.441351d-7*Tgas2(i) &
            + 6.129318d-11*Tgas3(i) &
            - 7.331764d-14*Tgas4(i) &
            - 5.511735d3*invTgas(i) &
            - 7.788779d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,285) = krate(i,32)*exp(-8.265292d-1*(lnTgas(i)-1d0) &
            + 3.497762d-4*Tgas(i) &
            - 2.375307d-8*Tgas2(i) &
            - 6.169716d-13*Tgas3(i) &
            + 9.741164d-17*Tgas4(i) &
            - 6.1791d3*invTgas(i) &
            + 4.461637d0)
      else
        krate(i,285) = 0d0
      end if
    end do

    !CH2O + OH + H2O -> CH3OOH + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,286) = krate(i,33)*exp(-6.086949d0*(lnTgas(i)-1d0) &
            + 1.47221d-2*Tgas(i) &
            - 6.426651d-6*Tgas2(i) &
            + 1.511144d-9*Tgas3(i) &
            - 7.562055d-14*Tgas4(i) &
            - 2.778346d4*invTgas(i) &
            + 1.162041d1)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,286) = krate(i,33)*exp(1.918815d0*(lnTgas(i)-1d0) &
            - 2.75695d-4*Tgas(i) &
            + 7.3772d-9*Tgas2(i) &
            + 6.857414d-13*Tgas3(i) &
            - 6.249903d-17*Tgas4(i) &
            - 2.613658d4*invTgas(i) &
            - 2.73239d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,286) = 0d0
      end if
    end do

    !CH3O2 + H2O -> CH3OOH + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,287) = krate(i,34)*exp(-2.272937d0*(lnTgas(i)-1d0) &
            + 1.121419d-2*Tgas(i) &
            - 7.392447d-6*Tgas2(i) &
            + 2.83605d-9*Tgas3(i) &
            - 4.596974d-13*Tgas4(i) &
            - 1.690218d4*invTgas(i) &
            + 9.304197d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,287) = krate(i,34)*exp(2.371567d0*(lnTgas(i)-1d0) &
            - 1.186566d-3*Tgas(i) &
            + 1.23036d-7*Tgas2(i) &
            - 8.536614d-12*Tgas3(i) &
            + 2.597795d-16*Tgas4(i) &
            - 1.63214d4*invTgas(i) &
            - 1.144528d1)
      else
        krate(i,287) = 0d0
      end if
    end do

    !CH2O + HO2 -> CH3O + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,288) = krate(i,35)*exp(-1.60126d0*(lnTgas(i)-1d0) &
            + 4.428043d-3*Tgas(i) &
            - 1.829648d-6*Tgas2(i) &
            + 4.347009d-10*Tgas3(i) &
            - 2.83473d-14*Tgas4(i) &
            - 1.435896d4*invTgas(i) &
            + 5.910624d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,288) = krate(i,35)*exp(1.076939d0*(lnTgas(i)-1d0) &
            + 1.170434d-5*Tgas(i) &
            - 4.022674d-8*Tgas2(i) &
            + 6.102407d-12*Tgas3(i) &
            - 2.907078d-16*Tgas4(i) &
            - 1.369182d4*invTgas(i) &
            - 7.551193d0)
      else
        krate(i,288) = 0d0
      end if
    end do

    !CHO + H2O -> CH2O + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,289) = krate(i,36)*exp(3.495261d-1*(lnTgas(i)-1d0) &
            - 3.476123d-3*Tgas(i) &
            + 3.569212d-6*Tgas2(i) &
            - 1.90798d-9*Tgas3(i) &
            + 4.197165d-13*Tgas4(i) &
            - 1.541102d4*invTgas(i) &
            - 1.96054d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,289) = krate(i,36)*exp(-5.889973d-1*(lnTgas(i)-1d0) &
            + 9.023218d-4*Tgas(i) &
            - 1.832985d-7*Tgas2(i) &
            + 1.733301d-11*Tgas3(i) &
            - 6.365431d-16*Tgas4(i) &
            - 1.538159d4*invTgas(i) &
            + 1.423705d0)
      else
        krate(i,289) = 0d0
      end if
    end do

    !CO + HO2 -> CHO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,290) = krate(i,37)*exp(1.386709d-1*(lnTgas(i)-1d0) &
            - 4.79006d-4*Tgas(i) &
            + 2.792039d-7*Tgas2(i) &
            + 2.200758d-11*Tgas3(i) &
            - 3.849687d-14*Tgas4(i) &
            - 1.688854d4*invTgas(i) &
            - 2.590471d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,290) = krate(i,37)*exp(3.602028d-1*(lnTgas(i)-1d0) &
            - 2.687266d-5*Tgas(i) &
            + 3.319591d-9*Tgas2(i) &
            + 2.323029d-12*Tgas3(i) &
            - 2.10765d-16*Tgas4(i) &
            - 1.667255d4*invTgas(i) &
            - 1.978641d0)
      else
        krate(i,290) = 0d0
      end if
    end do

    !CO2 + H + M -> CO + OH + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,291) = krate(i,38)*exp(2.714705d0*(lnTgas(i)-1d0) &
            - 5.997775d-3*Tgas(i) &
            + 2.12592d-6*Tgas2(i) &
            - 4.524548d-10*Tgas3(i) &
            + 3.00828d-14*Tgas4(i) &
            - 1.192312d4*invTgas(i) &
            - 6.04981d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,291) = krate(i,38)*exp(-1.249495d0*(lnTgas(i)-1d0) &
            - 1.41158d-4*Tgas(i) &
            + 3.601722d-8*Tgas2(i) &
            - 3.288595d-12*Tgas3(i) &
            + 1.020506d-16*Tgas4(i) &
            - 1.298294d4*invTgas(i) &
            + 1.424362d1)
      else
        krate(i,291) = 0d0
      end if
    end do

    !H2O2 + O2 + M -> HO2 + HO2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,292) = krate(i,39)*exp(5.059883d-1*(lnTgas(i)-1d0) &
            - 2.827059d-3*Tgas(i) &
            + 2.471474d-6*Tgas2(i) &
            - 1.349599d-9*Tgas3(i) &
            + 3.125637d-13*Tgas4(i) &
            - 1.929872d4*invTgas(i) &
            + 5.019155d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,292) = krate(i,39)*exp(1.038411d-1*(lnTgas(i)-1d0) &
            - 4.736367d-4*Tgas(i) &
            + 1.245071d-7*Tgas2(i) &
            - 1.498831d-11*Tgas3(i) &
            + 6.524264d-16*Tgas4(i) &
            - 1.92852d4*invTgas(i) &
            + 1.83502d0)
      else
        krate(i,292) = 0d0
      end if
    end do

    !HO2 + H2O -> H2O2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,293) = krate(i,40)*exp(-1.932975d-1*(lnTgas(i)-1d0) &
            + 1.768533d-3*Tgas(i) &
            - 9.035932d-7*Tgas2(i) &
            + 2.674051d-10*Tgas3(i) &
            - 3.057618d-14*Tgas4(i) &
            - 1.569186d4*invTgas(i) &
            + 3.020815d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,293) = krate(i,40)*exp(5.689771d-1*(lnTgas(i)-1d0) &
            + 1.531575d-4*Tgas(i) &
            - 7.873356d-8*Tgas2(i) &
            + 1.053183d-11*Tgas3(i) &
            - 4.863518d-16*Tgas4(i) &
            - 1.55455d4*invTgas(i) &
            - 3.33031d0)
      else
        krate(i,293) = 0d0
      end if
    end do

    !CO2 + SH -> COS + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,294) = krate(i,41)*exp(-2.775076d-1*(lnTgas(i)-1d0) &
            + 1.258706d-3*Tgas(i) &
            - 4.510004d-7*Tgas2(i) &
            - 6.102868d-11*Tgas3(i) &
            + 6.191498d-14*Tgas4(i) &
            - 1.770436d4*invTgas(i) &
            + 1.658291d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,294) = krate(i,41)*exp(5.450499d-1*(lnTgas(i)-1d0) &
            - 3.939921d-4*Tgas(i) &
            + 5.516733d-8*Tgas2(i) &
            - 4.211364d-12*Tgas3(i) &
            + 1.161815d-16*Tgas4(i) &
            - 1.759891d4*invTgas(i) &
            - 2.155117d0)
      else
        krate(i,294) = 0d0
      end if
    end do

    !CO + SO -> COS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,295) = krate(i,42)*exp(-2.257872d0*(lnTgas(i)-1d0) &
            + 8.400735d-3*Tgas(i) &
            - 5.554705d-6*Tgas2(i) &
            + 2.47746d-9*Tgas3(i) &
            - 4.967152d-13*Tgas4(i) &
            - 2.581411d4*invTgas(i) &
            + 5.859493d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,295) = krate(i,42)*exp(9.007697d-1*(lnTgas(i)-1d0) &
            + 1.738856d-4*Tgas(i) &
            - 5.041413d-8*Tgas2(i) &
            + 5.80007d-12*Tgas3(i) &
            - 2.538475d-16*Tgas4(i) &
            - 2.530287d4*invTgas(i) &
            - 8.614472d0)
      else
        krate(i,295) = 0d0
      end if
    end do

    !SH + COS -> CS2 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,296) = krate(i,43)*exp(7.076339d-1*(lnTgas(i)-1d0) &
            - 2.334753d-3*Tgas(i) &
            + 2.330059d-6*Tgas2(i) &
            - 1.405889d-9*Tgas3(i) &
            + 3.427434d-13*Tgas4(i) &
            - 1.840448d4*invTgas(i) &
            - 3.82013d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,296) = krate(i,43)*exp(3.814896d-1*(lnTgas(i)-1d0) &
            - 2.809352d-4*Tgas(i) &
            + 3.560127d-8*Tgas2(i) &
            - 2.764971d-12*Tgas3(i) &
            + 1.318082d-16*Tgas4(i) &
            - 1.84268d4*invTgas(i) &
            - 2.690905d0)
      else
        krate(i,296) = 0d0
      end if
    end do

    !CS + SO -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,297) = krate(i,44)*exp(-2.009268d0*(lnTgas(i)-1d0) &
            + 1.01334d-2*Tgas(i) &
            - 8.049007d-6*Tgas2(i) &
            + 4.063298d-9*Tgas3(i) &
            - 8.878359d-13*Tgas4(i) &
            - 9.967159d3*invTgas(i) &
            + 3.121036d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,297) = krate(i,44)*exp(7.541485d-1*(lnTgas(i)-1d0) &
            + 2.786439d-4*Tgas(i) &
            - 7.387907d-8*Tgas2(i) &
            + 8.863557d-12*Tgas3(i) &
            - 4.086194d-16*Tgas4(i) &
            - 9.721839d3*invTgas(i) &
            - 8.403234d0)
      else
        krate(i,297) = 0d0
      end if
    end do

    !COS + O -> CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,298) = krate(i,45)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,298) = krate(i,45)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,298) = 0d0
      end if
    end do

    !COS + O2 -> CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,299) = krate(i,46)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,299) = krate(i,46)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,299) = 0d0
      end if
    end do

    !CO + S -> CS + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,300) = krate(i,47)*exp(1.002725d0*(lnTgas(i)-1d0) &
            - 5.273593d-3*Tgas(i) &
            + 5.386224d-6*Tgas2(i) &
            - 3.07128d-9*Tgas3(i) &
            + 7.158338d-13*Tgas4(i) &
            - 4.340154d4*invTgas(i) &
            - 2.970349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,300) = krate(i,47)*exp(3.853828d-1*(lnTgas(i)-1d0) &
            - 6.850666d-5*Tgas(i) &
            - 2.520613d-9*Tgas2(i) &
            - 4.551522d-14*Tgas3(i) &
            + 7.495758d-17*Tgas4(i) &
            - 4.324062d4*invTgas(i) &
            - 1.65598d0)
      else
        krate(i,300) = 0d0
      end if
    end do

    !H2O + SH -> H2S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,301) = krate(i,48)*exp(2.289248d-1*(lnTgas(i)-1d0) &
            - 2.744914d-3*Tgas(i) &
            + 3.195679d-6*Tgas2(i) &
            - 1.867226d-9*Tgas3(i) &
            + 4.404765d-13*Tgas4(i) &
            - 1.407682d4*invTgas(i) &
            + 2.589349d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,301) = krate(i,48)*exp(1.087555d-1*(lnTgas(i)-1d0) &
            + 2.368912d-4*Tgas(i) &
            - 5.712312d-8*Tgas2(i) &
            + 6.875396d-12*Tgas3(i) &
            - 3.158447d-16*Tgas4(i) &
            - 1.386166d4*invTgas(i) &
            - 4.086126d-1)
      else
        krate(i,301) = 0d0
      end if
    end do

    !OH + SH -> H2S + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,302) = krate(i,49)*exp(-3.881414d-1*(lnTgas(i)-1d0) &
            - 3.001707d-3*Tgas(i) &
            + 3.8507d-6*Tgas2(i) &
            - 2.188698d-9*Tgas3(i) &
            + 4.983884d-13*Tgas4(i) &
            - 6.167561d3*invTgas(i) &
            + 1.669856d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,302) = krate(i,49)*exp(-3.476293d-1*(lnTgas(i)-1d0) &
            + 6.02411d-4*Tgas(i) &
            - 8.878291d-8*Tgas2(i) &
            + 8.146112d-12*Tgas3(i) &
            - 3.109825d-16*Tgas4(i) &
            - 5.80616d3*invTgas(i) &
            - 2.936611d-1)
      else
        krate(i,302) = 0d0
      end if
    end do

    !H2 + SH -> H2S + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,303) = krate(i,50)*exp(5.912447d-1*(lnTgas(i)-1d0) &
            - 6.552842d-3*Tgas(i) &
            + 6.759321d-6*Tgas2(i) &
            - 3.681057d-9*Tgas3(i) &
            + 8.29721d-13*Tgas4(i) &
            - 6.805796d3*invTgas(i) &
            - 1.615769d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,303) = krate(i,50)*exp(-4.856017d-1*(lnTgas(i)-1d0) &
            + 7.564716d-4*Tgas(i) &
            - 1.126842d-7*Tgas2(i) &
            + 9.954863d-12*Tgas3(i) &
            - 3.737095d-16*Tgas4(i) &
            - 6.564682d3*invTgas(i) &
            + 1.206637d0)
      else
        krate(i,303) = 0d0
      end if
    end do

    !H2O + HSO -> H2S + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,304) = krate(i,51)*exp(8.775649d-2*(lnTgas(i)-1d0) &
            - 4.496811d-4*Tgas(i) &
            + 3.892062d-7*Tgas2(i) &
            - 1.499713d-10*Tgas3(i) &
            + 2.458951d-14*Tgas4(i) &
            - 3.069931d4*invTgas(i) &
            + 2.097107d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,304) = krate(i,51)*exp(1.268016d-1*(lnTgas(i)-1d0) &
            - 1.406198d-5*Tgas(i) &
            + 2.514851d-8*Tgas2(i) &
            - 3.019166d-12*Tgas3(i) &
            + 1.213293d-16*Tgas4(i) &
            - 3.060978d4*invTgas(i) &
            - 3.04536d-1)
      else
        krate(i,304) = 0d0
      end if
    end do

    !H + SO -> SH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,305) = krate(i,52)*exp(7.343407d-1*(lnTgas(i)-1d0) &
            + 1.144254d-3*Tgas(i) &
            - 2.977785d-6*Tgas2(i) &
            + 2.086034d-9*Tgas3(i) &
            - 5.285474d-13*Tgas4(i) &
            - 2.003287d4*invTgas(i) &
            - 1.848608d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,305) = krate(i,52)*exp(-8.93775d-1*(lnTgas(i)-1d0) &
            + 4.267197d-4*Tgas(i) &
            - 6.956424d-8*Tgas2(i) &
            + 6.722839d-12*Tgas3(i) &
            - 2.679784d-16*Tgas4(i) &
            - 2.068689d4*invTgas(i) &
            + 7.784268d0)
      else
        krate(i,305) = 0d0
      end if
    end do

    !OH + SO -> SH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,306) = krate(i,53)*exp(-1.434543d-1*(lnTgas(i)-1d0) &
            + 2.486079d-3*Tgas(i) &
            - 3.213186d-6*Tgas2(i) &
            + 2.113195d-9*Tgas3(i) &
            - 5.401537d-13*Tgas4(i) &
            - 1.195143d4*invTgas(i) &
            - 5.855506d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,306) = krate(i,53)*exp(-1.149817d-1*(lnTgas(i)-1d0) &
            + 2.148543d-4*Tgas(i) &
            - 4.339076d-8*Tgas2(i) &
            + 4.519098d-12*Tgas3(i) &
            - 1.878125d-16*Tgas4(i) &
            - 1.202075d4*invTgas(i) &
            - 1.42932d-2)
      else
        krate(i,306) = 0d0
      end if
    end do

    !HSO + O2 -> SH + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,307) = krate(i,54)*exp(-8.260563d-1*(lnTgas(i)-1d0) &
            + 5.994522d-3*Tgas(i) &
            - 4.896531d-6*Tgas2(i) &
            + 2.362869d-9*Tgas3(i) &
            - 4.864896d-13*Tgas4(i) &
            - 3.665581d4*invTgas(i) &
            + 7.549152d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,307) = krate(i,54)*exp(7.35362d0*(lnTgas(i)-1d0) &
            - 6.932257d-3*Tgas(i) &
            + 1.445745d-6*Tgas2(i) &
            - 1.57388d-10*Tgas3(i) &
            + 6.675964d-15*Tgas4(i) &
            - 3.430647d4*invTgas(i) &
            - 4.130635d1)
      else
        krate(i,307) = 0d0
      end if
    end do

    !HSO + NO -> SH + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,308) = krate(i,55)*exp(-7.255499d-1*(lnTgas(i)-1d0) &
            + 4.996486d-3*Tgas(i) &
            - 4.62783d-6*Tgas2(i) &
            + 2.489118d-9*Tgas3(i) &
            - 5.608644d-13*Tgas4(i) &
            - 1.277888d4*invTgas(i) &
            + 1.61497d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,308) = krate(i,55)*exp(3.083307d-1*(lnTgas(i)-1d0) &
            - 1.471428d-4*Tgas(i) &
            + 2.44934d-8*Tgas2(i) &
            - 6.753068d-13*Tgas3(i) &
            - 1.668065d-17*Tgas4(i) &
            - 1.280997d4*invTgas(i) &
            - 3.495076d0)
      else
        krate(i,308) = 0d0
      end if
    end do

    !SO2 + O2 -> SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,309) = krate(i,56)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,309) = krate(i,56)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,309) = 0d0
      end if
    end do

    !SO2 + O -> SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,310) = krate(i,57)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,310) = krate(i,57)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,310) = 0d0
      end if
    end do

    !SO2 + H -> SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,311) = krate(i,58)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,311) = krate(i,58)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,311) = 0d0
      end if
    end do

    !SO2 + NO -> SO + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,312) = krate(i,59)*exp(-3.307801d-1*(lnTgas(i)-1d0) &
            - 7.751533d-4*Tgas(i) &
            + 1.465141d-6*Tgas2(i) &
            - 9.739932d-10*Tgas3(i) &
            + 2.437221d-13*Tgas4(i) &
            - 2.95164d4*invTgas(i) &
            + 2.427759d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,312) = krate(i,59)*exp(2.087491d-1*(lnTgas(i)-1d0) &
            - 1.603123d-4*Tgas(i) &
            + 4.013125d-8*Tgas2(i) &
            - 2.849403d-12*Tgas3(i) &
            + 7.811797d-17*Tgas4(i) &
            - 2.92732d4*invTgas(i) &
            - 9.201861d-1)
      else
        krate(i,312) = 0d0
      end if
    end do

    !SO + O -> S + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,313) = krate(i,60)*exp(-6.871497d-1*(lnTgas(i)-1d0) &
            + 3.692253d-3*Tgas(i) &
            - 3.775049d-6*Tgas2(i) &
            + 2.192748d-9*Tgas3(i) &
            - 5.221234d-13*Tgas4(i) &
            - 2.801316d3*invTgas(i) &
            + 1.303125d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,313) = krate(i,60)*exp(2.774641d-2*(lnTgas(i)-1d0) &
            - 1.023326d-4*Tgas(i) &
            + 1.819606d-8*Tgas2(i) &
            - 1.263844d-12*Tgas3(i) &
            + 2.380994d-17*Tgas4(i) &
            - 2.787962d3*invTgas(i) &
            - 1.260456d0)
      else
        krate(i,313) = 0d0
      end if
    end do

    !O2 + SO -> S + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,314) = krate(i,61)*exp(-1.676413d0*(lnTgas(i)-1d0) &
            + 6.076223d-3*Tgas(i) &
            - 3.642205d-6*Tgas2(i) &
            + 1.434695d-9*Tgas3(i) &
            - 2.528266d-13*Tgas4(i) &
            - 4.991596d4*invTgas(i) &
            + 4.322183d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,314) = krate(i,61)*exp(7.579753d0*(lnTgas(i)-1d0) &
            - 6.738596d-3*Tgas(i) &
            + 1.395783d-6*Tgas2(i) &
            - 1.51943d-10*Tgas3(i) &
            + 6.433537d-15*Tgas4(i) &
            - 4.712151d4*invTgas(i) &
            - 4.405122d1)
      else
        krate(i,314) = 0d0
      end if
    end do

    !H + SO -> S + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,315) = krate(i,62)*exp(1.906453d-1*(lnTgas(i)-1d0) &
            + 2.350427d-3*Tgas(i) &
            - 3.539649d-6*Tgas2(i) &
            + 2.165587d-9*Tgas3(i) &
            - 5.105171d-13*Tgas4(i) &
            - 1.088276d4*invTgas(i) &
            + 4.006756d-2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,315) = krate(i,62)*exp(-7.510469d-1*(lnTgas(i)-1d0) &
            + 1.095328d-4*Tgas(i) &
            - 7.977419d-9*Tgas2(i) &
            + 9.398975d-13*Tgas3(i) &
            - 5.635597d-17*Tgas4(i) &
            - 1.14541d4*invTgas(i) &
            + 6.538105d0)
      else
        krate(i,315) = 0d0
      end if
    end do

    !OH + SO3 -> SO2 + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,316) = krate(i,63)*exp(1.61001d0*(lnTgas(i)-1d0) &
            - 8.009682d-3*Tgas(i) &
            + 6.273806d-6*Tgas2(i) &
            - 3.072642d-9*Tgas3(i) &
            + 6.534202d-13*Tgas4(i) &
            - 8.876536d3*invTgas(i) &
            - 1.31498d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,316) = krate(i,63)*exp(-5.787838d-1*(lnTgas(i)-1d0) &
            - 1.413477d-4*Tgas(i) &
            + 6.323959d-8*Tgas2(i) &
            - 7.909567d-12*Tgas3(i) &
            + 3.569003d-16*Tgas4(i) &
            - 9.036184d3*invTgas(i) &
            + 7.706091d0)
      else
        krate(i,316) = 0d0
      end if
    end do

    !SO3 + NO -> SO2 + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,317) = krate(i,64)*exp(1.025629d0*(lnTgas(i)-1d0) &
            - 5.308428d-3*Tgas(i) &
            + 4.452449d-6*Tgas2(i) &
            - 2.300779d-9*Tgas3(i) &
            + 5.084427d-13*Tgas4(i) &
            - 5.032934d3*invTgas(i) &
            - 1.104258d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,317) = krate(i,64)*exp(-2.884992d-1*(lnTgas(i)-1d0) &
            - 3.753731d-5*Tgas(i) &
            + 5.46136d-9*Tgas2(i) &
            + 1.309689d-12*Tgas3(i) &
            - 9.695444d-17*Tgas4(i) &
            - 5.098039d3*invTgas(i) &
            + 4.106939d0)
      else
        krate(i,317) = 0d0
      end if
    end do

    !SO3 + O2 -> SO2 + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,318) = krate(i,65)*exp(9.251222d-1*(lnTgas(i)-1d0) &
            - 4.310392d-3*Tgas(i) &
            + 4.183748d-6*Tgas2(i) &
            - 2.427028d-9*Tgas3(i) &
            + 5.828176d-13*Tgas4(i) &
            - 2.890986d4*invTgas(i) &
            - 5.108402d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,318) = krate(i,65)*exp(6.75679d0*(lnTgas(i)-1d0) &
            - 6.822651d-3*Tgas(i) &
            + 1.426713d-6*Tgas2(i) &
            - 1.55403d-10*Tgas3(i) &
            + 6.59569d-15*Tgas4(i) &
            - 2.659453d4*invTgas(i) &
            - 3.370434d1)
      else
        krate(i,318) = 0d0
      end if
    end do

    !SO2 + OH -> HSO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,319) = krate(i,66)*exp(2.513155d-1*(lnTgas(i)-1d0) &
            - 3.28556d-3*Tgas(i) &
            + 2.879785d-6*Tgas2(i) &
            - 1.349916d-9*Tgas3(i) &
            + 2.644329d-13*Tgas4(i) &
            - 2.868895d4*invTgas(i) &
            + 1.680711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,319) = krate(i,66)*exp(-2.145633d-1*(lnTgas(i)-1d0) &
            + 2.016848d-4*Tgas(i) &
            - 2.77529d-8*Tgas2(i) &
            + 2.345002d-12*Tgas3(i) &
            - 9.301392d-17*Tgas4(i) &
            - 2.848398d4*invTgas(i) &
            + 2.560596d0)
      else
        krate(i,319) = 0d0
      end if
    end do

    !O2 + O2 + SH -> HSO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,320) = krate(i,67)*exp(-3.706548d0*(lnTgas(i)-1d0) &
            + 5.543697d-4*Tgas(i) &
            + 4.589081d-6*Tgas2(i) &
            - 3.664404d-9*Tgas3(i) &
            + 9.760036d-13*Tgas4(i) &
            + 1.734979d3*invTgas(i) &
            + 4.837009d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,320) = krate(i,67)*exp(6.324081d0*(lnTgas(i)-1d0) &
            - 5.98477d-3*Tgas(i) &
            + 1.287301d-6*Tgas2(i) &
            - 1.430812d-10*Tgas3(i) &
            + 6.126488d-15*Tgas4(i) &
            + 5.307368d3*invTgas(i) &
            - 5.070441d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,320) = 0d0
      end if
    end do

    !NO + HSO2 -> HSO + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,321) = krate(i,68)*exp(1.237552d-1*(lnTgas(i)-1d0) &
            - 4.973881d-3*Tgas(i) &
            + 4.997421d-6*Tgas2(i) &
            - 2.592509d-9*Tgas3(i) &
            + 5.552013d-13*Tgas4(i) &
            - 1.905388d4*invTgas(i) &
            - 1.89003d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,321) = krate(i,68)*exp(-4.243896d-1*(lnTgas(i)-1d0) &
            + 4.753954d-4*Tgas(i) &
            - 7.880742d-8*Tgas2(i) &
            + 9.082322d-12*Tgas3(i) &
            - 4.052156d-16*Tgas4(i) &
            - 1.872285d4*invTgas(i) &
            + 2.992486d-1)
      else
        krate(i,321) = 0d0
      end if
    end do

    !HO2 + SO2 -> HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,322) = krate(i,69)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,322) = krate(i,69)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,322) = 0d0
      end if
    end do

    !HO2 + SO3 -> HSO3 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,323) = krate(i,70)*exp(2.907618d-1*(lnTgas(i)-1d0) &
            + 2.512911d-3*Tgas(i) &
            - 2.879387d-6*Tgas2(i) &
            + 1.644321d-9*Tgas3(i) &
            - 3.807075d-13*Tgas4(i) &
            - 3.383301d3*invTgas(i) &
            - 1.13614d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,323) = krate(i,70)*exp(3.861018d-1*(lnTgas(i)-1d0) &
            - 9.114733d-5*Tgas(i) &
            - 1.33545d-8*Tgas2(i) &
            + 3.127919d-12*Tgas3(i) &
            - 1.627562d-16*Tgas4(i) &
            - 3.561354d3*invTgas(i) &
            + 4.998669d-1)
      else
        krate(i,323) = 0d0
      end if
    end do

    !H2SO4 -> SO3 + H2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,324) = krate(i,71)*exp(2.039365d0*(lnTgas(i)-1d0) &
            - 8.55842d-3*Tgas(i) &
            + 5.821711d-6*Tgas2(i) &
            - 2.687074d-9*Tgas3(i) &
            + 5.604965d-13*Tgas4(i) &
            - 1.132526d4*invTgas(i) &
            + 8.315681d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,324) = krate(i,71)*exp(-1.361725d0*(lnTgas(i)-1d0) &
            + 5.032744d-5*Tgas(i) &
            + 1.803258d-8*Tgas2(i) &
            - 2.660542d-12*Tgas3(i) &
            + 1.317101d-16*Tgas4(i) &
            - 1.191318d4*invTgas(i) &
            + 2.406728d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,324) = 0d0
      end if
    end do

    !SO3 + M -> SO2 + O + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,325) = krate(i,72)*exp(4.468463d0*(lnTgas(i)-1d0) &
            - 8.475314d-3*Tgas(i) &
            + 4.624041d-6*Tgas2(i) &
            - 1.883545d-9*Tgas3(i) &
            + 3.626003d-13*Tgas4(i) &
            - 4.110368d4*invTgas(i) &
            - 3.083707d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,325) = krate(i,72)*exp(6.310961d-1*(lnTgas(i)-1d0) &
            - 5.418875d-4*Tgas(i) &
            + 7.125405d-8*Tgas2(i) &
            - 5.612989d-12*Tgas3(i) &
            + 2.029648d-16*Tgas4(i) &
            - 4.192898d4*invTgas(i) &
            + 1.551566d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,325) = 0d0
      end if
    end do

    !HSO3 + M -> SO2 + OH + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,326) = krate(i,73)*exp(4.482077d0*(lnTgas(i)-1d0) &
            - 9.672905d-3*Tgas(i) &
            + 5.280526d-6*Tgas2(i) &
            - 2.124199d-9*Tgas3(i) &
            + 4.034084d-13*Tgas4(i) &
            - 1.063905d4*invTgas(i) &
            - 5.185011d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,326) = krate(i,73)*exp(2.856084d-2*(lnTgas(i)-1d0) &
            - 4.957808d-4*Tgas(i) &
            + 7.049484d-8*Tgas2(i) &
            - 5.55515d-12*Tgas3(i) &
            + 1.947841d-16*Tgas4(i) &
            - 1.159243d4*invTgas(i) &
            + 1.639613d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,326) = 0d0
      end if
    end do

    !SO2 -> CH3SCH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,327) = krate(i,74)*exp(5.597728d0*(lnTgas(i)-1d0) &
            - 1.118526d-3*Tgas(i) &
            + 6.816712d-6*Tgas2(i) &
            - 3.990816d-9*Tgas3(i) &
            + 9.087257d-13*Tgas4(i) &
            - 3.408447d4*invTgas(i) &
            - 6.022887d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,327) = krate(i,74)*exp(3.920635d0*(lnTgas(i)-1d0) &
            + 7.508924d-3*Tgas(i) &
            - 8.588814d-7*Tgas2(i) &
            + 6.733831d-11*Tgas3(i) &
            - 2.3796d-15*Tgas4(i) &
            - 3.395525d4*invTgas(i) &
            - 3.431551d-1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,327) = 0d0
      end if
    end do

    !SO2 -> CH3SCH3 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,328) = krate(i,75)*exp(4.774011d0*(lnTgas(i)-1d0) &
            - 1.557652d-3*Tgas(i) &
            + 7.154449d-6*Tgas2(i) &
            - 4.178225d-9*Tgas3(i) &
            + 9.461989d-13*Tgas4(i) &
            - 5.983783d4*invTgas(i) &
            - 3.866955d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,328) = krate(i,75)*exp(3.625742d0*(lnTgas(i)-1d0) &
            + 6.941559d-3*Tgas(i) &
            - 8.105797d-7*Tgas2(i) &
            + 6.424539d-11*Tgas3(i) &
            - 2.282433d-15*Tgas4(i) &
            - 5.948346d4*invTgas(i) &
            - 1.265807d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,328) = 0d0
      end if
    end do

    !SO2 + CH4O3S + M -> CH3SCH3 + OH + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,329) = krate(i,76)*exp(3.171767d-1*(lnTgas(i)-1d0) &
            - 2.342044d-3*Tgas(i) &
            - 6.420484d-7*Tgas2(i) &
            + 8.147534d-10*Tgas3(i) &
            - 2.199762d-13*Tgas4(i) &
            - 4.031441d4*invTgas(i) &
            - 8.072663d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,329) = krate(i,76)*exp(-2.545704d0*(lnTgas(i)-1d0) &
            - 2.859464d-4*Tgas(i) &
            + 5.634379d-8*Tgas2(i) &
            - 5.532956d-12*Tgas3(i) &
            + 2.123002d-16*Tgas4(i) &
            - 4.130451d4*invTgas(i) &
            + 7.676252d0)
      else
        krate(i,329) = 0d0
      end if
    end do

    !SO4 -> H2SO4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,330) = krate(i,77)*exp(3.530814d0*(lnTgas(i)-1d0) &
            - 2.295174d-3*Tgas(i) &
            + 1.201235d-6*Tgas2(i) &
            - 2.381115d-10*Tgas3(i) &
            - 4.367358d-15*Tgas4(i) &
            + 6.002935d4*invTgas(i) &
            - 1.51445d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,330) = krate(i,77)*exp(1.251161d0*(lnTgas(i)-1d0) &
            + 1.315445d-3*Tgas(i) &
            - 1.302866d-7*Tgas2(i) &
            + 9.052629d-12*Tgas3(i) &
            - 2.887276d-16*Tgas4(i) &
            + 5.945508d4*invTgas(i) &
            - 3.632933d0)
      else
        krate(i,330) = 0d0
      end if
    end do

    !CH3 + OH -> O(1D) + CH4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,331) = krate(i,78)*exp(-9.6876d-4*(lnTgas(i)-1d0) &
            - 6.7128d-3*Tgas(i) &
            + 6.549928d-6*Tgas2(i) &
            - 3.2248d-9*Tgas3(i) &
            + 6.590328d-13*Tgas4(i) &
            - 2.19596d4*invTgas(i) &
            - 1.585919d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,331) = krate(i,78)*exp(-1.41118d0*(lnTgas(i)-1d0) &
            + 1.355588d-3*Tgas(i) &
            - 1.873897d-7*Tgas2(i) &
            + 1.60158d-11*Tgas3(i) &
            - 5.895982d-16*Tgas4(i) &
            - 2.16921d4*invTgas(i) &
            + 2.565502d0)
      else
        krate(i,331) = 0d0
      end if
    end do

    !CH3O + H -> O(1D) + CH4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,332) = krate(i,79)*exp(1.43639d0*(lnTgas(i)-1d0) &
            - 5.447718d-3*Tgas(i) &
            + 1.95325d-6*Tgas2(i) &
            - 1.573045d-10*Tgas3(i) &
            - 8.246415d-14*Tgas4(i) &
            - 1.496982d4*invTgas(i) &
            - 6.142107d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,332) = krate(i,79)*exp(-2.852322d0*(lnTgas(i)-1d0) &
            + 1.087508d-3*Tgas(i) &
            - 1.161445d-7*Tgas2(i) &
            + 8.622238d-12*Tgas3(i) &
            - 2.889284d-16*Tgas4(i) &
            - 1.603562d4*invTgas(i) &
            + 1.554641d1)
      else
        krate(i,332) = 0d0
      end if
    end do

    !CH2O + H2 -> O(1D) + CH4
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,333) = krate(i,80)*exp(5.101409d-1*(lnTgas(i)-1d0) &
            - 5.886128d-3*Tgas(i) &
            + 5.255125d-6*Tgas2(i) &
            - 2.618629d-9*Tgas3(i) &
            + 5.604206d-13*Tgas4(i) &
            - 5.704834d4*invTgas(i) &
            - 1.30219d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,333) = krate(i,80)*exp(-1.696922d0*(lnTgas(i)-1d0) &
            + 1.298314d-3*Tgas(i) &
            - 1.661587d-7*Tgas2(i) &
            + 1.334764d-11*Tgas3(i) &
            - 4.714264d-16*Tgas4(i) &
            - 5.726117d4*invTgas(i) &
            + 8.115178d0)
      else
        krate(i,333) = 0d0
      end if
    end do

    !CH3O + CH3O + O2 -> CH3O2 + CH3O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,334) = krate(i,81)*exp(-1.262676d0*(lnTgas(i)-1d0) &
            - 9.905654d-4*Tgas(i) &
            + 1.718222d-6*Tgas2(i) &
            - 9.264509d-10*Tgas3(i) &
            + 1.941291d-13*Tgas4(i) &
            + 1.80955d3*invTgas(i) &
            - 1.117247d1)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,334) = krate(i,81)*exp(-2.065936d0*(lnTgas(i)-1d0) &
            + 1.352754d-3*Tgas(i) &
            - 1.569633d-7*Tgas2(i) &
            + 1.172223d-11*Tgas3(i) &
            - 3.880138d-16*Tgas4(i) &
            + 1.635689d3*invTgas(i) &
            - 7.464933d0)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,334) = 0d0
      end if
    end do

    !CH3OH + CH2O + O2 -> CH3O2 + CH3O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,335) = krate(i,82)*exp(-4.291299d0*(lnTgas(i)-1d0) &
            + 9.308139d-3*Tgas(i) &
            - 3.482714d-6*Tgas2(i) &
            + 6.69499d-10*Tgas3(i) &
            - 9.876835d-16*Tgas4(i) &
            - 4.079707d4*invTgas(i) &
            + 2.26688d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,335) = krate(i,82)*exp(7.528545d-1*(lnTgas(i)-1d0) &
            + 5.386365d-4*Tgas(i) &
            - 7.606508d-8*Tgas2(i) &
            + 6.118677d-12*Tgas3(i) &
            - 2.117243d-16*Tgas4(i) &
            - 3.969616d4*invTgas(i) &
            - 2.26082d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,335) = 0d0
      end if
    end do

    !CO2 + M -> O + CO + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,336) = krate(i,83)*exp(4.390988d0*(lnTgas(i)-1d0) &
            - 6.436901d-3*Tgas(i) &
            + 2.463657d-6*Tgas2(i) &
            - 6.398634d-10*Tgas3(i) &
            + 6.755604d-14*Tgas4(i) &
            - 6.315014d4*invTgas(i) &
            - 4.340561d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,336) = krate(i,83)*exp(9.556118d-1*(lnTgas(i)-1d0) &
            - 7.085225d-4*Tgas(i) &
            + 8.431887d-8*Tgas2(i) &
            - 6.381516d-12*Tgas3(i) &
            + 1.992179d-16*Tgas4(i) &
            - 6.39848d4*invTgas(i) &
            + 1.287429d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,336) = 0d0
      end if
    end do

    !CHO + M -> H + CO + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,337) = krate(i,84)*exp(1.841987d0*(lnTgas(i)-1d0) &
            + 1.355199d-3*Tgas(i) &
            - 2.164369d-6*Tgas2(i) &
            + 1.19425d-9*Tgas3(i) &
            - 2.639293d-13*Tgas4(i) &
            - 7.257162d3*invTgas(i) &
            - 2.466222d-1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,337) = krate(i,84)*exp(1.62847d0*(lnTgas(i)-1d0) &
            - 5.855326d-4*Tgas(i) &
            + 3.086835d-8*Tgas2(i) &
            - 2.230192d-12*Tgas3(i) &
            + 1.369954d-16*Tgas4(i) &
            - 7.554114d3*invTgas(i) &
            + 1.989644d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,337) = 0d0
      end if
    end do

    !H2 + CO -> H + CHO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,338) = krate(i,85)*exp(8.136815d-1*(lnTgas(i)-1d0) &
            - 5.34546d-3*Tgas(i) &
            + 5.410727d-6*Tgas2(i) &
            - 2.874018d-9*Tgas3(i) &
            + 6.327352d-13*Tgas4(i) &
            - 4.460809d4*invTgas(i) &
            - 1.329754d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,338) = krate(i,85)*exp(4.386638d-1*(lnTgas(i)-1d0) &
            + 1.722286d-4*Tgas(i) &
            - 6.467958d-9*Tgas2(i) &
            + 9.460219d-13*Tgas3(i) &
            - 1.025552d-16*Tgas4(i) &
            - 4.420627d4*invTgas(i) &
            - 1.858681d0)
      else
        krate(i,338) = 0d0
      end if
    end do

    !CH2O + CO -> CHO + CHO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,339) = krate(i,86)*exp(1.018356d-1*(lnTgas(i)-1d0) &
            + 1.938591d-3*Tgas(i) &
            - 1.722127d-6*Tgas2(i) &
            + 8.477932d-10*Tgas3(i) &
            - 1.762258d-13*Tgas4(i) &
            - 3.64681d4*invTgas(i) &
            + 2.50549d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,339) = krate(i,86)*exp(1.622018d0*(lnTgas(i)-1d0) &
            - 1.249674d-3*Tgas(i) &
            + 2.323916d-7*Tgas2(i) &
            - 1.946645d-11*Tgas3(i) &
            + 5.918528d-16*Tgas4(i) &
            - 3.612166d4*invTgas(i) &
            - 4.897636d0)
      else
        krate(i,339) = 0d0
      end if
    end do

    !H2O + CO -> OH + CHO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,340) = krate(i,87)*exp(4.513616d-1*(lnTgas(i)-1d0) &
            - 1.537532d-3*Tgas(i) &
            + 1.847085d-6*Tgas2(i) &
            - 1.060187d-9*Tgas3(i) &
            + 2.434907d-13*Tgas4(i) &
            - 5.187912d4*invTgas(i) &
            + 5.449499d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,340) = krate(i,87)*exp(1.033021d0*(lnTgas(i)-1d0) &
            - 3.473518d-4*Tgas(i) &
            + 4.909309d-8*Tgas2(i) &
            - 2.133445d-12*Tgas3(i) &
            - 4.469039d-17*Tgas4(i) &
            - 5.150325d4*invTgas(i) &
            - 3.473931d0)
      else
        krate(i,340) = 0d0
      end if
    end do

    !H + CO2 -> O + CHO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,341) = krate(i,88)*exp(2.549d0*(lnTgas(i)-1d0) &
            - 7.792101d-3*Tgas(i) &
            + 4.628026d-6*Tgas2(i) &
            - 1.834114d-9*Tgas3(i) &
            + 3.314854d-13*Tgas4(i) &
            - 5.589298d4*invTgas(i) &
            - 4.093939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,341) = krate(i,88)*exp(-6.728587d-1*(lnTgas(i)-1d0) &
            - 1.2299d-4*Tgas(i) &
            + 5.345052d-8*Tgas2(i) &
            - 4.151324d-12*Tgas3(i) &
            + 6.222247d-17*Tgas4(i) &
            - 5.643069d4*invTgas(i) &
            + 1.088464d1)
      else
        krate(i,341) = 0d0
      end if
    end do

    !OH + CO -> O + CHO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,342) = krate(i,89)*exp(-1.657045d-1*(lnTgas(i)-1d0) &
            - 1.794326d-3*Tgas(i) &
            + 2.502106d-6*Tgas2(i) &
            - 1.381659d-9*Tgas3(i) &
            + 3.014026d-13*Tgas4(i) &
            - 4.396986d4*invTgas(i) &
            + 1.955871d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,342) = krate(i,89)*exp(5.766362d-1*(lnTgas(i)-1d0) &
            + 1.8168d-5*Tgas(i) &
            + 1.74333d-8*Tgas2(i) &
            - 8.627292d-13*Tgas3(i) &
            - 3.982814d-17*Tgas4(i) &
            - 4.344775d4*invTgas(i) &
            - 3.358979d0)
      else
        krate(i,342) = 0d0
      end if
    end do

    !H2 + CHO -> H + CH2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,343) = krate(i,90)*exp(7.118459d-1*(lnTgas(i)-1d0) &
            - 7.284051d-3*Tgas(i) &
            + 7.132854d-6*Tgas2(i) &
            - 3.721811d-9*Tgas3(i) &
            + 8.08961d-13*Tgas4(i) &
            - 8.139988d3*invTgas(i) &
            - 3.835244d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,343) = krate(i,90)*exp(-1.183355d0*(lnTgas(i)-1d0) &
            + 1.421902d-3*Tgas(i) &
            - 2.388595d-7*Tgas2(i) &
            + 2.041247d-11*Tgas3(i) &
            - 6.944079d-16*Tgas4(i) &
            - 8.084613d3*invTgas(i) &
            + 3.038954d0)
      else
        krate(i,343) = 0d0
      end if
    end do

    !OH + CHO -> O + CH2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,344) = krate(i,91)*exp(-2.675401d-1*(lnTgas(i)-1d0) &
            - 3.732916d-3*Tgas(i) &
            + 4.224233d-6*Tgas2(i) &
            - 2.229452d-9*Tgas3(i) &
            + 4.776283d-13*Tgas4(i) &
            - 7.501754d3*invTgas(i) &
            - 5.496187d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,344) = krate(i,91)*exp(-1.045382d0*(lnTgas(i)-1d0) &
            + 1.267842d-3*Tgas(i) &
            - 2.149583d-7*Tgas2(i) &
            + 1.860372d-11*Tgas3(i) &
            - 6.316809d-16*Tgas4(i) &
            - 7.326092d3*invTgas(i) &
            + 1.538656d0)
      else
        krate(i,344) = 0d0
      end if
    end do

    !H + OH -> O(1D) + H2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,345) = krate(i,92)*exp(-1.647715d0*(lnTgas(i)-1d0) &
            + 5.19088d-3*Tgas(i) &
            - 4.015856d-6*Tgas2(i) &
            + 2.003062d-9*Tgas3(i) &
            - 4.369742d-13*Tgas4(i) &
            - 2.223604d4*invTgas(i) &
            + 5.850537d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,345) = krate(i,92)*exp(8.802017d-2*(lnTgas(i)-1d0) &
            - 1.335215d-4*Tgas(i) &
            + 2.292629d-8*Tgas2(i) &
            - 1.991642d-12*Tgas3(i) &
            + 7.662483d-17*Tgas4(i) &
            - 2.20141d4*invTgas(i) &
            - 1.772083d0)
      else
        krate(i,345) = 0d0
      end if
    end do

    !H + H2O -> OH + H2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,346) = krate(i,93)*exp(-3.623198d-1*(lnTgas(i)-1d0) &
            + 3.807928d-3*Tgas(i) &
            - 3.563642d-6*Tgas2(i) &
            + 1.813831d-9*Tgas3(i) &
            - 3.892445d-13*Tgas4(i) &
            - 7.271029d3*invTgas(i) &
            + 1.874704d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,346) = krate(i,93)*exp(5.943572d-1*(lnTgas(i)-1d0) &
            - 5.195803d-4*Tgas(i) &
            + 5.556105d-8*Tgas2(i) &
            - 3.079467d-12*Tgas3(i) &
            + 5.786477d-17*Tgas4(i) &
            - 7.296976d3*invTgas(i) &
            - 1.615249d0)
      else
        krate(i,346) = 0d0
      end if
    end do

    !SO2 + OH -> SO + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,347) = krate(i,94)*exp(2.536015d-1*(lnTgas(i)-1d0) &
            - 3.476407d-3*Tgas(i) &
            + 3.286498d-6*Tgas2(i) &
            - 1.745856d-9*Tgas3(i) &
            + 3.886996d-13*Tgas4(i) &
            - 3.336001d4*invTgas(i) &
            + 2.217038d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,347) = krate(i,94)*exp(-8.153549d-2*(lnTgas(i)-1d0) &
            - 2.641227d-4*Tgas(i) &
            + 9.790948d-8*Tgas2(i) &
            - 1.206866d-11*Tgas3(i) &
            + 5.319727d-16*Tgas4(i) &
            - 3.321134d4*invTgas(i) &
            + 2.678966d0)
      else
        krate(i,347) = 0d0
      end if
    end do

    !S2O2 + M -> SO + SO + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,348) = krate(i,95)*exp(4.218492d0*(lnTgas(i)-1d0) &
            - 7.750528d-3*Tgas(i) &
            + 5.972407d-6*Tgas2(i) &
            - 2.916257d-9*Tgas3(i) &
            + 6.215523d-13*Tgas4(i) &
            - 7.062465d3*invTgas(i) &
            + 5.561934d-1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,348) = krate(i,95)*exp(1.913866d0*(lnTgas(i)-1d0) &
            - 1.228811d-4*Tgas(i) &
            + 6.787752d-8*Tgas2(i) &
            - 7.869408d-12*Tgas3(i) &
            + 3.452773d-16*Tgas4(i) &
            - 7.308169d3*invTgas(i) &
            + 1.040676d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,348) = 0d0
      end if
    end do

    !SO2 + S2O -> SO + S2O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,349) = krate(i,96)*exp(-5.621238d-2*(lnTgas(i)-1d0) &
            - 2.302379d-3*Tgas(i) &
            + 5.295561d-7*Tgas2(i) &
            - 4.608726d-11*Tgas3(i) &
            - 7.753252d-15*Tgas4(i) &
            - 3.646489d4*invTgas(i) &
            - 1.603623d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,349) = krate(i,96)*exp(-1.415293d0*(lnTgas(i)-1d0) &
            - 6.510044d-4*Tgas(i) &
            + 1.066223d-7*Tgas2(i) &
            - 1.018498d-11*Tgas3(i) &
            + 4.020149d-16*Tgas4(i) &
            - 3.687813d4*invTgas(i) &
            + 5.566236d0)
      else
        krate(i,349) = 0d0
      end if
    end do

    !S + SO2 -> SO + SO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,350) = krate(i,97)*exp(1.245127d0*(lnTgas(i)-1d0) &
            - 5.85334d-3*Tgas(i) &
            + 4.838645d-6*Tgas2(i) &
            - 2.534938d-9*Tgas3(i) &
            + 5.709236d-13*Tgas4(i) &
            - 3.477366d3*invTgas(i) &
            - 1.301006d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,350) = krate(i,97)*exp(-3.257153d-1*(lnTgas(i)-1d0) &
            - 2.068308d-4*Tgas(i) &
            + 6.559971d-8*Tgas2(i) &
            - 7.619057d-12*Tgas3(i) &
            + 3.372259d-16*Tgas4(i) &
            - 3.64818d3*invTgas(i) &
            + 5.31976d0)
      else
        krate(i,350) = 0d0
      end if
    end do

    !SO2 + SO2 -> SO + SO3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,351) = krate(i,98)*exp(-1.356409d0*(lnTgas(i)-1d0) &
            + 4.533275d-3*Tgas(i) &
            - 2.987308d-6*Tgas2(i) &
            + 1.326786d-9*Tgas3(i) &
            - 2.647206d-13*Tgas4(i) &
            - 2.448347d4*invTgas(i) &
            + 3.532017d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,351) = krate(i,98)*exp(4.972483d-1*(lnTgas(i)-1d0) &
            - 1.22775d-4*Tgas(i) &
            + 3.466989d-8*Tgas2(i) &
            - 4.159092d-12*Tgas3(i) &
            + 1.750724d-16*Tgas4(i) &
            - 2.417516d4*invTgas(i) &
            - 5.027125d0)
      else
        krate(i,351) = 0d0
      end if
    end do

    !S + H2S -> SH + SH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,352) = krate(i,99)*exp(9.318368d-1*(lnTgas(i)-1d0) &
            + 1.795534d-3*Tgas(i) &
            - 3.288836d-6*Tgas2(i) &
            + 2.109145d-9*Tgas3(i) &
            - 5.164187d-13*Tgas4(i) &
            - 2.98255d3*invTgas(i) &
            - 3.558532d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,352) = krate(i,99)*exp(2.049013d-1*(lnTgas(i)-1d0) &
            - 2.852241d-4*Tgas(i) &
            + 2.719609d-8*Tgas2(i) &
            - 2.36317d-12*Tgas3(i) &
            + 9.936002d-17*Tgas4(i) &
            - 3.426632d3*invTgas(i) &
            + 1.539824d0)
      else
        krate(i,352) = 0d0
      end if
    end do

    !H2 + S -> SH + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,353) = krate(i,100)*exp(1.523081d0*(lnTgas(i)-1d0) &
            - 4.757308d-3*Tgas(i) &
            + 3.470485d-6*Tgas2(i) &
            - 1.571912d-9*Tgas3(i) &
            + 3.133023d-13*Tgas4(i) &
            - 9.788345d3*invTgas(i) &
            - 5.174301d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,353) = krate(i,100)*exp(-2.807005d-1*(lnTgas(i)-1d0) &
            + 4.712474d-4*Tgas(i) &
            - 8.548808d-8*Tgas2(i) &
            + 7.591693d-12*Tgas3(i) &
            - 2.743495d-16*Tgas4(i) &
            - 9.991314d3*invTgas(i) &
            + 2.746461d0)
      else
        krate(i,353) = 0d0
      end if
    end do

    !H2S + CHO -> SH + CH2O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,354) = krate(i,101)*exp(1.206012d-1*(lnTgas(i)-1d0) &
            - 7.312091d-4*Tgas(i) &
            + 3.73533d-7*Tgas2(i) &
            - 4.075394d-11*Tgas3(i) &
            - 2.076002d-14*Tgas4(i) &
            - 1.334192d3*invTgas(i) &
            - 2.219475d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,354) = krate(i,101)*exp(-6.977528d-1*(lnTgas(i)-1d0) &
            + 6.654306d-4*Tgas(i) &
            - 1.261754d-7*Tgas2(i) &
            + 1.045761d-11*Tgas3(i) &
            - 3.206984d-16*Tgas4(i) &
            - 1.519931d3*invTgas(i) &
            + 1.832318d0)
      else
        krate(i,354) = 0d0
      end if
    end do

    !S2 + M -> S + S + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,355) = krate(i,102)*exp(1.757146d0*(lnTgas(i)-1d0) &
            + 2.278676d-3*Tgas(i) &
            - 3.728805d-6*Tgas2(i) &
            + 2.355543d-9*Tgas3(i) &
            - 5.772374d-13*Tgas4(i) &
            - 5.057956d4*invTgas(i) &
            + 2.326109d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,355) = krate(i,102)*exp(1.926233d0*(lnTgas(i)-1d0) &
            - 9.555358d-4*Tgas(i) &
            + 1.277824d-7*Tgas2(i) &
            - 1.045632d-11*Tgas3(i) &
            + 3.530205d-16*Tgas4(i) &
            - 5.071914d4*invTgas(i) &
            + 2.632804d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,355) = 0d0
      end if
    end do

    !S3 + M -> S + S2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,356) = krate(i,103)*exp(2.520361d0*(lnTgas(i)-1d0) &
            - 4.394677d-3*Tgas(i) &
            + 2.278692d-6*Tgas2(i) &
            - 8.539568d-10*Tgas3(i) &
            + 1.531129d-13*Tgas4(i) &
            - 3.090909d4*invTgas(i) &
            + 2.134197d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,356) = krate(i,103)*exp(1.788388d-1*(lnTgas(i)-1d0) &
            - 5.55983d-5*Tgas(i) &
            + 3.147439d-8*Tgas2(i) &
            - 3.465545d-12*Tgas3(i) &
            + 1.516089d-16*Tgas4(i) &
            - 3.146634d4*invTgas(i) &
            + 1.37352d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,356) = 0d0
      end if
    end do

    !S4 + M -> S + S3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,357) = krate(i,104)*exp(3.370273d0*(lnTgas(i)-1d0) &
            - 6.808341d-3*Tgas(i) &
            + 3.51533d-6*Tgas2(i) &
            - 1.309608d-9*Tgas3(i) &
            + 2.285942d-13*Tgas4(i) &
            - 3.385096d4*invTgas(i) &
            + 2.158188d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,357) = krate(i,104)*exp(2.845701d-1*(lnTgas(i)-1d0) &
            - 4.678589d-4*Tgas(i) &
            + 7.040091d-8*Tgas2(i) &
            - 6.131962d-12*Tgas3(i) &
            + 2.242026d-16*Tgas4(i) &
            - 3.448909d4*invTgas(i) &
            + 1.705524d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,357) = 0d0
      end if
    end do

    !S4 + M -> S2 + S2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,358) = krate(i,105)*exp(4.133488d0*(lnTgas(i)-1d0) &
            - 1.348169d-2*Tgas(i) &
            + 9.522827d-6*Tgas2(i) &
            - 4.519108d-9*Tgas3(i) &
            + 9.589445d-13*Tgas4(i) &
            - 1.41805d4*invTgas(i) &
            + 1.966275d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,358) = krate(i,105)*exp(-1.462824d0*(lnTgas(i)-1d0) &
            + 4.320787d-4*Tgas(i) &
            - 2.590707d-8*Tgas2(i) &
            + 8.588096d-13*Tgas3(i) &
            + 2.279101d-17*Tgas4(i) &
            - 1.523629d4*invTgas(i) &
            + 2.815763d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,358) = 0d0
      end if
    end do

    !S8 + M -> S4 + S4 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,359) = krate(i,106)*exp(-8.890915d-1*(lnTgas(i)-1d0) &
            - 1.019551d-2*Tgas(i) &
            + 1.12212d-5*Tgas2(i) &
            - 7.091641d-9*Tgas3(i) &
            + 1.763928d-12*Tgas4(i) &
            - 2.117277d4*invTgas(i) &
            + 2.742703d1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,359) = krate(i,106)*exp(-2.175131d0*(lnTgas(i)-1d0) &
            - 1.67668d-3*Tgas(i) &
            + 3.655858d-7*Tgas2(i) &
            - 3.938726d-11*Tgas3(i) &
            + 1.677997d-15*Tgas4(i) &
            - 2.154344d4*invTgas(i) &
            + 3.244209d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,359) = 0d0
      end if
    end do

    !S + S + M -> S2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,360) = krate(i,107)*exp(-1.757146d0*(lnTgas(i)-1d0) &
            - 2.278676d-3*Tgas(i) &
            + 3.728805d-6*Tgas2(i) &
            - 2.355543d-9*Tgas3(i) &
            + 5.772374d-13*Tgas4(i) &
            + 5.057956d4*invTgas(i) &
            - 2.326109d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,360) = krate(i,107)*exp(-1.926233d0*(lnTgas(i)-1d0) &
            + 9.555358d-4*Tgas(i) &
            - 1.277824d-7*Tgas2(i) &
            + 1.045632d-11*Tgas3(i) &
            - 3.530205d-16*Tgas4(i) &
            + 5.071914d4*invTgas(i) &
            - 2.632804d0)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,360) = 0d0
      end if
    end do

    !S + SO -> S2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,361) = krate(i,108)*exp(1.097821d-1*(lnTgas(i)-1d0) &
            - 3.673751d-4*Tgas(i) &
            + 5.268929d-7*Tgas2(i) &
            - 3.77365d-10*Tgas3(i) &
            + 1.041936d-13*Tgas4(i) &
            - 1.153022d4*invTgas(i) &
            - 5.767928d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,361) = krate(i,108)*exp(-4.721737d-1*(lnTgas(i)-1d0) &
            + 4.977041d-4*Tgas(i) &
            - 8.745815d-8*Tgas2(i) &
            + 8.303293d-12*Tgas3(i) &
            - 3.122093d-16*Tgas4(i) &
            - 1.173682d4*invTgas(i) &
            + 2.535966d0)
      else
        krate(i,361) = 0d0
      end if
    end do

    !CH3O -> O + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,362) = krate(i,109)*exp(3.113642d0*(lnTgas(i)-1d0) &
            + 8.259561d-4*Tgas(i) &
            - 4.258941d-6*Tgas2(i) &
            + 2.880087d-9*Tgas3(i) &
            - 7.040237d-13*Tgas4(i) &
            - 4.423725d4*invTgas(i) &
            - 2.84694d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,362) = krate(i,109)*exp(7.639652d-1*(lnTgas(i)-1d0) &
            - 8.354445d-4*Tgas(i) &
            + 1.195469d-7*Tgas2(i) &
            - 1.048648d-11*Tgas3(i) &
            + 3.97837d-16*Tgas4(i) &
            - 4.534539d4*invTgas(i) &
            + 1.161157d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,362) = 0d0
      end if
    end do

    !CH2O + H -> O + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,363) = krate(i,110)*exp(-4.682763d-1*(lnTgas(i)-1d0) &
            + 4.377806d-3*Tgas(i) &
            - 4.203424d-6*Tgas2(i) &
            + 2.09853d-9*Tgas3(i) &
            - 4.299448d-13*Tgas4(i) &
            - 3.445051d4*invTgas(i) &
            + 3.569354d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,363) = krate(i,110)*exp(-1.477691d-1*(lnTgas(i)-1d0) &
            - 2.113349d-4*Tgas(i) &
            + 4.513223d-8*Tgas2(i) &
            - 4.476912d-12*Tgas3(i) &
            + 1.808988d-16*Tgas4(i) &
            - 3.481055d4*invTgas(i) &
            + 4.049378d0)
      else
        krate(i,363) = 0d0
      end if
    end do

    !CH4 + M -> H + CH3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,364) = krate(i,111)*exp(1.008922d0*(lnTgas(i)-1d0) &
            + 7.913419d-3*Tgas(i) &
            - 7.319426d-6*Tgas2(i) &
            + 3.548094d-9*Tgas3(i) &
            - 7.272011d-13*Tgas4(i) &
            - 5.21417d4*invTgas(i) &
            + 5.86008d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,364) = krate(i,111)*exp(3.566335d0*(lnTgas(i)-1d0) &
            - 1.902414d-3*Tgas(i) &
            + 2.347164d-7*Tgas2(i) &
            - 1.929161d-11*Tgas3(i) &
            + 7.006632d-16*Tgas4(i) &
            - 5.208239d4*invTgas(i) &
            - 4.206622d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,364) = 0d0
      end if
    end do

    !CH3O + O2 -> O3 + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,365) = krate(i,112)*exp(-4.296995d-1*(lnTgas(i)-1d0) &
            + 4.990878d-3*Tgas(i) &
            - 4.699234d-6*Tgas2(i) &
            + 2.336604d-9*Tgas3(i) &
            - 4.838065d-13*Tgas4(i) &
            - 3.204343d4*invTgas(i) &
            - 2.740732d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,365) = krate(i,112)*exp(6.889659d0*(lnTgas(i)-1d0) &
            - 7.116208d-3*Tgas(i) &
            + 1.475006d-6*Tgas2(i) &
            - 1.602765d-10*Tgas3(i) &
            + 6.790563d-15*Tgas4(i) &
            - 3.001093d4*invTgas(i) &
            - 3.760842d1)
      else
        krate(i,365) = 0d0
      end if
    end do

    !CH4 + HO2 -> H2O2 + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,366) = krate(i,113)*exp(-1.477724d0*(lnTgas(i)-1d0) &
            + 9.864285d-3*Tgas(i) &
            - 7.905735d-6*Tgas2(i) &
            + 3.681436d-9*Tgas3(i) &
            - 7.373387d-13*Tgas4(i) &
            - 8.697276d3*invTgas(i) &
            + 5.863834d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,366) = krate(i,113)*exp(1.47382d0*(lnTgas(i)-1d0) &
            - 8.163719d-4*Tgas(i) &
            + 7.60214d-8*Tgas2(i) &
            - 4.396141d-12*Tgas3(i) &
            + 1.220064d-16*Tgas4(i) &
            - 8.570528d3*invTgas(i) &
            - 6.052645d0)
      else
        krate(i,366) = 0d0
      end if
    end do

    !CH3O + H -> OH + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,367) = krate(i,114)*exp(1.437359d0*(lnTgas(i)-1d0) &
            + 1.265082d-3*Tgas(i) &
            - 4.596678d-6*Tgas2(i) &
            + 3.067495d-9*Tgas3(i) &
            - 7.41497d-13*Tgas4(i) &
            + 6.98977d3*invTgas(i) &
            - 4.556189d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,367) = krate(i,114)*exp(-1.441141d0*(lnTgas(i)-1d0) &
            - 2.680799d-4*Tgas(i) &
            + 7.124526d-8*Tgas2(i) &
            - 7.393561d-12*Tgas3(i) &
            + 3.006697d-16*Tgas4(i) &
            + 5.656478d3*invTgas(i) &
            + 1.298091d1)
      else
        krate(i,367) = 0d0
      end if
    end do

    !CH4 + O -> OH + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,368) = krate(i,115)*exp(-6.673605d-1*(lnTgas(i)-1d0) &
            + 8.352546d-3*Tgas(i) &
            - 7.657163d-6*Tgas2(i) &
            + 3.735503d-9*Tgas3(i) &
            - 7.646744d-13*Tgas4(i) &
            - 9.146774d2*invTgas(i) &
            + 4.150831d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,368) = krate(i,115)*exp(1.361228d0*(lnTgas(i)-1d0) &
            - 1.335049d-3*Tgas(i) &
            + 1.864147d-7*Tgas2(i) &
            - 1.619869d-11*Tgas3(i) &
            + 6.03496d-16*Tgas4(i) &
            - 1.080523d3*invTgas(i) &
            - 2.837287d0)
      else
        krate(i,368) = 0d0
      end if
    end do

    !CH3OH + M -> OH + CH3 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,369) = krate(i,116)*exp(1.990653d0*(lnTgas(i)-1d0) &
            + 8.011937d-3*Tgas(i) &
            - 9.853131d-6*Tgas2(i) &
            + 5.445002d-9*Tgas3(i) &
            - 1.210693d-12*Tgas4(i) &
            - 4.540359d4*invTgas(i) &
            + 2.466867d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,369) = krate(i,116)*exp(2.289383d0*(lnTgas(i)-1d0) &
            - 1.706307d-3*Tgas(i) &
            + 2.265582d-7*Tgas2(i) &
            - 1.900669d-11*Tgas3(i) &
            + 6.938974d-16*Tgas4(i) &
            - 4.62102d4*invTgas(i) &
            + 5.399839d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,369) = 0d0
      end if
    end do

    !CH3O + OH -> HO2 + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,370) = krate(i,117)*exp(2.551885d-1*(lnTgas(i)-1d0) &
            + 1.291588d-3*Tgas(i) &
            - 2.609176d-6*Tgas2(i) &
            + 1.69099d-9*Tgas3(i) &
            - 4.132039d-13*Tgas4(i) &
            - 1.201011d4*invTgas(i) &
            - 1.078213d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,370) = krate(i,117)*exp(-4.459147d-1*(lnTgas(i)-1d0) &
            - 4.349047d-4*Tgas(i) &
            + 1.115324d-7*Tgas2(i) &
            - 1.278306d-11*Tgas3(i) &
            + 5.517725d-16*Tgas4(i) &
            - 1.245259d4*invTgas(i) &
            + 3.80201d0)
      else
        krate(i,370) = 0d0
      end if
    end do

    !CH4 + O2 -> HO2 + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,371) = krate(i,118)*exp(-9.717359d-1*(lnTgas(i)-1d0) &
            + 7.037226d-3*Tgas(i) &
            - 5.434261d-6*Tgas2(i) &
            + 2.331836d-9*Tgas3(i) &
            - 4.247749d-13*Tgas4(i) &
            - 2.7996d4*invTgas(i) &
            + 6.365749d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,371) = krate(i,118)*exp(1.577661d0*(lnTgas(i)-1d0) &
            - 1.290009d-3*Tgas(i) &
            + 2.005285d-7*Tgas2(i) &
            - 1.938445d-11*Tgas3(i) &
            + 7.744328d-16*Tgas4(i) &
            - 2.785572d4*invTgas(i) &
            - 4.217625d0)
      else
        krate(i,371) = 0d0
      end if
    end do

    !CH4 + CO -> CHO + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,372) = krate(i,119)*exp(-8.33065d-1*(lnTgas(i)-1d0) &
            + 6.55822d-3*Tgas(i) &
            - 5.155057d-6*Tgas2(i) &
            + 2.353844d-9*Tgas3(i) &
            - 4.632718d-13*Tgas4(i) &
            - 4.488454d4*invTgas(i) &
            + 6.106702d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,372) = krate(i,119)*exp(1.937864d0*(lnTgas(i)-1d0) &
            - 1.316881d-3*Tgas(i) &
            + 2.03848d-7*Tgas2(i) &
            - 1.706142d-11*Tgas3(i) &
            + 5.636678d-16*Tgas4(i) &
            - 4.452827d4*invTgas(i) &
            - 6.196266d0)
      else
        krate(i,372) = 0d0
      end if
    end do

    !CH4 + CH2 -> CH3 + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,373) = krate(i,120)*exp(-2.027211d0*(lnTgas(i)-1d0) &
            + 1.014224d-2*Tgas(i) &
            - 7.769155d-6*Tgas2(i) &
            + 3.549074d-9*Tgas3(i) &
            - 7.005324d-13*Tgas4(i) &
            + 7.27547d3*invTgas(i) &
            + 8.727032d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,373) = krate(i,120)*exp(9.094383d-1*(lnTgas(i)-1d0) &
            - 4.514574d-4*Tgas(i) &
            + 4.156442d-8*Tgas2(i) &
            - 3.147845d-12*Tgas3(i) &
            + 1.229239d-16*Tgas4(i) &
            + 7.385811d3*invTgas(i) &
            - 3.097765d0)
      else
        krate(i,373) = 0d0
      end if
    end do

    !H + CH2 -> CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,374) = krate(i,121)*exp(-3.036134d0*(lnTgas(i)-1d0) &
            + 2.228825d-3*Tgas(i) &
            - 4.497294d-7*Tgas2(i) &
            + 9.799592d-13*Tgas3(i) &
            + 2.666877d-14*Tgas4(i) &
            + 5.941717d4*invTgas(i) &
            + 2.866953d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,374) = krate(i,121)*exp(-2.656896d0*(lnTgas(i)-1d0) &
            + 1.450956d-3*Tgas(i) &
            - 1.93152d-7*Tgas2(i) &
            + 1.614377d-11*Tgas3(i) &
            - 5.777394d-16*Tgas4(i) &
            + 5.94682d4*invTgas(i) &
            + 1.108857d0)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,374) = 0d0
      end if
    end do

    !H2 + CH -> CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,375) = krate(i,122)*exp(-2.17691d0*(lnTgas(i)-1d0) &
            - 3.089122d-3*Tgas(i) &
            + 4.437752d-6*Tgas2(i) &
            - 2.494846d-9*Tgas3(i) &
            + 5.624003d-13*Tgas4(i) &
            + 5.33201d4*invTgas(i) &
            - 1.093759d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,375) = krate(i,122)*exp(-2.475682d0*(lnTgas(i)-1d0) &
            + 1.60294d-3*Tgas(i) &
            - 2.279502d-7*Tgas2(i) &
            + 1.938318d-11*Tgas3(i) &
            - 6.940581d-16*Tgas4(i) &
            + 5.36723d4*invTgas(i) &
            - 1.658374d0)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,375) = 0d0
      end if
    end do

    !CH2O + CH4 -> CH3O + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,376) = krate(i,123)*exp(-2.572996d0*(lnTgas(i)-1d0) &
            + 1.146527d-2*Tgas(i) &
            - 7.263909d-6*Tgas2(i) &
            + 2.766537d-9*Tgas3(i) &
            - 4.531222d-13*Tgas4(i) &
            - 4.235496d4*invTgas(i) &
            + 1.227637d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,376) = krate(i,123)*exp(2.6546d0*(lnTgas(i)-1d0) &
            - 1.278304d-3*Tgas(i) &
            + 1.603017d-7*Tgas2(i) &
            - 1.328204d-11*Tgas3(i) &
            + 4.83725d-16*Tgas4(i) &
            - 4.154755d4*invTgas(i) &
            - 1.176882d1)
      else
        krate(i,376) = 0d0
      end if
    end do

    !CH2O + CH4 -> CH2OH + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,377) = krate(i,124)*exp(-1.806478d0*(lnTgas(i)-1d0) &
            + 1.219224d-2*Tgas(i) &
            - 8.898364d-6*Tgas2(i) &
            + 3.668243d-9*Tgas3(i) &
            - 6.465305d-13*Tgas4(i) &
            - 3.752247d4*invTgas(i) &
            + 9.013084d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,377) = krate(i,124)*exp(2.989928d0*(lnTgas(i)-1d0) &
            - 2.025224d-3*Tgas(i) &
            + 2.656495d-7*Tgas2(i) &
            - 2.287236d-11*Tgas3(i) &
            + 8.607853d-16*Tgas4(i) &
            - 3.709928d4*invTgas(i) &
            - 1.164892d1)
      else
        krate(i,377) = 0d0
      end if
    end do

    !CH3O + CH3O -> CH3O2 + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,378) = krate(i,125)*exp(1.205265d0*(lnTgas(i)-1d0) &
            + 1.221149d-3*Tgas(i) &
            - 3.686399d-6*Tgas2(i) &
            + 2.524146d-9*Tgas3(i) &
            - 6.31499d-13*Tgas4(i) &
            - 1.367824d4*invTgas(i) &
            - 8.65627d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,378) = krate(i,125)*exp(-9.821593d-1*(lnTgas(i)-1d0) &
            + 1.868195d-5*Tgas(i) &
            + 3.000115d-8*Tgas2(i) &
            - 4.180776d-12*Tgas3(i) &
            + 1.953294d-16*Tgas4(i) &
            - 1.469354d4*invTgas(i) &
            + 4.664495d0)
      else
        krate(i,378) = 0d0
      end if
    end do

    !CH4 + H -> H2 + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,379) = krate(i,126)*exp(-1.646746d0*(lnTgas(i)-1d0) &
            + 1.190368d-2*Tgas(i) &
            - 1.056578d-5*Tgas2(i) &
            + 5.227862d-9*Tgas3(i) &
            - 1.096007d-12*Tgas4(i) &
            - 2.76443d2*invTgas(i) &
            + 7.436456d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,379) = krate(i,126)*exp(1.4992d0*(lnTgas(i)-1d0) &
            - 1.48911d-3*Tgas(i) &
            + 2.10316d-7*Tgas2(i) &
            - 1.800744d-11*Tgas3(i) &
            + 6.66223d-16*Tgas4(i) &
            - 3.22001d2*invTgas(i) &
            - 4.337585d0)
      else
        krate(i,379) = 0d0
      end if
    end do

    !CHO + H -> O + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,380) = krate(i,127)*exp(6.240342d-1*(lnTgas(i)-1d0) &
            - 1.144809d-3*Tgas(i) &
            + 1.328003d-7*Tgas2(i) &
            + 5.550612d-11*Tgas3(i) &
            - 1.645852d-14*Tgas4(i) &
            - 5.014241d4*invTgas(i) &
            - 1.556467d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,380) = krate(i,127)*exp(-7.413616d-1*(lnTgas(i)-1d0) &
            + 1.729149d-4*Tgas(i) &
            - 2.49757d-8*Tgas2(i) &
            + 1.075964d-12*Tgas3(i) &
            + 2.979003d-17*Tgas4(i) &
            - 5.060297d4*invTgas(i) &
            + 5.848513d0)
      else
        krate(i,380) = 0d0
      end if
    end do

    !H + H + CO -> O + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,381) = krate(i,128)*exp(-1.217953d0*(lnTgas(i)-1d0) &
            - 2.500008d-3*Tgas(i) &
            + 2.297169d-6*Tgas2(i) &
            - 1.138744d-9*Tgas3(i) &
            + 2.474708d-13*Tgas4(i) &
            - 4.288525d4*invTgas(i) &
            - 1.309844d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,381) = krate(i,128)*exp(-2.369832d0*(lnTgas(i)-1d0) &
            + 7.584475d-4*Tgas(i) &
            - 5.584406d-8*Tgas2(i) &
            + 3.306156d-12*Tgas3(i) &
            - 1.072054d-16*Tgas4(i) &
            - 4.304886d4*invTgas(i) &
            + 3.858869d0)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,381) = 0d0
      end if
    end do

    !H2 + CO -> O + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,382) = krate(i,129)*exp(1.437716d0*(lnTgas(i)-1d0) &
            - 6.490269d-3*Tgas(i) &
            + 5.543528d-6*Tgas2(i) &
            - 2.818512d-9*Tgas3(i) &
            + 6.162767d-13*Tgas4(i) &
            - 9.47505d4*invTgas(i) &
            - 2.88622d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,382) = krate(i,129)*exp(-3.026978d-1*(lnTgas(i)-1d0) &
            + 3.451434d-4*Tgas(i) &
            - 3.144366d-8*Tgas2(i) &
            + 2.021986d-12*Tgas3(i) &
            - 7.276513d-17*Tgas4(i) &
            - 9.480924d4*invTgas(i) &
            + 3.989832d0)
      else
        krate(i,382) = 0d0
      end if
    end do

    !H2 + CH -> H + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,383) = krate(i,130)*exp(8.592238d-1*(lnTgas(i)-1d0) &
            - 5.317947d-3*Tgas(i) &
            + 4.887482d-6*Tgas2(i) &
            - 2.495826d-9*Tgas3(i) &
            + 5.357315d-13*Tgas4(i) &
            - 6.097065d3*invTgas(i) &
            - 3.960711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,383) = krate(i,130)*exp(1.812142d-1*(lnTgas(i)-1d0) &
            + 1.519837d-4*Tgas(i) &
            - 3.479819d-8*Tgas2(i) &
            + 3.239415d-12*Tgas3(i) &
            - 1.163187d-16*Tgas4(i) &
            - 5.795898d3*invTgas(i) &
            - 2.767231d0)
      else
        krate(i,383) = 0d0
      end if
    end do

    !H + H + CO2 -> O2 + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,384) = krate(i,131)*exp(6.189566d-1*(lnTgas(i)-1d0) &
            - 7.155958d-3*Tgas(i) &
            + 4.187688d-6*Tgas2(i) &
            - 1.564038d-9*Tgas3(i) &
            + 2.659473d-13*Tgas4(i) &
            - 4.672693d4*invTgas(i) &
            - 6.096596d0)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,384) = krate(i,131)*exp(-2.840534d0*(lnTgas(i)-1d0) &
            + 4.054241d-4*Tgas(i) &
            + 6.346646d-9*Tgas2(i) &
            - 2.186181d-12*Tgas3(i) &
            + 7.501115d-17*Tgas4(i) &
            - 4.736566d4*invTgas(i) &
            + 1.030393d1)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,384) = 0d0
      end if
    end do

    !H2 + CO2 -> O2 + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,385) = krate(i,132)*exp(3.274625d0*(lnTgas(i)-1d0) &
            - 1.114622d-2*Tgas(i) &
            + 7.434047d-6*Tgas2(i) &
            - 3.243805d-9*Tgas3(i) &
            + 6.347532d-13*Tgas4(i) &
            - 9.859219d4*invTgas(i) &
            - 7.672972d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,385) = krate(i,132)*exp(-7.733993d-1*(lnTgas(i)-1d0) &
            - 7.879927d-6*Tgas(i) &
            + 3.074704d-8*Tgas2(i) &
            - 3.470351d-12*Tgas3(i) &
            + 1.094514d-16*Tgas4(i) &
            - 9.912604d4*invTgas(i) &
            + 1.043489d1)
      else
        krate(i,385) = 0d0
      end if
    end do

    !CO + H2O -> O2 + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,386) = krate(i,133)*exp(1.976009d-1*(lnTgas(i)-1d0) &
            - 1.340515d-3*Tgas(i) &
            + 1.744485d-6*Tgas2(i) &
            - 9.775195d-10*Tgas3(i) &
            + 2.154258d-13*Tgas4(i) &
            - 9.394009d4*invTgas(i) &
            + 2.515411d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,386) = krate(i,133)*exp(1.070453d0*(lnTgas(i)-1d0) &
            - 3.863023d-4*Tgas(i) &
            + 5.029087d-8*Tgas2(i) &
            - 3.261222d-12*Tgas3(i) &
            + 6.526555d-17*Tgas4(i) &
            - 9.344008d4*invTgas(i) &
            - 5.423979d0)
      else
        krate(i,386) = 0d0
      end if
    end do

    !O + CH2O -> O2 + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,387) = krate(i,134)*exp(1.377939d-2*(lnTgas(i)-1d0) &
            + 3.929933d-3*Tgas(i) &
            - 4.326833d-6*Tgas2(i) &
            + 2.31212d-9*Tgas3(i) &
            - 5.056932d-13*Tgas4(i) &
            - 3.455922d4*invTgas(i) &
            + 2.562099d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,387) = krate(i,134)*exp(1.082814d0*(lnTgas(i)-1d0) &
            - 1.306792d-3*Tgas(i) &
            + 2.16156d-7*Tgas2(i) &
            - 1.97315d-11*Tgas3(i) &
            + 7.416368d-16*Tgas4(i) &
            - 3.461074d4*invTgas(i) &
            - 3.488705d0)
      else
        krate(i,387) = 0d0
      end if
    end do

    !H + CH2O -> OH + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,388) = krate(i,135)*exp(8.915744d-1*(lnTgas(i)-1d0) &
            + 2.588107d-3*Tgas(i) &
            - 4.091432d-6*Tgas2(i) &
            + 2.284958d-9*Tgas3(i) &
            - 4.940869d-13*Tgas4(i) &
            - 4.264066d4*invTgas(i) &
            - 1.006848d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,388) = krate(i,135)*exp(3.040205d-1*(lnTgas(i)-1d0) &
            - 1.094927d-3*Tgas(i) &
            + 1.899826d-7*Tgas2(i) &
            - 1.752776d-11*Tgas3(i) &
            + 6.614709d-16*Tgas4(i) &
            - 4.327688d4*invTgas(i) &
            + 4.309857d0)
      else
        krate(i,388) = 0d0
      end if
    end do

    !CO + CH3 -> CHO + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,389) = krate(i,136)*exp(1.194146d0*(lnTgas(i)-1d0) &
            - 3.584024d-3*Tgas(i) &
            + 2.614098d-6*Tgas2(i) &
            - 1.19523d-9*Tgas3(i) &
            + 2.372606d-13*Tgas4(i) &
            - 5.216001d4*invTgas(i) &
            - 2.62033d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,389) = krate(i,136)*exp(1.028426d0*(lnTgas(i)-1d0) &
            - 8.654238d-4*Tgas(i) &
            + 1.622836d-7*Tgas2(i) &
            - 1.391357d-11*Tgas3(i) &
            + 4.40744d-16*Tgas4(i) &
            - 5.191408d4*invTgas(i) &
            - 3.098501d0)
      else
        krate(i,389) = 0d0
      end if
    end do

    !CH2O + CH3O -> CH3O2 + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,390) = krate(i,137)*exp(6.594806d-1*(lnTgas(i)-1d0) &
            + 2.544174d-3*Tgas(i) &
            - 3.181153d-6*Tgas2(i) &
            + 1.741609d-9*Tgas3(i) &
            - 3.840889d-13*Tgas4(i) &
            - 6.330867d4*invTgas(i) &
            - 5.106929d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,390) = krate(i,137)*exp(7.630027d-1*(lnTgas(i)-1d0) &
            - 8.081649d-4*Tgas(i) &
            + 1.487385d-7*Tgas2(i) &
            - 1.431497d-11*Tgas3(i) &
            + 5.561306d-16*Tgas4(i) &
            - 6.36269d4*invTgas(i) &
            - 4.006558d0)
      else
        krate(i,390) = 0d0
      end if
    end do

    !CH2O + CO -> CO2 + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,391) = krate(i,138)*exp(-1.82313d0*(lnTgas(i)-1d0) &
            + 8.585882d-3*Tgas(i) &
            - 6.217352d-6*Tgas2(i) &
            + 2.737413d-9*Tgas3(i) &
            - 5.241697d-13*Tgas4(i) &
            - 3.071753d4*invTgas(i) &
            + 5.042962d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,391) = krate(i,138)*exp(1.553515d0*(lnTgas(i)-1d0) &
            - 9.537688d-4*Tgas(i) &
            + 1.539653d-7*Tgas2(i) &
            - 1.423916d-11*Tgas3(i) &
            + 5.594203d-16*Tgas4(i) &
            - 3.029395d4*invTgas(i) &
            - 9.933766d0)
      else
        krate(i,391) = 0d0
      end if
    end do

    !H + CO -> O + CH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,392) = krate(i,139)*exp(5.784919d-1*(lnTgas(i)-1d0) &
            - 1.172322d-3*Tgas(i) &
            + 6.560458d-7*Tgas2(i) &
            - 3.226858d-10*Tgas3(i) &
            + 8.05452d-14*Tgas4(i) &
            - 8.865344d4*invTgas(i) &
            + 1.074491d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,392) = krate(i,139)*exp(-4.83912d-1*(lnTgas(i)-1d0) &
            + 1.931598d-4*Tgas(i) &
            + 3.354526d-9*Tgas2(i) &
            - 1.217429d-12*Tgas3(i) &
            + 4.355354d-17*Tgas4(i) &
            - 8.901335d4*invTgas(i) &
            + 6.757063d0)
      else
        krate(i,392) = 0d0
      end if
    end do

    !CHO + NO -> CH + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,393) = krate(i,140)*exp(-1.022356d0*(lnTgas(i)-1d0) &
            + 3.349763d-3*Tgas(i) &
            - 1.679916d-6*Tgas2(i) &
            + 4.543309d-10*Tgas3(i) &
            - 3.754174d-14*Tgas4(i) &
            - 5.983986d4*invTgas(i) &
            + 2.807317d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,393) = krate(i,140)*exp(2.249631d-1*(lnTgas(i)-1d0) &
            + 1.119774d-4*Tgas(i) &
            - 3.156981d-8*Tgas2(i) &
            + 3.475057d-12*Tgas3(i) &
            - 1.193703d-16*Tgas4(i) &
            - 5.973652d4*invTgas(i) &
            - 2.66201d0)
      else
        krate(i,393) = 0d0
      end if
    end do

    !O + CHO -> O2 + CH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,394) = krate(i,141)*exp(-1.335985d-1*(lnTgas(i)-1d0) &
            + 1.963829d-3*Tgas(i) &
            - 2.081461d-6*Tgas2(i) &
            + 1.086134d-9*Tgas3(i) &
            - 2.324637d-13*Tgas4(i) &
            - 3.660214d4*invTgas(i) &
            + 3.816777d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,394) = krate(i,141)*exp(-2.817548d-1*(lnTgas(i)-1d0) &
            - 3.687364d-5*Tgas(i) &
            + 1.209471d-8*Tgas2(i) &
            - 2.558441d-12*Tgas3(i) &
            + 1.635476d-16*Tgas4(i) &
            - 3.689946d4*invTgas(i) &
            + 2.317481d0)
      else
        krate(i,394) = 0d0
      end if
    end do

    !OH + CO -> O2 + CH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,395) = krate(i,142)*exp(-2.993031d-1*(lnTgas(i)-1d0) &
            + 1.695038d-4*Tgas(i) &
            + 4.206454d-7*Tgas2(i) &
            - 2.955247d-10*Tgas3(i) &
            + 6.893888d-14*Tgas4(i) &
            - 8.0572d4*invTgas(i) &
            + 2.337549d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,395) = krate(i,142)*exp(2.948813d-1*(lnTgas(i)-1d0) &
            - 1.870564d-5*Tgas(i) &
            + 2.952801d-8*Tgas2(i) &
            - 3.421171d-12*Tgas3(i) &
            + 1.237195d-16*Tgas4(i) &
            - 8.034721d4*invTgas(i) &
            - 1.041499d0)
      else
        krate(i,395) = 0d0
      end if
    end do

    !H + CH2O -> H2O + CH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,396) = krate(i,143)*exp(3.946704d-1*(lnTgas(i)-1d0) &
            + 4.098127d-3*Tgas(i) &
            - 5.415272d-6*Tgas2(i) &
            + 2.966953d-9*Tgas3(i) &
            - 6.405738d-13*Tgas4(i) &
            - 2.927256d4*invTgas(i) &
            + 1.07916d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,396) = krate(i,143)*exp(-4.715508d-1*(lnTgas(i)-1d0) &
            - 7.2733d-4*Tgas(i) &
            + 1.692197d-7*Tgas2(i) &
            - 1.768771d-11*Tgas3(i) &
            + 7.199248d-16*Tgas4(i) &
            - 3.018401d4*invTgas(i) &
            + 8.692337d0)
      else
        krate(i,396) = 0d0
      end if
    end do

    !H + CH2 -> H2 + CH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,397) = krate(i,144)*exp(-8.592238d-1*(lnTgas(i)-1d0) &
            + 5.317947d-3*Tgas(i) &
            - 4.887482d-6*Tgas2(i) &
            + 2.495826d-9*Tgas3(i) &
            - 5.357315d-13*Tgas4(i) &
            + 6.097065d3*invTgas(i) &
            + 3.960711d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,397) = krate(i,144)*exp(-1.812142d-1*(lnTgas(i)-1d0) &
            - 1.519837d-4*Tgas(i) &
            + 3.479819d-8*Tgas2(i) &
            - 3.239415d-12*Tgas3(i) &
            + 1.163187d-16*Tgas4(i) &
            + 5.795898d3*invTgas(i) &
            + 2.767231d0)
      else
        krate(i,397) = 0d0
      end if
    end do

    !CH3 -> H2 + CH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,398) = krate(i,145)*exp(2.17691d0*(lnTgas(i)-1d0) &
            + 3.089122d-3*Tgas(i) &
            - 4.437752d-6*Tgas2(i) &
            + 2.494846d-9*Tgas3(i) &
            - 5.624003d-13*Tgas4(i) &
            - 5.33201d4*invTgas(i) &
            + 1.093759d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,398) = krate(i,145)*exp(2.475682d0*(lnTgas(i)-1d0) &
            - 1.60294d-3*Tgas(i) &
            + 2.279502d-7*Tgas2(i) &
            - 1.938318d-11*Tgas3(i) &
            + 6.940581d-16*Tgas4(i) &
            - 5.36723d4*invTgas(i) &
            + 1.658374d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,398) = 0d0
      end if
    end do

    !CH3O + CH3 -> CH3OH + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,399) = krate(i,146)*exp(2.482839d0*(lnTgas(i)-1d0) &
            - 8.975679d-3*Tgas(i) &
            + 5.706182d-6*Tgas2(i) &
            - 2.378487d-9*Tgas3(i) &
            + 4.425269d-13*Tgas4(i) &
            - 7.02381d3*invTgas(i) &
            - 9.890009d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,399) = krate(i,146)*exp(-1.073628d0*(lnTgas(i)-1d0) &
            - 1.272959d-5*Tgas(i) &
            + 3.783904d-8*Tgas2(i) &
            - 4.530641d-12*Tgas3(i) &
            + 1.845117d-16*Tgas4(i) &
            - 7.601515d3*invTgas(i) &
            + 6.472214d0)
      else
        krate(i,399) = 0d0
      end if
    end do

    !CH2OH + CH3 -> CH3OH + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,400) = krate(i,147)*exp(1.716321d0*(lnTgas(i)-1d0) &
            - 9.702647d-3*Tgas(i) &
            + 7.340637d-6*Tgas2(i) &
            - 3.280193d-9*Tgas3(i) &
            + 6.359351d-13*Tgas4(i) &
            - 1.18563d4*invTgas(i) &
            - 6.62672d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,400) = krate(i,147)*exp(-1.408956d0*(lnTgas(i)-1d0) &
            + 7.3419d-4*Tgas(i) &
            - 6.750871d-8*Tgas2(i) &
            + 5.059676d-12*Tgas3(i) &
            - 1.925485d-16*Tgas4(i) &
            - 1.204979d4*invTgas(i) &
            + 6.352319d0)
      else
        krate(i,400) = 0d0
      end if
    end do

    !CH3O + OH -> CH3OH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,401) = krate(i,148)*exp(1.122988d0*(lnTgas(i)-1d0) &
            - 7.185981d-3*Tgas(i) &
            + 5.59419d-6*Tgas2(i) &
            - 2.564916d-9*Tgas3(i) &
            + 5.066689d-13*Tgas4(i) &
            + 1.166338d3*invTgas(i) &
            - 5.313807d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,401) = krate(i,148)*exp(-1.525418d0*(lnTgas(i)-1d0) &
            + 8.708622d-4*Tgas(i) &
            - 1.070113d-7*Tgas2(i) &
            + 8.520205d-12*Tgas3(i) &
            - 2.960604d-16*Tgas4(i) &
            + 8.648186d2*invTgas(i) &
            + 6.211735d0)
      else
        krate(i,401) = 0d0
      end if
    end do

    !CH2OH + OH -> CH3OH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,402) = krate(i,149)*exp(3.564702d-1*(lnTgas(i)-1d0) &
            - 7.912949d-3*Tgas(i) &
            + 7.228645d-6*Tgas2(i) &
            - 3.466621d-9*Tgas3(i) &
            + 7.000771d-13*Tgas4(i) &
            - 3.666155d3*invTgas(i) &
            - 2.050518d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,402) = krate(i,149)*exp(-1.860746d0*(lnTgas(i)-1d0) &
            + 1.617782d-3*Tgas(i) &
            - 2.12359d-7*Tgas2(i) &
            + 1.811052d-11*Tgas3(i) &
            - 6.731206d-16*Tgas4(i) &
            - 3.583453d3*invTgas(i) &
            + 6.091841d0)
      else
        krate(i,402) = 0d0
      end if
    end do

    !CH3 + H2O -> CH3OH + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,403) = krate(i,150)*exp(3.026956d-1*(lnTgas(i)-1d0) &
            - 8.194269d-3*Tgas(i) &
            + 9.535848d-6*Tgas2(i) &
            - 5.310939d-9*Tgas3(i) &
            + 1.190254d-12*Tgas4(i) &
            - 1.37327d4*invTgas(i) &
            - 2.16854d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,403) = krate(i,150)*exp(3.721085d-1*(lnTgas(i)-1d0) &
            + 7.734223d-4*Tgas(i) &
            - 1.465968d-7*Tgas2(i) &
            + 1.464305d-11*Tgas3(i) &
            - 6.015924d-16*Tgas4(i) &
            - 1.284716d4*invTgas(i) &
            - 6.884126d0)
      else
        krate(i,403) = 0d0
      end if
    end do

    !CH3O + H2 -> CH3OH + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,404) = krate(i,151)*exp(2.102374d0*(lnTgas(i)-1d0) &
            - 1.073711d-2*Tgas(i) &
            + 8.502812d-6*Tgas2(i) &
            - 4.057274d-9*Tgas3(i) &
            + 8.380015d-13*Tgas4(i) &
            + 5.281033d2*invTgas(i) &
            - 8.599432d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,404) = krate(i,151)*exp(-1.66339d0*(lnTgas(i)-1d0) &
            + 1.024923d-3*Tgas(i) &
            - 1.309126d-7*Tgas2(i) &
            + 1.032896d-11*Tgas3(i) &
            - 3.587874d-16*Tgas4(i) &
            + 1.06297d2*invTgas(i) &
            + 7.712033d0)
      else
        krate(i,404) = 0d0
      end if
    end do

    !CH2OH + H2 -> CH3OH + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,405) = krate(i,152)*exp(1.335856d0*(lnTgas(i)-1d0) &
            - 1.146408d-2*Tgas(i) &
            + 1.013727d-5*Tgas2(i) &
            - 4.95898d-9*Tgas3(i) &
            + 1.03141d-12*Tgas4(i) &
            - 4.304389d3*invTgas(i) &
            - 5.336143d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,405) = krate(i,152)*exp(-1.998718d0*(lnTgas(i)-1d0) &
            + 1.771842d-3*Tgas(i) &
            - 2.362603d-7*Tgas2(i) &
            + 1.991927d-11*Tgas3(i) &
            - 7.358477d-16*Tgas4(i) &
            - 4.341974d3*invTgas(i) &
            + 7.592139d0)
      else
        krate(i,405) = 0d0
      end if
    end do

    !CH3O + H2O -> CH3OH + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,406) = krate(i,153)*exp(1.740055d0*(lnTgas(i)-1d0) &
            - 6.929187d-3*Tgas(i) &
            + 4.93917d-6*Tgas2(i) &
            - 2.243443d-9*Tgas3(i) &
            + 4.48757d-13*Tgas4(i) &
            - 6.742926d3*invTgas(i) &
            - 6.724729d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,406) = krate(i,153)*exp(-1.069033d0*(lnTgas(i)-1d0) &
            + 5.053424d-4*Tgas(i) &
            - 7.53515d-8*Tgas2(i) &
            + 7.249489d-12*Tgas3(i) &
            - 3.009227d-16*Tgas4(i) &
            - 7.19068d3*invTgas(i) &
            + 6.096784d0)
      else
        krate(i,406) = 0d0
      end if
    end do

    !CH2OH + H2O -> CH3OH + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,407) = krate(i,154)*exp(9.735364d-1*(lnTgas(i)-1d0) &
            - 7.656155d-3*Tgas(i) &
            + 6.573624d-6*Tgas2(i) &
            - 3.145149d-9*Tgas3(i) &
            + 6.421652d-13*Tgas4(i) &
            - 1.157542d4*invTgas(i) &
            - 3.46144d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,407) = krate(i,154)*exp(-1.404361d0*(lnTgas(i)-1d0) &
            + 1.252262d-3*Tgas(i) &
            - 1.806993d-7*Tgas2(i) &
            + 1.683981d-11*Tgas3(i) &
            - 6.779829d-16*Tgas4(i) &
            - 1.163895d4*invTgas(i) &
            + 5.97689d0)
      else
        krate(i,407) = 0d0
      end if
    end do

    !CH2O + H2O + H -> CH3OH + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,408) = krate(i,155)*exp(-1.841864d0*(lnTgas(i)-1d0) &
            - 3.377337d-3*Tgas(i) &
            + 4.994686d-6*Tgas2(i) &
            - 3.025d-9*Tgas3(i) &
            + 7.228359d-13*Tgas4(i) &
            + 3.043814d3*invTgas(i) &
            - 3.084352d-1)*(1.3806488d-22*Tgas(i))**(1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,408) = krate(i,155)*exp(-1.980767d0*(lnTgas(i)-1d0) &
            + 1.129452d-3*Tgas(i) &
            - 1.497662d-7*Tgas2(i) &
            + 1.325906d-11*Tgas3(i) &
            - 5.178609d-16*Tgas4(i) &
            + 3.344158d3*invTgas(i) &
            - 1.465413d0)*(1.3806488d-22*Tgas(i))**(1)
      else
        krate(i,408) = 0d0
      end if
    end do

    !CH4 + CH3O -> CH3OH + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,409) = krate(i,156)*exp(4.556279d-1*(lnTgas(i)-1d0) &
            + 1.166565d-3*Tgas(i) &
            - 2.062972d-6*Tgas2(i) &
            + 1.170587d-9*Tgas3(i) &
            - 2.580055d-13*Tgas4(i) &
            + 2.516602d2*invTgas(i) &
            - 1.162976d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,409) = krate(i,156)*exp(-1.641898d-1*(lnTgas(i)-1d0) &
            - 4.64187d-4*Tgas(i) &
            + 7.940346d-8*Tgas2(i) &
            - 7.678486d-12*Tgas3(i) &
            + 3.074356d-16*Tgas4(i) &
            - 2.15704d2*invTgas(i) &
            + 3.374448d0)
      else
        krate(i,409) = 0d0
      end if
    end do

    !CH4 + CH2OH -> CH3OH + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,410) = krate(i,157)*exp(-3.108903d-1*(lnTgas(i)-1d0) &
            + 4.395969d-4*Tgas(i) &
            - 4.28518d-7*Tgas2(i) &
            + 2.688817d-10*Tgas3(i) &
            - 6.459725d-14*Tgas4(i) &
            - 4.580832d3*invTgas(i) &
            + 2.100312d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,410) = krate(i,157)*exp(-4.995178d-1*(lnTgas(i)-1d0) &
            + 2.827326d-4*Tgas(i) &
            - 2.59443d-8*Tgas2(i) &
            + 1.911831d-12*Tgas3(i) &
            - 6.962466d-17*Tgas4(i) &
            - 4.663975d3*invTgas(i) &
            + 3.254554d0)
      else
        krate(i,410) = 0d0
      end if
    end do

    !CH2O + CH3 -> CH2OH + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,411) = krate(i,158)*exp(2.207336d-1*(lnTgas(i)-1d0) &
            + 2.049993d-3*Tgas(i) &
            - 1.129209d-6*Tgas2(i) &
            + 1.191687d-10*Tgas3(i) &
            + 5.40019d-14*Tgas4(i) &
            - 4.479794d4*invTgas(i) &
            + 2.860521d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,411) = krate(i,158)*exp(2.08049d0*(lnTgas(i)-1d0) &
            - 1.573766d-3*Tgas(i) &
            + 2.240851d-7*Tgas2(i) &
            - 1.972451d-11*Tgas3(i) &
            + 7.378614d-16*Tgas4(i) &
            - 4.448509d4*invTgas(i) &
            - 8.551159d0)
      else
        krate(i,411) = 0d0
      end if
    end do

    !CH2O + OH -> CH2OH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,412) = krate(i,159)*exp(-1.139117d0*(lnTgas(i)-1d0) &
            + 3.839692d-3*Tgas(i) &
            - 1.241201d-6*Tgas2(i) &
            - 6.725993d-11*Tgas3(i) &
            + 1.181439d-13*Tgas4(i) &
            - 3.660779d4*invTgas(i) &
            + 4.862254d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,412) = krate(i,159)*exp(1.6287d0*(lnTgas(i)-1d0) &
            - 6.901746d-4*Tgas(i) &
            + 7.923472d-8*Tgas2(i) &
            - 6.673668d-12*Tgas3(i) &
            + 2.572893d-16*Tgas4(i) &
            - 3.601875d4*invTgas(i) &
            - 8.811637d0)
      else
        krate(i,412) = 0d0
      end if
    end do

    !CH3 + OH -> CH2OH + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,413) = krate(i,160)*exp(-6.708408d-1*(lnTgas(i)-1d0) &
            - 5.381141d-4*Tgas(i) &
            + 2.962224d-6*Tgas2(i) &
            - 2.16579d-9*Tgas3(i) &
            + 5.480888d-13*Tgas4(i) &
            - 2.157278d3*invTgas(i) &
            + 1.2929d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,413) = krate(i,160)*exp(1.776469d0*(lnTgas(i)-1d0) &
            - 4.788397d-4*Tgas(i) &
            + 3.410249d-8*Tgas2(i) &
            - 2.196756d-12*Tgas3(i) &
            + 7.639048d-17*Tgas4(i) &
            - 1.208207d3*invTgas(i) &
            - 1.286102d1)
      else
        krate(i,413) = 0d0
      end if
    end do

    !CH3OH -> CH2OH + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,414) = krate(i,161)*exp(1.319813d0*(lnTgas(i)-1d0) &
            + 7.473823d-3*Tgas(i) &
            - 6.890907d-6*Tgas2(i) &
            + 3.279213d-9*Tgas3(i) &
            - 6.626039d-13*Tgas4(i) &
            - 4.756087d4*invTgas(i) &
            + 3.759767d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,414) = krate(i,161)*exp(4.065852d0*(lnTgas(i)-1d0) &
            - 2.185146d-3*Tgas(i) &
            + 2.606607d-7*Tgas2(i) &
            - 2.120344d-11*Tgas3(i) &
            + 7.702879d-16*Tgas4(i) &
            - 4.741841d4*invTgas(i) &
            - 7.461176d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,414) = 0d0
      end if
    end do

    !CH2O + H2 -> CH2OH + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,415) = krate(i,162)*exp(-1.597311d-1*(lnTgas(i)-1d0) &
            + 2.885578d-4*Tgas(i) &
            + 1.66742d-6*Tgas2(i) &
            - 1.559619d-9*Tgas3(i) &
            + 4.494766d-13*Tgas4(i) &
            - 3.724602d4*invTgas(i) &
            + 1.576629d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,415) = krate(i,162)*exp(1.490728d0*(lnTgas(i)-1d0) &
            - 5.361141d-4*Tgas(i) &
            + 5.533346d-8*Tgas2(i) &
            - 4.864917d-12*Tgas3(i) &
            + 1.945623d-16*Tgas4(i) &
            - 3.677728d4*invTgas(i) &
            - 7.311339d0)
      else
        krate(i,415) = 0d0
      end if
    end do

    !CH3OH + HO2 -> CH2OH + H2O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,416) = krate(i,163)*exp(-1.166834d0*(lnTgas(i)-1d0) &
            + 9.424688d-3*Tgas(i) &
            - 7.477217d-6*Tgas2(i) &
            + 3.412554d-9*Tgas3(i) &
            - 6.727414d-13*Tgas4(i) &
            - 4.116444d3*invTgas(i) &
            + 3.763521d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,416) = krate(i,163)*exp(1.973338d0*(lnTgas(i)-1d0) &
            - 1.099105d-3*Tgas(i) &
            + 1.019657d-7*Tgas2(i) &
            - 6.307972d-12*Tgas3(i) &
            + 1.916311d-16*Tgas4(i) &
            - 3.906553d3*invTgas(i) &
            - 9.307199d0)
      else
        krate(i,416) = 0d0
      end if
    end do

    !CH2O + H2O -> CH2OH + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,417) = krate(i,164)*exp(-5.220509d-1*(lnTgas(i)-1d0) &
            + 4.096486d-3*Tgas(i) &
            - 1.896222d-6*Tgas2(i) &
            + 2.542123d-10*Tgas3(i) &
            + 6.023202d-14*Tgas4(i) &
            - 4.451705d4*invTgas(i) &
            + 3.451332d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,417) = krate(i,164)*exp(2.085085d0*(lnTgas(i)-1d0) &
            - 1.055694d-3*Tgas(i) &
            + 1.108945d-7*Tgas2(i) &
            - 7.944383d-12*Tgas3(i) &
            + 2.52427d-16*Tgas4(i) &
            - 4.407425d4*invTgas(i) &
            - 8.926589d0)
      else
        krate(i,417) = 0d0
      end if
    end do

    !CH2O + H2O2 -> CH2OH + HO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,418) = krate(i,165)*exp(-3.287534d-1*(lnTgas(i)-1d0) &
            + 2.327953d-3*Tgas(i) &
            - 9.926283d-7*Tgas2(i) &
            - 1.319275d-11*Tgas3(i) &
            + 9.08082d-14*Tgas4(i) &
            - 2.882519d4*invTgas(i) &
            + 3.149251d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,418) = krate(i,165)*exp(1.516108d0*(lnTgas(i)-1d0) &
            - 1.208852d-3*Tgas(i) &
            + 1.896281d-7*Tgas2(i) &
            - 1.847622d-11*Tgas3(i) &
            + 7.387788d-16*Tgas4(i) &
            - 2.852875d4*invTgas(i) &
            - 5.596279d0)
      else
        krate(i,418) = 0d0
      end if
    end do

    !CH3OH + CO -> CH2OH + CHO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,419) = krate(i,166)*exp(-5.221747d-1*(lnTgas(i)-1d0) &
            + 6.118623d-3*Tgas(i) &
            - 4.726539d-6*Tgas2(i) &
            + 2.084962d-9*Tgas3(i) &
            - 3.986745d-13*Tgas4(i) &
            - 4.03037d4*invTgas(i) &
            + 4.00639d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,419) = krate(i,166)*exp(2.437382d0*(lnTgas(i)-1d0) &
            - 1.599614d-3*Tgas(i) &
            + 2.297923d-7*Tgas2(i) &
            - 1.897325d-11*Tgas3(i) &
            + 6.332925d-16*Tgas4(i) &
            - 3.98643d4*invTgas(i) &
            - 9.45082d0)
      else
        krate(i,419) = 0d0
      end if
    end do

    !CH2O + CH2O -> CH2OH + CHO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,420) = krate(i,167)*exp(-8.71577d-1*(lnTgas(i)-1d0) &
            + 7.572608d-3*Tgas(i) &
            - 5.465433d-6*Tgas2(i) &
            + 2.162192d-9*Tgas3(i) &
            - 3.594844d-13*Tgas4(i) &
            - 2.910604d4*invTgas(i) &
            + 5.411872d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,420) = krate(i,167)*exp(2.674082d0*(lnTgas(i)-1d0) &
            - 1.958016d-3*Tgas(i) &
            + 2.94193d-7*Tgas2(i) &
            - 2.527739d-11*Tgas3(i) &
            + 8.889702d-16*Tgas4(i) &
            - 2.869266d4*invTgas(i) &
            - 1.035029d1)
      else
        krate(i,420) = 0d0
      end if
    end do

    !CH2O + CH3OH -> CH2OH + CH2OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,421) = krate(i,168)*exp(-1.495587d0*(lnTgas(i)-1d0) &
            + 1.175264d-2*Tgas(i) &
            - 8.469846d-6*Tgas2(i) &
            + 3.399361d-9*Tgas3(i) &
            - 5.819332d-13*Tgas4(i) &
            - 3.294163d4*invTgas(i) &
            + 6.912772d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,421) = krate(i,168)*exp(3.489446d0*(lnTgas(i)-1d0) &
            - 2.307956d-3*Tgas(i) &
            + 2.915938d-7*Tgas2(i) &
            - 2.478419d-11*Tgas3(i) &
            + 9.304099d-16*Tgas4(i) &
            - 3.24353d4*invTgas(i) &
            - 1.490348d1)
      else
        krate(i,421) = 0d0
      end if
    end do

    !O + NO -> N + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,422) = krate(i,169)*exp(-1.10441d0*(lnTgas(i)-1d0) &
            + 2.461233d-3*Tgas(i) &
            - 1.306678d-6*Tgas2(i) &
            + 4.822772d-10*Tgas3(i) &
            - 8.3724d-14*Tgas4(i) &
            - 1.607334d4*invTgas(i) &
            + 3.519041d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,422) = krate(i,169)*exp(2.725542d-1*(lnTgas(i)-1d0) &
            - 1.662193d-4*Tgas(i) &
            + 2.88566d-8*Tgas2(i) &
            - 1.966411d-12*Tgas3(i) &
            + 5.886389d-17*Tgas4(i) &
            - 1.577035d4*invTgas(i) &
            - 3.226327d0)
      else
        krate(i,422) = 0d0
      end if
    end do

    !N2 + O -> N + NO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,423) = krate(i,170)*exp(1.932658d-2*(lnTgas(i)-1d0) &
            - 6.184507d-4*Tgas(i) &
            + 8.173734d-7*Tgas2(i) &
            - 4.706496d-10*Tgas3(i) &
            + 1.050851d-13*Tgas4(i) &
            - 3.787445d4*invTgas(i) &
            + 1.455115d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,423) = krate(i,170)*exp(1.804419d-1*(lnTgas(i)-1d0) &
            - 1.841076d-6*Tgas(i) &
            - 8.554073d-9*Tgas2(i) &
            + 1.343213d-12*Tgas3(i) &
            - 4.909747d-17*Tgas4(i) &
            - 3.775314d4*invTgas(i) &
            + 2.244325d-1)
      else
        krate(i,423) = 0d0
      end if
    end do

    !NO + OH -> H + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,424) = krate(i,171)*exp(-1.766552d0*(lnTgas(i)-1d0) &
            + 2.727759d-3*Tgas(i) &
            + 1.661445d-7*Tgas2(i) &
            - 6.046423d-10*Tgas3(i) &
            + 1.833156d-13*Tgas4(i) &
            - 1.515628d4*invTgas(i) &
            + 3.688698d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,424) = krate(i,171)*exp(1.285511d0*(lnTgas(i)-1d0) &
            - 6.301437d-5*Tgas(i) &
            - 1.749104d-8*Tgas2(i) &
            + 3.829757d-12*Tgas3(i) &
            - 2.02752d-16*Tgas4(i) &
            - 1.417092d4*invTgas(i) &
            - 1.277805d1)
      else
        krate(i,424) = 0d0
      end if
    end do

    !O2 + NO2 -> O + NO3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,425) = krate(i,172)*exp(-2.384627d0*(lnTgas(i)-1d0) &
            + 5.896556d-3*Tgas(i) &
            - 1.469131d-6*Tgas2(i) &
            - 3.43961d-10*Tgas3(i) &
            + 2.346135d-13*Tgas4(i) &
            - 3.510249d4*invTgas(i) &
            + 6.684475d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,425) = krate(i,172)*exp(1.481399d0*(lnTgas(i)-1d0) &
            - 1.391785d-4*Tgas(i) &
            - 7.404981d-9*Tgas2(i) &
            - 6.5501d-14*Tgas3(i) &
            + 3.075374d-17*Tgas4(i) &
            - 3.42554d4*invTgas(i) &
            - 1.253746d1)
      else
        krate(i,425) = 0d0
      end if
    end do

    !N2H4 + M -> NH2 + NH2 + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,426) = krate(i,173)*exp(4.549239d0*(lnTgas(i)-1d0) &
            - 1.721463d-3*Tgas(i) &
            - 4.054954d-6*Tgas2(i) &
            + 3.297731d-9*Tgas3(i) &
            - 8.612204d-13*Tgas4(i) &
            - 3.228326d4*invTgas(i) &
            - 5.933677d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,426) = krate(i,173)*exp(2.456874d-1*(lnTgas(i)-1d0) &
            - 8.9825d-4*Tgas(i) &
            + 1.38093d-7*Tgas2(i) &
            - 1.404944d-11*Tgas3(i) &
            + 7.901018d-16*Tgas4(i) &
            - 3.386481d4*invTgas(i) &
            + 1.85057d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,426) = 0d0
      end if
    end do

    !N2H3 + H2 -> N2H4 + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,427) = krate(i,174)*exp(5.691353d-1*(lnTgas(i)-1d0) &
            - 4.989333d-3*Tgas(i) &
            + 5.802848d-6*Tgas2(i) &
            - 3.354613d-9*Tgas3(i) &
            + 7.807214d-13*Tgas4(i) &
            - 1.066109d4*invTgas(i) &
            - 3.208736d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,427) = krate(i,174)*exp(4.618722d-1*(lnTgas(i)-1d0) &
            + 3.06131d-4*Tgas(i) &
            - 6.188837d-8*Tgas2(i) &
            + 5.677946d-12*Tgas3(i) &
            - 2.153649d-16*Tgas4(i) &
            - 1.024524d4*invTgas(i) &
            - 5.000986d0)
      else
        krate(i,427) = 0d0
      end if
    end do

    !NH2 + NH2 -> N2H3 + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,428) = krate(i,175)*exp(-2.462705d0*(lnTgas(i)-1d0) &
            + 2.720536d-3*Tgas(i) &
            + 1.498464d-6*Tgas2(i) &
            - 1.622886d-9*Tgas3(i) &
            + 4.493049d-13*Tgas4(i) &
            - 8.920898d3*invTgas(i) &
            + 7.566037d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,428) = krate(i,175)*exp(1.359575d0*(lnTgas(i)-1d0) &
            + 1.78815d-4*Tgas(i) &
            - 5.180425d-8*Tgas2(i) &
            + 7.087325d-12*Tgas3(i) &
            - 5.402967d-16*Tgas4(i) &
            - 7.650338d3*invTgas(i) &
            - 1.337376d1)
      else
        krate(i,428) = 0d0
      end if
    end do

    !N2 + OH -> NH + NO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,429) = krate(i,176)*exp(1.885598d-1*(lnTgas(i)-1d0) &
            - 9.01679d-4*Tgas(i) &
            + 9.06933d-7*Tgas2(i) &
            - 4.51252d-10*Tgas3(i) &
            + 9.077285d-14*Tgas4(i) &
            - 4.962915d4*invTgas(i) &
            + 1.265488d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,429) = krate(i,176)*exp(2.533321d-1*(lnTgas(i)-1d0) &
            + 8.27847d-6*Tgas(i) &
            - 1.121273d-8*Tgas2(i) &
            + 2.260558d-12*Tgas3(i) &
            - 1.253509d-16*Tgas4(i) &
            - 4.949377d4*invTgas(i) &
            + 3.930197d-1)
      else
        krate(i,429) = 0d0
      end if
    end do

    !N + OH -> NH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,430) = krate(i,177)*exp(1.692332d-1*(lnTgas(i)-1d0) &
            - 2.832283d-4*Tgas(i) &
            + 8.955956d-8*Tgas2(i) &
            + 1.939757d-11*Tgas3(i) &
            - 1.431222d-14*Tgas4(i) &
            - 1.17547d4*invTgas(i) &
            - 1.896271d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,430) = krate(i,177)*exp(7.289018d-2*(lnTgas(i)-1d0) &
            + 1.011955d-5*Tgas(i) &
            - 2.658662d-9*Tgas2(i) &
            + 9.173455d-13*Tgas3(i) &
            - 7.625348d-17*Tgas4(i) &
            - 1.174062d4*invTgas(i) &
            + 1.685872d-1)
      else
        krate(i,430) = 0d0
      end if
    end do

    !N2 + H2O -> NH2 + NO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,431) = krate(i,178)*exp(6.809386d-1*(lnTgas(i)-1d0) &
            - 2.262923d-3*Tgas(i) &
            + 1.950754d-6*Tgas2(i) &
            - 9.614172d-10*Tgas3(i) &
            + 1.999146d-13*Tgas4(i) &
            - 6.237213d4*invTgas(i) &
            + 7.167012d-2)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,431) = krate(i,178)*exp(2.237276d-1*(lnTgas(i)-1d0) &
            + 1.488827d-4*Tgas(i) &
            - 4.09064d-8*Tgas2(i) &
            + 3.813014d-12*Tgas3(i) &
            - 4.544088d-17*Tgas4(i) &
            - 6.230501d4*invTgas(i) &
            + 1.520221d0)
      else
        krate(i,431) = 0d0
      end if
    end do

    !NH + OH -> NH2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,432) = krate(i,179)*exp(-1.246873d-1*(lnTgas(i)-1d0) &
            - 1.618038d-3*Tgas(i) &
            + 1.698842d-6*Tgas2(i) &
            - 8.316375d-10*Tgas3(i) &
            + 1.670537d-13*Tgas4(i) &
            - 4.833717d3*invTgas(i) &
            + 2.171037d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,432) = krate(i,179)*exp(-4.859893d-1*(lnTgas(i)-1d0) &
            + 5.06124d-4*Tgas(i) &
            - 6.135346d-8*Tgas2(i) &
            + 2.823172d-12*Tgas3(i) &
            + 8.477232d-17*Tgas4(i) &
            - 4.755741d3*invTgas(i) &
            + 1.242153d0)
      else
        krate(i,432) = 0d0
      end if
    end do

    !NH2 + OH -> NH3 + O(1D)
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,433) = krate(i,180)*exp(-1.223275d0*(lnTgas(i)-1d0) &
            - 6.202754d-4*Tgas(i) &
            + 1.6411d-6*Tgas2(i) &
            - 9.296993d-10*Tgas3(i) &
            + 1.98787d-13*Tgas4(i) &
            - 2.073377d4*invTgas(i) &
            + 3.466815d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,433) = krate(i,180)*exp(-8.418093d-1*(lnTgas(i)-1d0) &
            + 7.885067d-4*Tgas(i) &
            - 1.061022d-7*Tgas2(i) &
            + 9.390193d-12*Tgas3(i) &
            - 4.574504d-16*Tgas4(i) &
            - 2.041764d4*invTgas(i) &
            + 4.956503d-1)
      else
        krate(i,433) = 0d0
      end if
    end do

    !NH2 + H2O -> NH3 + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,434) = krate(i,181)*exp(6.212039d-2*(lnTgas(i)-1d0) &
            - 2.003227d-3*Tgas(i) &
            + 2.093314d-6*Tgas2(i) &
            - 1.11893d-9*Tgas3(i) &
            + 2.465167d-13*Tgas4(i) &
            - 5.768761d3*invTgas(i) &
            - 5.090189d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,434) = krate(i,181)*exp(-3.354723d-1*(lnTgas(i)-1d0) &
            + 4.024479d-4*Tgas(i) &
            - 7.346739d-8*Tgas2(i) &
            + 8.302369d-12*Tgas3(i) &
            - 4.762104d-16*Tgas4(i) &
            - 5.700516d3*invTgas(i) &
            + 6.524838d-1)
      else
        krate(i,434) = 0d0
      end if
    end do

    !NH3 + M -> NH2 + H + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,435) = krate(i,182)*exp(2.231229d0*(lnTgas(i)-1d0) &
            + 1.820895d-3*Tgas(i) &
            - 2.410598d-6*Tgas2(i) &
            + 1.252994d-9*Tgas3(i) &
            - 2.669553d-13*Tgas4(i) &
            - 5.336752d4*invTgas(i) &
            + 8.073466d-1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,435) = krate(i,182)*exp(2.996964d0*(lnTgas(i)-1d0) &
            - 1.335332d-3*Tgas(i) &
            + 1.534288d-7*Tgas2(i) &
            - 1.266601d-11*Tgas3(i) &
            + 5.685154d-16*Tgas4(i) &
            - 5.335685d4*invTgas(i) &
            - 2.13677d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,435) = 0d0
      end if
    end do

    !N2O + H -> NH + NO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,436) = krate(i,183)*exp(2.954399d0*(lnTgas(i)-1d0) &
            - 7.816407d-3*Tgas(i) &
            + 3.871093d-6*Tgas2(i) &
            - 1.378405d-9*Tgas3(i) &
            + 2.350279d-13*Tgas4(i) &
            - 1.773564d4*invTgas(i) &
            - 6.182349d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,436) = krate(i,183)*exp(-1.278634d0*(lnTgas(i)-1d0) &
            - 5.307744d-5*Tgas(i) &
            + 1.743342d-8*Tgas2(i) &
            - 1.017461d-12*Tgas3(i) &
            + 1.188801d-17*Tgas4(i) &
            - 1.872056d4*invTgas(i) &
            + 1.475826d1)
      else
        krate(i,436) = 0d0
      end if
    end do

    !NO + H -> NH + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,437) = krate(i,184)*exp(-5.738149d-2*(lnTgas(i)-1d0) &
            + 8.361791d-4*Tgas(i) &
            - 9.817179d-7*Tgas2(i) &
            + 4.745136d-10*Tgas3(i) &
            - 8.642991d-14*Tgas4(i) &
            - 3.590947d4*invTgas(i) &
            + 2.066356d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,437) = krate(i,184)*exp(-4.333489d-1*(lnTgas(i)-1d0) &
            + 5.576564d-5*Tgas(i) &
            + 2.446297d-11*Tgas2(i) &
            + 1.154676d-12*Tgas3(i) &
            - 9.75555d-17*Tgas4(i) &
            - 3.617712d4*invTgas(i) &
            + 4.740821d0)
      else
        krate(i,437) = 0d0
      end if
    end do

    !CH4 + SH -> CH3 + H2S
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,438) = krate(i,185)*exp(-1.055502d0*(lnTgas(i)-1d0) &
            + 5.350838d-3*Tgas(i) &
            - 3.806463d-6*Tgas2(i) &
            + 1.546805d-9*Tgas3(i) &
            - 2.66286d-13*Tgas4(i) &
            - 7.082239d3*invTgas(i) &
            + 5.820687d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,438) = krate(i,185)*exp(1.013599d0*(lnTgas(i)-1d0) &
            - 7.326382d-4*Tgas(i) &
            + 9.763184d-8*Tgas2(i) &
            - 8.052579d-12*Tgas3(i) &
            + 2.925135d-16*Tgas4(i) &
            - 6.886683d3*invTgas(i) &
            - 3.130948d0)
      else
        krate(i,438) = 0d0
      end if
    end do

    !CO + SH -> COS + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,439) = krate(i,186)*exp(-2.992212d0*(lnTgas(i)-1d0) &
            + 7.256481d-3*Tgas(i) &
            - 2.57692d-6*Tgas2(i) &
            + 3.914262d-10*Tgas3(i) &
            + 3.183218d-14*Tgas4(i) &
            - 5.781238d3*invTgas(i) &
            + 7.708101d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,439) = krate(i,186)*exp(1.794545d0*(lnTgas(i)-1d0) &
            - 2.528341d-4*Tgas(i) &
            + 1.915011d-8*Tgas2(i) &
            - 9.22769d-13*Tgas3(i) &
            + 1.413092d-17*Tgas4(i) &
            - 4.615975d3*invTgas(i) &
            - 1.639874d1)
      else
        krate(i,439) = 0d0
      end if
    end do

    !CO + S2 -> COS + S
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,440) = krate(i,187)*exp(-2.367654d0*(lnTgas(i)-1d0) &
            + 8.76811d-3*Tgas(i) &
            - 6.081598d-6*Tgas2(i) &
            + 2.854825d-9*Tgas3(i) &
            - 6.009088d-13*Tgas4(i) &
            - 1.428389d4*invTgas(i) &
            + 6.436286d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,440) = krate(i,187)*exp(1.372943d0*(lnTgas(i)-1d0) &
            - 3.238185d-4*Tgas(i) &
            + 3.704401d-8*Tgas2(i) &
            - 2.503223d-12*Tgas3(i) &
            + 5.836173d-17*Tgas4(i) &
            - 1.356605d4*invTgas(i) &
            - 1.115044d1)
      else
        krate(i,440) = 0d0
      end if
    end do

    !COS + NO -> CS + NO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,441) = krate(i,188)*exp(1.68469d0*(lnTgas(i)-1d0) &
            - 8.59614d-3*Tgas(i) &
            + 7.567425d-6*Tgas2(i) &
            - 3.987796d-9*Tgas3(i) &
            + 8.853476d-13*Tgas4(i) &
            - 4.362647d4*invTgas(i) &
            - 5.101077d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,441) = krate(i,188)*exp(1.90774d-2*(lnTgas(i)-1d0) &
            - 1.958738d-4*Tgas(i) &
            + 2.242506d-8*Tgas2(i) &
            - 1.075931d-12*Tgas3(i) &
            + 6.969717d-17*Tgas4(i) &
            - 4.356277d4*invTgas(i) &
            + 7.185453d-1)
      else
        krate(i,441) = 0d0
      end if
    end do

    !COS + H -> CO + SH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,442) = krate(i,189)*exp(2.992212d0*(lnTgas(i)-1d0) &
            - 7.256481d-3*Tgas(i) &
            + 2.57692d-6*Tgas2(i) &
            - 3.914262d-10*Tgas3(i) &
            - 3.183218d-14*Tgas4(i) &
            + 5.781238d3*invTgas(i) &
            - 7.708101d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,442) = krate(i,189)*exp(-1.794545d0*(lnTgas(i)-1d0) &
            + 2.528341d-4*Tgas(i) &
            - 1.915011d-8*Tgas2(i) &
            + 9.22769d-13*Tgas3(i) &
            - 1.413092d-17*Tgas4(i) &
            + 4.615974d3*invTgas(i) &
            + 1.639874d1)
      else
        krate(i,442) = 0d0
      end if
    end do

    !CO + S2 -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,443) = krate(i,190)*exp(-1.116324d0*(lnTgas(i)-1d0) &
            + 5.227182d-3*Tgas(i) &
            - 3.189675d-6*Tgas2(i) &
            + 1.369383d-9*Tgas3(i) &
            - 2.761957d-13*Tgas4(i) &
            - 4.183848d4*invTgas(i) &
            + 7.274798d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,443) = krate(i,190)*exp(1.611705d0*(lnTgas(i)-1d0) &
            - 2.875669d-4*Tgas(i) &
            + 1.105846d-8*Tgas2(i) &
            + 5.147485d-13*Tgas3(i) &
            - 2.145253d-17*Tgas4(i) &
            - 4.122564d4*invTgas(i) &
            - 1.259518d1)
      else
        krate(i,443) = 0d0
      end if
    end do

    !COS + S -> CS2 + O
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,444) = krate(i,191)*exp(1.251329d0*(lnTgas(i)-1d0) &
            - 3.540927d-3*Tgas(i) &
            + 2.891923d-6*Tgas2(i) &
            - 1.485442d-9*Tgas3(i) &
            + 3.247131d-13*Tgas4(i) &
            - 2.75546d4*invTgas(i) &
            - 5.708806d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,444) = krate(i,191)*exp(2.387615d-1*(lnTgas(i)-1d0) &
            + 3.625165d-5*Tgas(i) &
            - 2.598555d-8*Tgas2(i) &
            + 3.017971d-12*Tgas3(i) &
            - 7.981426d-17*Tgas4(i) &
            - 2.765959d4*invTgas(i) &
            - 1.444742d0)
      else
        krate(i,444) = 0d0
      end if
    end do

    !H2O + NH -> OH + NH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,445) = krate(i,192)*exp(4.923788d-1*(lnTgas(i)-1d0) &
            - 1.361244d-3*Tgas(i) &
            + 1.043821d-6*Tgas2(i) &
            - 5.101652d-10*Tgas3(i) &
            + 1.091418d-13*Tgas4(i) &
            - 1.274298d4*invTgas(i) &
            - 1.193818d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,445) = krate(i,192)*exp(-2.960452d-2*(lnTgas(i)-1d0) &
            + 1.406042d-4*Tgas(i) &
            - 2.969366d-8*Tgas2(i) &
            + 1.552456d-12*Tgas3(i) &
            + 7.991007d-17*Tgas4(i) &
            - 1.281124d4*invTgas(i) &
            + 1.127201d0)
      else
        krate(i,445) = 0d0
      end if
    end do

    !NH2 + N -> NH + NH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,446) = krate(i,193)*exp(2.939206d-1*(lnTgas(i)-1d0) &
            + 1.33481d-3*Tgas(i) &
            - 1.609282d-6*Tgas2(i) &
            + 8.51035d-10*Tgas3(i) &
            - 1.813659d-13*Tgas4(i) &
            - 6.920978d3*invTgas(i) &
            - 4.067308d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,446) = krate(i,193)*exp(5.588795d-1*(lnTgas(i)-1d0) &
            - 4.960044d-4*Tgas(i) &
            + 5.86948d-8*Tgas2(i) &
            - 1.905826d-12*Tgas3(i) &
            - 1.610258d-16*Tgas4(i) &
            - 6.984882d3*invTgas(i) &
            - 1.073566d0)
      else
        krate(i,446) = 0d0
      end if
    end do

    !NH3 + N -> NH2 + NH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,447) = krate(i,194)*exp(7.24179d-1*(lnTgas(i)-1d0) &
            + 1.976793d-3*Tgas(i) &
            - 2.658776d-6*Tgas2(i) &
            + 1.4598d-9*Tgas3(i) &
            - 3.187408d-13*Tgas4(i) &
            - 1.38952d4*invTgas(i) &
            - 1.09153d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.3d3) then
        krate(i,447) = krate(i,194)*exp(8.647473d-1*(lnTgas(i)-1d0) &
            - 7.578481d-4*Tgas(i) &
            + 1.024685d-7*Tgas2(i) &
            - 8.655739d-12*Tgas3(i) &
            + 3.950947d-16*Tgas4(i) &
            - 1.409561d4*invTgas(i) &
            - 5.988482d-1)
      else
        krate(i,447) = 0d0
      end if
    end do

    !NO + M -> O + N + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,448) = krate(i,195)*exp(1.449668d0*(lnTgas(i)-1d0) &
            + 6.802812d-4*Tgas(i) &
            - 7.335402d-7*Tgas2(i) &
            + 2.677074d-10*Tgas3(i) &
            - 3.464445d-14*Tgas4(i) &
            - 7.53818d4*invTgas(i) &
            + 3.965232d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,448) = krate(i,195)*exp(1.698868d0*(lnTgas(i)-1d0) &
            - 5.217185d-4*Tgas(i) &
            + 5.098478d-8*Tgas2(i) &
            - 2.85559d-12*Tgas3(i) &
            + 7.586524d-17*Tgas4(i) &
            - 7.543836d4*invTgas(i) &
            + 3.202899d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,448) = 0d0
      end if
    end do

    !NH + M -> H + N + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,449) = krate(i,196)*exp(1.50705d0*(lnTgas(i)-1d0) &
            - 1.558979d-4*Tgas(i) &
            + 2.481777d-7*Tgas2(i) &
            - 2.068062d-10*Tgas3(i) &
            + 5.178546d-14*Tgas4(i) &
            - 3.947233d4*invTgas(i) &
            + 1.898876d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,449) = krate(i,196)*exp(2.132216d0*(lnTgas(i)-1d0) &
            - 5.774841d-4*Tgas(i) &
            + 5.096031d-8*Tgas2(i) &
            - 4.010267d-12*Tgas3(i) &
            + 1.734207d-16*Tgas4(i) &
            - 3.926124d4*invTgas(i) &
            - 1.537922d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,449) = 0d0
      end if
    end do

    !N2O + O -> NO2 + N
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,450) = krate(i,197)*exp(1.018614d0*(lnTgas(i)-1d0) &
            - 4.805419d-3*Tgas(i) &
            + 3.947678d-6*Tgas2(i) &
            - 2.002445d-9*Tgas3(i) &
            + 4.326557d-13*Tgas4(i) &
            - 2.113722d4*invTgas(i) &
            - 2.304025d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,450) = krate(i,197)*exp(-6.601297d-2*(lnTgas(i)-1d0) &
            - 1.262114d-4*Tgas(i) &
            + 2.601039d-9*Tgas2(i) &
            + 1.89495d-12*Tgas3(i) &
            - 1.146105d-16*Tgas4(i) &
            - 2.115086d4*invTgas(i) &
            + 1.811619d0)
      else
        krate(i,450) = 0d0
      end if
    end do

    !O2 + M -> O + O + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,451) = krate(i,198)*exp(2.554078d0*(lnTgas(i)-1d0) &
            - 1.780952d-3*Tgas(i) &
            + 5.731377d-7*Tgas2(i) &
            - 2.145698d-10*Tgas3(i) &
            + 4.907955d-14*Tgas4(i) &
            - 5.930846d4*invTgas(i) &
            + 4.461912d-1)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,451) = krate(i,198)*exp(1.426313d0*(lnTgas(i)-1d0) &
            - 3.554992d-4*Tgas(i) &
            + 2.212817d-8*Tgas2(i) &
            - 8.891797d-13*Tgas3(i) &
            + 1.700135d-17*Tgas4(i) &
            - 5.9668d4*invTgas(i) &
            + 6.429226d0)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,451) = 0d0
      end if
    end do

    !HOCO + M -> OH + CO + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,452) = krate(i,199)*exp(4.676756d0*(lnTgas(i)-1d0) &
            - 5.867759d-3*Tgas(i) &
            + 1.02064d-6*Tgas2(i) &
            + 2.931537d-10*Tgas3(i) &
            - 1.529634d-13*Tgas4(i) &
            - 1.243612d4*invTgas(i) &
            - 7.834728d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,452) = krate(i,199)*exp(4.815028d-1*(lnTgas(i)-1d0) &
            - 7.296351d-4*Tgas(i) &
            + 9.880808d-8*Tgas2(i) &
            - 8.102431d-12*Tgas3(i) &
            + 2.902821d-16*Tgas4(i) &
            - 1.330317d4*invTgas(i) &
            + 1.407306d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,452) = 0d0
      end if
    end do

    !CO2 + OH -> HOCO + O(3P)
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,453) = krate(i,200)*exp(-2.857682d-1*(lnTgas(i)-1d0) &
            - 5.691425d-4*Tgas(i) &
            + 1.443017d-6*Tgas2(i) &
            - 9.330172d-10*Tgas3(i) &
            + 2.205194d-13*Tgas4(i) &
            - 5.071402d4*invTgas(i) &
            + 3.494167d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,453) = krate(i,200)*exp(4.74109d-1*(lnTgas(i)-1d0) &
            + 2.111259d-5*Tgas(i) &
            - 1.448921d-8*Tgas2(i) &
            + 1.720914d-12*Tgas3(i) &
            - 9.106425d-17*Tgas4(i) &
            - 5.068163d4*invTgas(i) &
            - 1.198772d0)
      else
        krate(i,453) = 0d0
      end if
    end do

    !CO2 + H2O -> HOCO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,454) = krate(i,201)*exp(3.312979d-1*(lnTgas(i)-1d0) &
            - 3.123488d-4*Tgas(i) &
            + 7.879958d-7*Tgas2(i) &
            - 6.115449d-10*Tgas3(i) &
            + 1.626076d-13*Tgas4(i) &
            - 5.862328d4*invTgas(i) &
            + 2.083246d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,454) = krate(i,201)*exp(9.304938d-1*(lnTgas(i)-1d0) &
            - 3.444072d-4*Tgas(i) &
            + 1.717059d-8*Tgas2(i) &
            + 4.501985d-13*Tgas3(i) &
            - 9.59265d-17*Tgas4(i) &
            - 5.873713d4*invTgas(i) &
            - 1.313723d0)
      else
        krate(i,454) = 0d0
      end if
    end do

    !H2O + CH2CO -> HOCO + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,455) = krate(i,202)*exp(2.20895d-1*(lnTgas(i)-1d0) &
            - 2.622424d-3*Tgas(i) &
            + 2.642832d-6*Tgas2(i) &
            - 1.414461d-9*Tgas3(i) &
            + 3.114513d-13*Tgas4(i) &
            - 3.045322d4*invTgas(i) &
            + 3.808863d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,455) = krate(i,202)*exp(-5.211934d-2*(lnTgas(i)-1d0) &
            + 1.959205d-4*Tgas(i) &
            - 5.248337d-8*Tgas2(i) &
            + 5.72463d-12*Tgas3(i) &
            - 2.494741d-16*Tgas4(i) &
            - 3.060926d4*invTgas(i) &
            + 5.938191d-1)
      else
        krate(i,455) = 0d0
      end if
    end do

    !CH4 + CO2 -> HOCO + CH3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,456) = krate(i,203)*exp(-9.531287d-1*(lnTgas(i)-1d0) &
            + 7.783403d-3*Tgas(i) &
            - 6.214146d-6*Tgas2(i) &
            + 2.802486d-9*Tgas3(i) &
            - 5.441549d-13*Tgas4(i) &
            - 5.16287d4*invTgas(i) &
            + 7.644998d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,456) = krate(i,203)*exp(1.835337d0*(lnTgas(i)-1d0) &
            - 1.313937d-3*Tgas(i) &
            + 1.719255d-7*Tgas2(i) &
            - 1.447778d-11*Tgas3(i) &
            + 5.124317d-16*Tgas4(i) &
            - 5.176215d4*invTgas(i) &
            - 4.036058d0)
      else
        krate(i,456) = 0d0
      end if
    end do

    !H2O + CO -> HOCO + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,457) = krate(i,204)*exp(-2.383407d0*(lnTgas(i)-1d0) &
            + 5.685426d-3*Tgas(i) &
            - 1.337924d-6*Tgas2(i) &
            - 1.590901d-10*Tgas3(i) &
            + 1.325248d-13*Tgas4(i) &
            - 4.670016d4*invTgas(i) &
            + 8.133056d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,457) = krate(i,204)*exp(2.179989d0*(lnTgas(i)-1d0) &
            - 2.032492d-4*Tgas(i) &
            - 1.884663d-8*Tgas2(i) &
            + 3.738794d-12*Tgas3(i) &
            - 1.979771d-16*Tgas4(i) &
            - 4.575419d4*invTgas(i) &
            - 1.555735d1)
      else
        krate(i,457) = 0d0
      end if
    end do

    !H2 + CO2 -> HOCO + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,458) = krate(i,205)*exp(6.936178d-1*(lnTgas(i)-1d0) &
            - 4.120277d-3*Tgas(i) &
            + 4.351638d-6*Tgas2(i) &
            - 2.425376d-9*Tgas3(i) &
            + 5.518521d-13*Tgas4(i) &
            - 5.135225d4*invTgas(i) &
            + 2.085422d-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,458) = krate(i,205)*exp(3.361366d-1*(lnTgas(i)-1d0) &
            + 1.751731d-4*Tgas(i) &
            - 3.839046d-8*Tgas2(i) &
            + 3.529665d-12*Tgas3(i) &
            - 1.537913d-16*Tgas4(i) &
            - 5.144015d4*invTgas(i) &
            + 3.015261d-1)
      else
        krate(i,458) = 0d0
      end if
    end do

    !COCOOH -> HOCO + CO
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,459) = krate(i,206)*exp(2.73141d0*(lnTgas(i)-1d0) &
            - 5.824448d-3*Tgas(i) &
            + 2.173293d-6*Tgas2(i) &
            - 4.655507d-10*Tgas3(i) &
            + 2.265984d-14*Tgas4(i) &
            + 4.515092d4*invTgas(i) &
            + 3.743024d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,459) = krate(i,206)*exp(-6.202894d-1*(lnTgas(i)-1d0) &
            - 3.377647d-4*Tgas(i) &
            + 5.702193d-8*Tgas2(i) &
            - 5.111457d-12*Tgas3(i) &
            + 1.943231d-16*Tgas4(i) &
            + 4.402031d4*invTgas(i) &
            + 2.054176d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,459) = 0d0
      end if
    end do

    !H2O2 + M -> OH + OH + M
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,460) = krate(i,207)*exp(3.668817d0*(lnTgas(i)-1d0) &
            - 1.977371d-3*Tgas(i) &
            - 1.401192d-6*Tgas2(i) &
            + 1.243164d-9*Tgas3(i) &
            - 3.181556d-13*Tgas4(i) &
            - 2.444454d4*invTgas(i) &
            - 3.48173d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,460) = krate(i,207)*exp(1.097288d0*(lnTgas(i)-1d0) &
            - 9.192171d-4*Tgas(i) &
            + 1.184078d-7*Tgas2(i) &
            - 9.505971d-12*Tgas3(i) &
            + 3.275541d-16*Tgas4(i) &
            - 2.540279d4*invTgas(i) &
            + 1.102492d1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,460) = 0d0
      end if
    end do

    !O(3P) + CO2 -> O(1D) + CO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,461) = krate(i,208)*exp(-6.683292d-1*(lnTgas(i)-1d0) &
            + 1.639745d-3*Tgas(i) &
            - 1.107235d-6*Tgas2(i) &
            + 5.107031d-10*Tgas3(i) &
            - 1.056415d-13*Tgas4(i) &
            - 2.287427d4*invTgas(i) &
            + 2.564912d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,461) = krate(i,208)*exp(-4.995222d-2*(lnTgas(i)-1d0) &
            + 2.053902d-5*Tgas(i) &
            - 9.749684d-10*Tgas2(i) &
            - 1.828914d-13*Tgas3(i) &
            + 1.389781d-17*Tgas4(i) &
            - 2.277262d4*invTgas(i) &
            - 2.717851d-1)
      else
        krate(i,461) = 0d0
      end if
    end do

    !O(3P) + N2 -> O(1D) + N2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,462) = krate(i,209)*exp(-6.683292d-1*(lnTgas(i)-1d0) &
            + 1.639745d-3*Tgas(i) &
            - 1.107235d-6*Tgas2(i) &
            + 5.107031d-10*Tgas3(i) &
            - 1.056415d-13*Tgas4(i) &
            - 2.287427d4*invTgas(i) &
            + 2.564912d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,462) = krate(i,209)*exp(-4.995222d-2*(lnTgas(i)-1d0) &
            + 2.053902d-5*Tgas(i) &
            - 9.749684d-10*Tgas2(i) &
            - 1.828914d-13*Tgas3(i) &
            + 1.389781d-17*Tgas4(i) &
            - 2.277262d4*invTgas(i) &
            - 2.717851d-1)
      else
        krate(i,462) = 0d0
      end if
    end do

    !O(3P) + SO2 -> O(1D) + SO2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,463) = krate(i,210)*exp(-6.683292d-1*(lnTgas(i)-1d0) &
            + 1.639745d-3*Tgas(i) &
            - 1.107235d-6*Tgas2(i) &
            + 5.107031d-10*Tgas3(i) &
            - 1.056415d-13*Tgas4(i) &
            - 2.287427d4*invTgas(i) &
            + 2.564912d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,463) = krate(i,210)*exp(-4.995222d-2*(lnTgas(i)-1d0) &
            + 2.053902d-5*Tgas(i) &
            - 9.749684d-10*Tgas2(i) &
            - 1.828914d-13*Tgas3(i) &
            + 1.389781d-17*Tgas4(i) &
            - 2.277262d4*invTgas(i) &
            - 2.717851d-1)
      else
        krate(i,463) = 0d0
      end if
    end do

    !CH3 + CH3 -> CH4 + CH2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,464) = krate(i,211)*exp(2.027211d0*(lnTgas(i)-1d0) &
            - 1.014224d-2*Tgas(i) &
            + 7.769155d-6*Tgas2(i) &
            - 3.549074d-9*Tgas3(i) &
            + 7.005324d-13*Tgas4(i) &
            - 7.27547d3*invTgas(i) &
            - 8.727032d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,464) = krate(i,211)*exp(-9.094383d-1*(lnTgas(i)-1d0) &
            + 4.514574d-4*Tgas(i) &
            - 4.156442d-8*Tgas2(i) &
            + 3.147845d-12*Tgas3(i) &
            - 1.229239d-16*Tgas4(i) &
            - 7.385811d3*invTgas(i) &
            + 3.097765d0)
      else
        krate(i,464) = 0d0
      end if
    end do

    !OH + H -> O + H2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,465) = krate(i,212)*exp(-9.79386d-1*(lnTgas(i)-1d0) &
            + 3.551134d-3*Tgas(i) &
            - 2.908621d-6*Tgas2(i) &
            + 1.492359d-9*Tgas3(i) &
            - 3.313326d-13*Tgas4(i) &
            + 6.382343d2*invTgas(i) &
            + 3.285625d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,465) = krate(i,212)*exp(1.379724d-1*(lnTgas(i)-1d0) &
            - 1.540606d-4*Tgas(i) &
            + 2.390126d-8*Tgas2(i) &
            - 1.808751d-12*Tgas3(i) &
            + 6.272702d-17*Tgas4(i) &
            + 7.585217d2*invTgas(i) &
            - 1.500298d0)
      else
        krate(i,465) = 0d0
      end if
    end do

    !H2 -> H + H
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,466) = krate(i,213)*exp(2.655669d0*(lnTgas(i)-1d0) &
            - 3.99026d-3*Tgas(i) &
            + 3.246358d-6*Tgas2(i) &
            - 1.679767d-9*Tgas3(i) &
            + 3.688059d-13*Tgas4(i) &
            - 5.186526d4*invTgas(i) &
            - 1.576376d0)*(1.3806488d-22*Tgas(i))**(-1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,466) = krate(i,213)*exp(2.067134d0*(lnTgas(i)-1d0) &
            - 4.13304d-4*Tgas(i) &
            + 2.440039d-8*Tgas2(i) &
            - 1.28417d-12*Tgas3(i) &
            + 3.444024d-17*Tgas4(i) &
            - 5.176039d4*invTgas(i) &
            + 1.309629d-1)*(1.3806488d-22*Tgas(i))**(-1)
      else
        krate(i,466) = 0d0
      end if
    end do

    !HO2 + CO2 -> HOCO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,467) = krate(i,214)*exp(1.860719d-2*(lnTgas(i)-1d0) &
            + 7.46177d-4*Tgas(i) &
            - 7.798854d-7*Tgas2(i) &
            + 4.706494d-10*Tgas3(i) &
            - 1.1938d-13*Tgas4(i) &
            - 2.36327d4*invTgas(i) &
            + 1.279249d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,467) = krate(i,214)*exp(2.576756d-1*(lnTgas(i)-1d0) &
            - 2.392807d-5*Tgas(i) &
            - 2.860292d-8*Tgas2(i) &
            + 4.906672d-12*Tgas3(i) &
            - 2.620011d-16*Tgas4(i) &
            - 2.390643d4*invTgas(i) &
            + 1.815663d-1)
      else
        krate(i,467) = 0d0
      end if
    end do

    do i=1,cellsNumber
      krate(i,459) = 1.57d13*(Tgas(i)/298)*exp(-4883/Tgas(i))
    end do

  end subroutine computeReverseRates

end module patmo_reverseRates

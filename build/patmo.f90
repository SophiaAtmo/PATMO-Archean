module patmo
contains

  !*************
  !initialize all
  subroutine patmo_init()
    use patmo_photo
    use patmo_parameters
    use patmo_utils

    !load photo metrics (i.e. binning)
    call loadPhotoMetric("xsecs/photoMetric.dat")
    !load photo cross-sections
    call loadAllPhotoXsecs()

    !init default photoflux
    photoFlux(:) = 0d0

    !init cumulative flux for entropy production
    ! (i.e. integrated with time)
    cumulativeFlux(:,:) = 0d0

    !set reactions rate to zero by default
    krate(:,:) = 0d0

    !load verbatim reactions
    call loadReactionsVerbatim()

  end subroutine patmo_init

  !**************
  !run model for a time-step
  subroutine patmo_run(dt)
    !use dvode_f90_m
    use patmo_parameters
    use patmo_commons
    use patmo_ode
    use patmo_jacobian
    use patmo_sparsity
    use patmo_rates
    use patmo_photoRates
    use patmo_reverseRates
    use patmo_utils
    implicit none
    real*8,intent(in)::dt
    real*8::atol(neqAll),rtol(neqAll)
    real*8::tstart,tend,n(neqAll)
    real*8::sumxn(photoBinsNumber),m(speciesNumber)
    !type(VODE_OPTS)::OPTIONS
    integer::istate,itask,i,j

    integer,parameter::meth=2
    integer,parameter::lwm=2*neqAll**2 + 2*neqAll &
         + (neqAll**2+10*neqAll)/2
    integer::itol,iopt,mf,lrw,liw
    integer::iwork(20+9*neqAll+LWM),neq(1)
    real*8::rwork(20+neqAll*6+3*neqAll+lwm)

    lrw = size(rwork)
    liw = size(iwork)

    iwork(:) = 0
    rwork(:) = 0d0

    atol(:) = 1d-10 !absolute tolerances (array)
    rtol(:) = 1d-4 !relative tolerances (array)

    !computes sparsity if not already done
    if(nonZeroElements==0) then
       call computeSparsity()
    end if

    itol = 4
    istate = 1
    itask = 1
    iopt = 0
    MF = 222
    call xsetf(0)

    !set solver options (DVODE_f90)
    !OPTIONS = SET_OPTS(SPARSE_J=.true., ABSERR_VECTOR=ATOL(:), &
    !     RELERR_VECTOR=RTOL(:), MXSTEP=100000, &
    !     USER_SUPPLIED_SPARSITY = .true., &
    !     MA28_RPS = .true., &
    !     USER_SUPPLIED_JACOBIAN = .false.)

    !set the sparsity structure (DVODE_f90)
    !CALL USERSETS_IAJA(iaSparsity, size(iaSparsity), &
    !     jaSparsity, size(jaSparsity))

    tstart = 0d0
    tend = dt

    !upper layer opacity is zero
    tauAll(:,cellsNumber) = 0d0
    !loop on cells
    do j=cellsNumber-1,1,-1
       sumxn(:) = 0d0
       !loop on reactions
       do i=1,photoReactionsNumber
          sumxn(:) = sumxn(:) + xsecAll(:,i) * nall(j,photoPartnerIndex(i))
       end do
       tauAll(:,j) = tauAll(:,j+1) + gridSpace(j) * sumxn(:)
    end do

    !unroll chemistry
    do i=1,speciesNumber
       n((i-1)*cellsNumber+1:(i*cellsNumber)) = nall(:,i)
    end do
    !unroll Tgas
    n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber)) &
         = TgasAll(:)

    !compute rates
    call computeRates(TgasAll(:))
    call computePhotoRates(tauAll(:,:))
    call computeReverseRates(TgasAll(:))
    call computeHescape() 
    !compute tot density
    ntotAll(:) = sum(nall(:,1:chemSpeciesNumber),2)

    !compute mean molecular mass of the whole atmosphere
    ! (averaged between layers)
    m(:) = getSpeciesMass()
    meanMolecularMass = 0d0
    do i=1,cellsNumber
       meanMolecularMass = meanMolecularMass &
            + sum(m(1:chemSpeciesNumber) &
            * nAll(i,1:chemSpeciesNumber)) &
            / ntotAll(i) / cellsNumber
    end do

    !call the solver (DVODE_f90)
    !CALL DVODE_F90(fex, &
    !     neqAll, n(:), &
    !     tstart, tend, ITASK, ISTATE, OPTIONS, &
    !     jex)

    neq(:) = neqAll

    !loop until istate=2 or istate=error
    do
       CALL DLSODES(fex, NEQ(:), n(:), tstart, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)
       print *, istate
       !recompute sparsity if required
       if(istate==-5.or.istate==-4) then
          istate = 3
          cycle
       end if
       !loop when max iteration reached
       if(istate/=-1) exit
       istate = 1
    end do

    !check output state
    if(istate/=2) then
       print *,"ERROR: istate=",istate
       stop
    end if

    !avoid negative species
    do i=1,neqAll
       n(i) = max(n(i),0d0)
    end do

    !roll chemistry
    do i=1,speciesNumber
       nall(:,i) = n((i-1)*cellsNumber+1:(i*cellsNumber))
    end do
    !roll Tgas
    TgasAll(:) = n((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))

  end subroutine patmo_run

  !****************
  !dump the histogram of the connection degree to ifile, for the
  ! cell icell, at a given time (or any other independent varibles)
  subroutine patmo_dumpWeightedDegreeHistogram(ifile,icell,time)
    use patmo_utils
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::time
    integer,intent(in)::ifile,icell
    integer::hist(speciesNumber),i

    !get histogram
    hist(:) = getDegreeHistogram(nAll(:,:),icell)
    !write histogram to file
    do i=1,speciesNumber
       write(ifile,*) time, i, hist(i)
    end do
    write(ifile,*)

  end subroutine patmo_dumpWeightedDegreeHistogram

  !**************
  subroutine patmo_printBestFluxes(icell,bestFluxesNumber)
    use patmo_commons
    use patmo_utils
    use patmo_parameters
    implicit none
    integer,intent(in)::bestFluxesNumber,icell
    integer::idx(bestFluxesNumber),i
    real*8::flux(reactionsNumber)

    !get fluxes
    flux(:) = getFlux(nAll(:,:),icell)

    idx(:) = getBestFluxIdx(icell,bestFluxesNumber)
    print *,"*************"
    do i=1,bestFluxesNumber
       print *,idx(i),trim(reactionsVerbatim(idx(i))),flux(idx(i))
    end do

  end subroutine patmo_printBestFluxes

  !**************
  !compute cumulative flux for entropy production
  subroutine patmo_computeEntropyProductionFlux(dt)
    use patmo_utils
    implicit none
    real*8,intent(in)::dt

    call computeEntropyProductionFlux(dt)

  end subroutine patmo_computeEntropyProductionFlux

  !***************
  function patmo_getEntropyProduction(timeInterval)
    use patmo_utils
    implicit none
    real*8,intent(in)::timeInterval
    real*8::patmo_getEntropyProduction

    patmo_getEntropyProduction = getEntropyProduction(timeInterval)

  end function patmo_getEntropyProduction

  !***************
  !Assume a black-body flux, with starTbb (K), starRadius (Rsun)
  ! starDistance (AU). Default is Sun at 1AU
  subroutine patmo_setFluxBB(starTbb,starRadius,starDistance)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    use patmo_photo
    implicit none
    real*8,optional,intent(in)::starTbb,starRadius,starDistance
    real*8,parameter::AU2cm=1.496d13 !AU->cm
    real*8,parameter::Rsun2cm=6.963d10 !Rsun->cm
    real*8::Tbb,rstar,dstar
    integer::i

    !default is Sun
    Tbb = 5.777d3 !K
    rstar = Rsun2cm !cm
    dstar = AU2cm !cm

    !check optional parameters
    if(present(starTbb)) Tbb = starTbb
    if(present(starRadius)) rstar = starRadius*Rsun2cm
    if(present(starTbb)) dstar = starDistance*AU2cm

    !integrate flux
    do i=1,photoBinsNumber
       photoFlux(i) = (fluxBB(energyLeft(i),Tbb)+fluxBB(energyRight(i),Tbb)) &
            * energySpan(i)/2d0
    end do

    !scale geometric flux
    photoFlux(:) =  pi*rstar**2/dstar**2 * photoFlux(:)
    open(58,file="solar_flux.txt",status="old")
    do i=1,photoBinsNumber
        read(58,*) photoFlux(i)
    end do
    close(58)

  end subroutine patmo_setFluxBB

  !***************
  !dump opacity to file fname using unitEenergy as unit for
  ! energy (eV or mircron, eV default).
  ! File format is energy,layer,opacity
  subroutine patmo_dumpOpacity(fname,unitEnergy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    character(len=*),optional,intent(in)::unitEnergy
    character(len=*),intent(in)::fname
    character(len=100)::unitE
    integer::i,j

    unitE = "eV"
    if(present(unitEnergy)) unitE = trim(unitEnergy)

    open(22,file=trim(fname),status="replace")
    if(trim(unitE)=="eV") then
       !loop on energy
       do i=1,photoBinsNumber
          !loop on cells
          do j=1,cellsNumber
             write(22,*) energyMid(i),height(j)/1d5,tauAll(i,j)
          end do
          write(22,*)
       end do
    else if(trim(unitE)=="micron") then
       !loop on energy
       do i=photoBinsNumber,1,-1
          !loop on cells
          do j=1,cellsNumber
             !h*c/E -> cm -> micron
             write(22,*) 1d4*planck_eV*clight/energyMid(i),height(j)/1d5,&
                  tauAll(i,j)
          end do
          write(22,*)
       end do
    else
       print *,"ERROR: unknown unit "//trim(unitEnergy)
       stop
    end if
    close(22)

    print *, "Opacity dumped in ", trim(fname)

  end subroutine patmo_dumpOpacity

  !***************
  !find hydrostatic equilbrium knowing the pressure at ground
  ! (pground), using dp/dz = -mu*p*g/k/T
  ! Pressure unit is defined in unitP, default dyne
  subroutine patmo_hydrostaticEquilibrium(pground,unitP)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    real*8,intent(in)::pground
    real*8::p,zold,dz,ntot,n(speciesNumber)
    integer::i
    character(len=*),optional,intent(in)::unitP
    character(len=50)::units

    !optional argument units
    if(present(unitP)) then
       units = trim(unitP)
    else
       !default
       units = "dyne"
    end if

    !convert initial pressure to dyne/cm2 if necessary
    if(trim(units)=="dyne") then
       p = pground
    elseif(trim(units)=="dyne/cm2") then
       p = pground
    elseif(trim(units)=="atm") then
       p = pground*1.013250d6
    elseif(trim(units)=="mbar") then
       p = pground*1d3
    elseif(trim(units)=="bar") then
       p = pground*1d6
    else
       !error if units unknonw
       print *,"ERROR: unknown pressure unit for hydrostatic eq ",trim(units)
       stop
    end if

    !initial conditions
    zold = 0d0
    !loop on cells
    do i=1,cellsNumber
       !difference in height
       dz = height(i)-zold
       !temp array
       n(:) = nall(i,:)
       !compute difference in pressure
       p = p - getMeanMass(n(:)) * p / kboltzmann &
            / TgasAll(i) * gravity * dz
       !total number density p=n*k*T
       ntot = p/kboltzmann/TgasAll(i)
       !resacale abundances depending on total pressure
       nall(i,1:chemSpeciesNumber) = &
            n(1:chemSpeciesNumber) &
            / sum(n(1:chemSpeciesNumber))*ntot
       !store old height
       zold = height(i)
    end do

  end subroutine patmo_hydrostaticEquilibrium

  !***************
  !dump hydrostatic pressure profile to fname.
  ! Format: h/km, p/mbar, Tgas/K
  subroutine patmo_dumpHydrostaticProfile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    integer::i
    character(len=*),intent(in)::fname
    real*8::ntot

    open(22,file=trim(fname),status="replace")
    write(22,*) "#hydrostatic equilibrium dump"
    write(22,*) "#alt/km p/mbar Tgas/K"
    !loop on cells
    do i=1,cellsNumber
       ntot = sum(nall(i,1:chemSpeciesNumber))
       write(22,*) height(i)/1d5,ntot*kboltzmann*TgasAll(i)/1d3,TgasAll(i)
    end do
    close(33)
    print *,"Hydrostatic equilibrium dumped in ",trim(fname)

  end subroutine patmo_dumpHydrostaticProfile

  !****************
  !load initial profile (density, etc...) from fname.
  ! Height in unitH, species in unitX
  subroutine patmo_loadInitialProfile(fname,unitH,unitX,defaultDensity)
    use patmo_parameters
    use patmo_commons
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=*),optional,intent(in)::unitH,unitX
    character(len=50)::units,unitsX
    character(len=50),allocatable::arow(:)
    real*8,optional,intent(in)::defaultDensity
    real*8,allocatable::x(:),rout(:)
    real*8::zold,defaultN
    integer::ios,i,idx,j,nonZero,offset
    logical::firstRow

    units = "cm"
    !optional argument units height
    if(present(unitH)) units = trim(unitH)

    unitsX = "1/cm3"
    !optional argument units chemical species
    if(present(unitX)) unitsX = trim(unitX)

    defaultN = 0d0
    !optional argument for default density
    if(present(defaultDensity)) defaultN = defaultDensity

    !read file
    print *,"reading ",trim(fname)
    open(22,file=trim(fname),status="old",iostat=ios)
    !check for file opening
    if(ios/=0) then
       print *,"ERROR: problem while opening ",trim(fname)
       stop
    end if

    !read until comment is found
    do
       read(22,*,iostat=ios) offset,nonZero
       if(ios==0) exit
    end do
    allocate(arow(offset+nonZero))
    allocate(x(nonZero))
    allocate(rout(offset))
    read(22,*) arow(:)

    !set default abundance
    nall(:,1:chemSpeciesNumber) = defaultN
    !set not chemial species to zero
    nall(:,chemSpeciesNumber+1:speciesNumber) = 0d0
    !loop on cells (file lines have to be the same number)
    do j=1,cellsNumber
       !read data+chemistry
       read(22,*,iostat=ios) rout(:),x(:)
       if(ios/=0) then
          print *,"ERROR: problem while reading ",trim(fname)
          if(j>1) print *,&
               "(could be less file lines than declared lines number)"
          stop
       end if

       !loop on data accoding to header
       do i=1,offset
          if(trim(arow(i))=="alt") then
             height(j) = rout(i)
          elseif(trim(arow(i))=="Tgas") then
             TgasAll(j) = rout(i)
          elseif(trim(arow(i))=="Dzz") then
             diffusionDzz(j) = rout(i)
          elseif(trim(arow(i))=="Kzz") then
             eddyKzz(j) = rout(i)
          elseif(trim(arow(i))=="index") then
             continue
          elseif(trim(arow(i))=="idx") then
             continue
          elseif(trim(arow(i))=="dummy") then
             continue
          else
             print *,"ERROR: unknown header element: ",trim(arow(i))
             stop
          end if
       end do

       !load species into common array
       do i=1,nonZero
          idx = getSpeciesIndex(arow(offset+i),error=.false.)
          if(idx/=-1) nall(j,idx) = x(i)
       end do

       !convert units if necessary
       if(trim(unitsX)=="ppbv") then
          nall(j,:) = nall(j,:)/sum(nall(j,1:chemSpeciesNumber))
       elseif(trim(unitsX)=="1/cm3") then
          continue
       else
          print *,&
               "ERROR: unknown chemical abundance units while reading profile",&
               trim(fname),trim(units)
          stop
       end if

    end do
    close(22)
    deallocate(arow)
    deallocate(x)
    deallocate(rout)

    !convert units if necessary
    if(trim(units)=="km") then
       height(:) = height(:)*1d5
    elseif(trim(units)=="cm") then
       continue
    else
       print *,"ERROR: unknown units while reading profile", &
            trim(fname),trim(units)
       stop
    end if

    !store inverse grid space squared, 1/dz**2, and dz
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/(height(j)-zold)**2
       gridSpace(j) = (height(j)-zold)
       zold = height(j)
    end do

  end subroutine patmo_loadInitialProfile

  !****************
  !return total mass in g/cm3
  function patmo_getTotalMass()
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    integer::icell
    real*8::patmo_getTotalMass
    real*8::m(speciesNumber)

    m(:) = getSpeciesMass()

    patmo_getTotalMass = 0d0
    do icell=1,cellsNumber
       patmo_getTotalMass = patmo_getTotalMass &
            + sum(m(1:chemSpeciesNumber) &
            * nall(icell,1:chemSpeciesNumber))
    end do

  end function patmo_getTotalMass

!***************************
function patmo_getTotalMassNuclei_H()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_H

patmo_getTotalMassNuclei_H = getTotalMassNuclei_H() 

end function

!***************************
function patmo_getTotalMassNuclei_O()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_O

patmo_getTotalMassNuclei_O = getTotalMassNuclei_O() 

end function

!***************************
function patmo_getTotalMassNuclei_C()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_C

patmo_getTotalMassNuclei_C = getTotalMassNuclei_C() 

end function

!***************************
function patmo_getTotalMassNuclei_M()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_M

patmo_getTotalMassNuclei_M = getTotalMassNuclei_M() 

end function

!***************************
function patmo_getTotalMassNuclei_N()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_N

patmo_getTotalMassNuclei_N = getTotalMassNuclei_N() 

end function

!***************************
function patmo_getTotalMassNuclei_S()
 use patmo_utils
 implicit none
 real*8::patmo_getTotalMassNuclei_S

patmo_getTotalMassNuclei_S = getTotalMassNuclei_S() 

end function



  !***************
  !set uniform grid spacing, cm
 subroutine patmo_setGridSpacing(dz)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::dz
    real*8::zold
    integer::j

    gridSpace(:) = dz
    !store inverse grid space squared, 1/dz**2, and height
    zold = 0d0
    do j=1,cellsNumber
       idh2(j) = 1d0/gridSpace(j)**2
       height(j) = zold
       zold = zold + gridSpace(j)
    end do

 end subroutine patmo_setGridSpacing

  !***************
  !set thermal diffusion
 subroutine patmo_setThermalDiffusion(alpha)
    use patmo_parameters
    implicit none
    real*8,intent(in)::alpha

    thermalDiffusionFactor = alpha

 end subroutine patmo_setThermalDiffusion

  !***************
  !set eddy Kzz coefficient of icell layer
 subroutine patmo_setEddyKzz(icell,kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz
    integer,intent(in)::icell

    eddyKzz(icell) = kzz

 end subroutine patmo_setEddyKzz

  !***************
  !set eddy Kzz, same for all layers
 subroutine patmo_setEddyKzzAll(kzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::kzz

    eddyKzz(:) = kzz

 end subroutine patmo_setEddyKzzAll

  !***************
  !set diffusion Dzz for layer icell
 subroutine patmo_setDiffusionDzz(icell,dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz
    integer,intent(in)::icell

    diffusionDzz(icell) = dzz

 end subroutine patmo_setDiffusionDzz

  !***************
  !set diffusion Dzz, same for all layers
 subroutine patmo_setDiffusionDzzAll(dzz)
    use patmo_parameters
    implicit none
    real*8,intent(in)::dzz

    diffusionDzz(:) = dzz

 end subroutine patmo_setDiffusionDzzAll

  !***************
  !append density of species idx to file number ifile
 subroutine patmo_dumpDensityToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx)
    end do
    write(ifile,*)

 end subroutine patmo_dumpDensityToFile

  !****************
  !dump all mixing rations to file (one column one species)
  ! first column is layer number
 subroutine patmo_dumpAllMixingRatioToFile(fname)
    use patmo_commons
    use patmo_parameters
    use patmo_utils
    implicit none
    character(len=*),intent(in)::fname
    character(len=500)::header
    character(len=maxNameLength)::names(speciesNumber)
    integer::i

    names(:) = getSpeciesNames()
    !prepare header (species names)
    header = "#layer"
    do i=1,chemSpeciesNumber
       header = trim(header)//" "//names(i)
    end do

    !open file to dump mixing ratios
    open(67,file=trim(fname),status="replace")
    !write file header (species names)
    write(67,*) trim(header)
    !write mixing ratios
    do i=1,cellsNumber
       write(67,'(I5,999E17.8e3)') i,nall(i,1:chemSpeciesNumber)
    end do
    close(67)

 end subroutine patmo_dumpAllMixingRatioToFile

  !***************
  !append mixing ration of species idx to file number ifile
 subroutine patmo_dumpMixingRatioToFile(ifile,time,idx)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::ifile,idx
    real*8,intent(in)::time
    integer::i,j

    do i=1,cellsNumber
       write(ifile,'(E17.8,I8,E17.8)') time, i, nall(i,idx) &
            / sum(nall(i,1:chemSpeciesNumber))
    end do
    write(ifile,*)

 end subroutine patmo_dumpMixingRatioToFile

  !****************
  !set gravity in cm/s2
 subroutine patmo_setGravity(g)
    use patmo_parameters
    implicit none
    real*8,intent(in)::g

    gravity = g

 end subroutine patmo_setGravity

  !****************
  !set chemistry of layer icell
 subroutine patmo_setChemistry(icell,n)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::n(speciesNumber)

    nall(icell,:) = n(:)

 end subroutine patmo_setChemistry

  !****************
  !set the same chemistry for all the layers
 subroutine patmo_setChemistryAll(n)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::n(speciesNumber)
    integer::icell

    do icell=1,cellsNumber
       nall(icell,:) = n(:)
    end do

 end subroutine patmo_setChemistryAll

  !**************
  !set Tgas for layer icell
 subroutine patmo_setTgas(icell,Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8,intent(in)::Tgas

    TgasAll(icell) = Tgas

 end subroutine patmo_setTgas

  !**************
  !set the same Tgas for all layers
 subroutine patmo_setTgasAll(Tgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::Tgas

    TgasAll(:) = Tgas

 end subroutine patmo_setTgasAll

  !**************
  !get density of species idx_species at layer icell
 function patmo_getDensity(icell,idx_species)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell,idx_species
    real*8::patmo_getDensity

    patmo_getDensity = nall(icell,idx_species)

 end function patmo_getDensity

  !**************
  !get Tgas of layer icell
 function patmo_getTgas(icell)
    use patmo_commons
    use patmo_parameters
    implicit none
    integer,intent(in)::icell
    real*8::patmo_getTgas

    patmo_getTgas = TgasAll(icell)

 end function patmo_getTgas

  !**************
  !dump J-Values
 subroutine patmo_dumpJValue(fname)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    integer::i

    open(22,file=trim(fname),status="replace")
    write(22,*) "altitude/km, O2->O+O, O3->O2+O(1D), O3->O2+O, OH->O+H, OH->O(1D)+H, HO2->OH+O, H2O->OH+H, H2O->H2+O, H2->H+H, N2O->N2+O(1D), NO2->NO+O, NO3->NO+O2, NO3->O+NO2, N2O5->NO2+NO3, N2O5->O+NO+NO3, HNO3->OH+NO2, HNO3->H+NO3, CH4->CH3+H, CH3OOH->CH3O+OH, CH2O->H+CHO, CH2O->H2+CO, CHO->H+CO, CO2->CO+O, H2O2->OH+OH, H2O2->H+HO2, COS->CO+S, SO->S+O, CS2->CS+S, H2S->SH+H, SO2->SO+O, SO3->SO2+O, H2SO4->SO2+OH+OH, CH3OH->CH3+OH, CH3OH->CH3O+H, S2O2->SO+SO, S2O->SO+S, N2H4->H+N2H3, NH3->H+NH2, NH3->H2+NH"
    !loop on cells
    do i=1,cellsNumber
        write(22,*) i, krate(i,215), krate(i,216), krate(i,217), krate(i,218), krate(i,219), krate(i,220), krate(i,221), krate(i,222), krate(i,223), krate(i,224), krate(i,225), krate(i,226), krate(i,227), krate(i,228), krate(i,229), krate(i,230), krate(i,231), krate(i,232), krate(i,233), krate(i,234), krate(i,235), krate(i,236), krate(i,237), krate(i,238), krate(i,239), krate(i,240), krate(i,241), krate(i,242), krate(i,243), krate(i,244), krate(i,245), krate(i,246), krate(i,247), krate(i,248), krate(i,249), krate(i,250), krate(i,251), krate(i,252), krate(i,253)
    end do
    write(22,*)
    close(22)

 end subroutine patmo_dumpJValue

  !**************
  !dump all reaction rates
 subroutine patmo_dumpAllRates(fname)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_rates  !Trieu added
    implicit none
    character(len=*),intent(in)::fname
    integer::i
    real*8::n(cellsNumber,speciesNumber)  !Trieu added

    call computeTotalDensityExcludingM(n) !Trieu added
    nall(:,patmo_idx_M) = n(:, patmo_idx_M)  !Trieu added
    
    open(22,file=trim(fname),status="replace")
    write(22,*) "altitude/km, O+O2+M->O3+M, O+O3->O2+O2, O(1D)+O3->O2+O2, O(1D)+O3->O2+O+O, O(1D)+N2->O+N2, O(1D)+O2->O+O2, OH+O3->HO2+O2, HO2+O3->OH+O2+O2, OH+HO2->H2O+O2, O(1D)+H2O->OH+OH, H2O+H->OH+H2, O(1D)+N2->N2O, O(1D)+N2O->N2+O2, O(1D)+N2O->NO+NO, O+NO2->NO+O2, NO+O3->NO2+O2, NO2+O3->NO3+O2, NO2+NO3+M->N2O5+M, NO2+OH+M->HNO3+M, HNO3+OH+M->NO3+H2O+M, HO2+NO->OH+NO2, H+O3->OH+O2, O+OH->H+O2, H+O2+M->HO2+M, O+HO2->OH+O2, H+HO2->OH+OH, H+HO2->O+H2O, H+HO2->H2+O2, CH4+OH->CH3+H2O, CH3+O2+M->CH3O2+M, CH3O2+HO2->CH3OOH+O2, CH3O2+NO->CH3O+NO2, CH3OOH+OH->CH2O+OH+H2O, CH3OOH+OH->CH3O2+H2O, CH3O+O2->CH2O+HO2, CH2O+OH->CHO+H2O, CHO+O2->CO+HO2, CO+OH+M->CO2+H+M, HO2+HO2+M->H2O2+O2+M, H2O2+OH->HO2+H2O, COS+OH->CO2+SH, COS+O->CO+SO, CS2+OH->SH+COS, CS2+O->CS+SO, CS+O2->COS+O, CS+O3->COS+O2, CS+O->CO+S, H2S+OH->H2O+SH, H2S+O->OH+SH, H2S+H->H2+SH, H2S+HO2->H2O+HSO, SH+O->H+SO, SH+O2->OH+SO, SH+O3->HSO+O2, SH+NO2->HSO+NO, SO+O3->SO2+O2, SO+O2->SO2+O, SO+OH->SO2+H, SO+NO2->SO2+NO, S+O2->SO+O, S+O3->O2+SO, S+OH->H+SO, SO2+HO2->OH+SO3, SO2+NO2->SO3+NO, SO2+O3->SO3+O2, HSO+O2->SO2+OH, HSO+O3->O2+O2+SH, HSO+NO2->NO+HSO2, HSO2+O2->HO2+SO2, HSO3+O2->HO2+SO3, SO3+H2O->H2SO4, SO2+O+M->SO3+M, SO2+OH+M->HSO3+M, CH3SCH3+OH->SO2, CH3SCH3+O->SO2, CH3SCH3+OH+M->SO2+CH4O3S+M, H2SO4->SO4, O(1D)+CH4->CH3+OH, O(1D)+CH4->CH3O+H, O(1D)+CH4->CH2O+H2, CH3O2+CH3O2->CH3O+CH3O+O2, CH3O2+CH3O2->CH3OH+CH2O+O2, O+CO+M->CO2+M, H+CO+M->CHO+M, H+CHO->H2+CO, CHO+CHO->CH2O+CO, OH+CHO->H2O+CO, O+CHO->H+CO2, O+CHO->OH+CO, H+CH2O->H2+CHO, O+CH2O->OH+CHO, O(1D)+H2->H+OH, OH+H2->H+H2O, SO+HO2->SO2+OH, SO+SO+M->S2O2+M, SO+S2O2->SO2+S2O, SO+SO->S+SO2, SO+SO3->SO2+SO2, SH+SH->S+H2S, SH+H->H2+S, SH+CH2O->H2S+CHO, S+S+M->S2+M, S+S2+M->S3+M, S+S3+M->S4+M, S2+S2+M->S4+M, S4+S4+M->S8+M, S2+M->S+S+M, S2+O->S+SO, O+CH3->CH3O, O+CH3->CH2O+H, H+CH3+M->CH4+M, O3+CH3->CH3O+O2, H2O2+CH3->CH4+HO2, OH+CH3->CH3O+H, OH+CH3->CH4+O, OH+CH3+M->CH3OH+M, HO2+CH3->CH3O+OH, HO2+CH3->CH4+O2, CHO+CH3->CH4+CO, CH3+CH3->CH4+CH2, CH3->H+CH2, CH3->H2+CH, CH3O+CH3->CH2O+CH4, CH2OH+CH3->CH2O+CH4, CH3O2+CH3->CH3O+CH3O, H2+CH3->CH4+H, O+CH2->CHO+H, O+CH2->H+H+CO, O+CH2->H2+CO, H+CH2->H2+CH, O2+CH2->H+H+CO2, O2+CH2->H2+CO2, O2+CH2->CO+H2O, O2+CH2->O+CH2O, OH+CH2->H+CH2O, CHO+CH2->CO+CH3, CH3O2+CH2->CH2O+CH3O, CO2+CH2->CH2O+CO, O+CH->H+CO, CH+NO2->CHO+NO, O2+CH->O+CHO, O2+CH->OH+CO, H2O+CH->H+CH2O, H2+CH->H+CH2, H2+CH->CH3, CH3OH+CH2->CH3O+CH3, CH3OH+CH2->CH2OH+CH3, CH3OH+O->CH3O+OH, CH3OH+O->CH2OH+OH, CH3OH+H->CH3+H2O, CH3OH+H->CH3O+H2, CH3OH+H->CH2OH+H2, CH3OH+OH->CH3O+H2O, CH3OH+OH->CH2OH+H2O, CH3OH+OH->CH2O+H2O+H, CH3OH+CH3->CH4+CH3O, CH3OH+CH3->CH4+CH2OH, CH2OH+CH2->CH2O+CH3, CH2OH+O->CH2O+OH, CH2OH+H->CH3+OH, CH2OH+H->CH3OH, CH2OH+H->CH2O+H2, CH2OH+H2O2->CH3OH+HO2, CH2OH+OH->CH2O+H2O, CH2OH+HO2->CH2O+H2O2, CH2OH+CHO->CH3OH+CO, CH2OH+CHO->CH2O+CH2O, CH2OH+CH2OH->CH2O+CH3OH, N+O2->O+NO, N+NO->N2+O, H+NO2->NO+OH, O+NO3->O2+NO2, NH2+NH2+M->N2H4+M, N2H4+H->N2H3+H2, N2H3+H->NH2+NH2, NH+NO->N2+OH, NH+O->N+OH, NH2+NO->N2+H2O, NH2+O->NH+OH, NH3+O(1D)->NH2+OH, NH3+OH->NH2+H2O, NH2+H+M->NH3+M, NH+NO->N2O+H, NH+O->NO+H, CH3+H2S->CH4+SH, COS+H->CO+SH, COS+S->CO+S2, CS+NO2->COS+NO, CO+SH->COS+H, CS2+O->CO+S2, CS2+O->COS+S, OH+NH2->H2O+NH, NH+NH->NH2+N, NH2+NH->NH3+N, O+N+M->NO+M, H+N+M->NH+M, NO2+N->N2O+O, O+O+M->O2+M, OH+CO+M->HOCO+M, HOCO+O(3P)->CO2+OH, HOCO+OH->CO2+H2O, HOCO+CH3->H2O+CH2CO, HOCO+CH3->CH4+CO2, HOCO+H->H2O+CO, HOCO+H->H2+CO2, HOCO+CO->COCOOH, OH+OH+M->H2O2+M, O(1D)+CO2->O(3P)+CO2, O(1D)+N2->O(3P)+N2, O(1D)+SO2->O(3P)+SO2, CH4+CH2->CH3+CH3, O+H2->OH+H, H+H->H2, HOCO+O2->HO2+CO2, O2->O+O, O3->O2+O(1D), O3->O2+O, OH->O+H, OH->O(1D)+H, HO2->OH+O, H2O->OH+H, H2O->H2+O, H2->H+H, N2O->N2+O(1D), NO2->NO+O, NO3->NO+O2, NO3->O+NO2, N2O5->NO2+NO3, N2O5->O+NO+NO3, HNO3->OH+NO2, HNO3->H+NO3, CH4->CH3+H, CH3OOH->CH3O+OH, CH2O->H+CHO, CH2O->H2+CO, CHO->H+CO, CO2->CO+O, H2O2->OH+OH, H2O2->H+HO2, COS->CO+S, SO->S+O, CS2->CS+S, H2S->SH+H, SO2->SO+O, SO3->SO2+O, H2SO4->SO2+OH+OH, CH3OH->CH3+OH, CH3OH->CH3O+H, S2O2->SO+SO, S2O->SO+S, N2H4->H+N2H3, NH3->H+NH2, NH3->H2+NH, O3+M->O+O2+M, O2+O2->O+O3, O2+O2->O(1D)+O3, O2+O+O->O(1D)+O3, O+N2->O(1D)+N2, O+O2->O(1D)+O2, HO2+O2->OH+O3, OH+O2+O2->HO2+O3, H2O+O2->OH+HO2, OH+OH->O(1D)+H2O, OH+H2->H2O+H, N2O->O(1D)+N2, N2+O2->O(1D)+N2O, NO+NO->O(1D)+N2O, NO+O2->O+NO2, NO2+O2->NO+O3, NO3+O2->NO2+O3, N2O5+M->NO2+NO3+M, HNO3+M->NO2+OH+M, NO3+H2O+M->HNO3+OH+M, OH+NO2->HO2+NO, OH+O2->H+O3, H+O2->O+OH, HO2+M->H+O2+M, OH+O2->O+HO2, OH+OH->H+HO2, O+H2O->H+HO2, H2+O2->H+HO2, CH3+H2O->CH4+OH, CH3O2+M->CH3+O2+M, CH3OOH+O2->CH3O2+HO2, CH3O+NO2->CH3O2+NO, CH2O+OH+H2O->CH3OOH+OH, CH3O2+H2O->CH3OOH+OH, CH2O+HO2->CH3O+O2, CHO+H2O->CH2O+OH, CO+HO2->CHO+O2, CO2+H+M->CO+OH+M, H2O2+O2+M->HO2+HO2+M, HO2+H2O->H2O2+OH, CO2+SH->COS+OH, CO+SO->COS+O, SH+COS->CS2+OH, CS+SO->CS2+O, COS+O->CS+O2, COS+O2->CS+O3, CO+S->CS+O, H2O+SH->H2S+OH, OH+SH->H2S+O, H2+SH->H2S+H, H2O+HSO->H2S+HO2, H+SO->SH+O, OH+SO->SH+O2, HSO+O2->SH+O3, HSO+NO->SH+NO2, SO2+O2->SO+O3, SO2+O->SO+O2, SO2+H->SO+OH, SO2+NO->SO+NO2, SO+O->S+O2, O2+SO->S+O3, H+SO->S+OH, OH+SO3->SO2+HO2, SO3+NO->SO2+NO2, SO3+O2->SO2+O3, SO2+OH->HSO+O2, O2+O2+SH->HSO+O3, NO+HSO2->HSO+NO2, HO2+SO2->HSO2+O2, HO2+SO3->HSO3+O2, H2SO4->SO3+H2O, SO3+M->SO2+O+M, HSO3+M->SO2+OH+M, SO2->CH3SCH3+OH, SO2->CH3SCH3+O, SO2+CH4O3S+M->CH3SCH3+OH+M, SO4->H2SO4, CH3+OH->O(1D)+CH4, CH3O+H->O(1D)+CH4, CH2O+H2->O(1D)+CH4, CH3O+CH3O+O2->CH3O2+CH3O2, CH3OH+CH2O+O2->CH3O2+CH3O2, CO2+M->O+CO+M, CHO+M->H+CO+M, H2+CO->H+CHO, CH2O+CO->CHO+CHO, H2O+CO->OH+CHO, H+CO2->O+CHO, OH+CO->O+CHO, H2+CHO->H+CH2O, OH+CHO->O+CH2O, H+OH->O(1D)+H2, H+H2O->OH+H2, SO2+OH->SO+HO2, S2O2+M->SO+SO+M, SO2+S2O->SO+S2O2, S+SO2->SO+SO, SO2+SO2->SO+SO3, S+H2S->SH+SH, H2+S->SH+H, H2S+CHO->SH+CH2O, S2+M->S+S+M, S3+M->S+S2+M, S4+M->S+S3+M, S4+M->S2+S2+M, S8+M->S4+S4+M, S+S+M->S2+M, S+SO->S2+O, CH3O->O+CH3, CH2O+H->O+CH3, CH4+M->H+CH3+M, CH3O+O2->O3+CH3, CH4+HO2->H2O2+CH3, CH3O+H->OH+CH3, CH4+O->OH+CH3, CH3OH+M->OH+CH3+M, CH3O+OH->HO2+CH3, CH4+O2->HO2+CH3, CH4+CO->CHO+CH3, CH4+CH2->CH3+CH3, H+CH2->CH3, H2+CH->CH3, CH2O+CH4->CH3O+CH3, CH2O+CH4->CH2OH+CH3, CH3O+CH3O->CH3O2+CH3, CH4+H->H2+CH3, CHO+H->O+CH2, H+H+CO->O+CH2, H2+CO->O+CH2, H2+CH->H+CH2, H+H+CO2->O2+CH2, H2+CO2->O2+CH2, CO+H2O->O2+CH2, O+CH2O->O2+CH2, H+CH2O->OH+CH2, CO+CH3->CHO+CH2, CH2O+CH3O->CH3O2+CH2, CH2O+CO->CO2+CH2, H+CO->O+CH, CHO+NO->CH+NO2, O+CHO->O2+CH, OH+CO->O2+CH, H+CH2O->H2O+CH, H+CH2->H2+CH, CH3->H2+CH, CH3O+CH3->CH3OH+CH2, CH2OH+CH3->CH3OH+CH2, CH3O+OH->CH3OH+O, CH2OH+OH->CH3OH+O, CH3+H2O->CH3OH+H, CH3O+H2->CH3OH+H, CH2OH+H2->CH3OH+H, CH3O+H2O->CH3OH+OH, CH2OH+H2O->CH3OH+OH, CH2O+H2O+H->CH3OH+OH, CH4+CH3O->CH3OH+CH3, CH4+CH2OH->CH3OH+CH3, CH2O+CH3->CH2OH+CH2, CH2O+OH->CH2OH+O, CH3+OH->CH2OH+H, CH3OH->CH2OH+H, CH2O+H2->CH2OH+H, CH3OH+HO2->CH2OH+H2O2, CH2O+H2O->CH2OH+OH, CH2O+H2O2->CH2OH+HO2, CH3OH+CO->CH2OH+CHO, CH2O+CH2O->CH2OH+CHO, CH2O+CH3OH->CH2OH+CH2OH, O+NO->N+O2, N2+O->N+NO, NO+OH->H+NO2, O2+NO2->O+NO3, N2H4+M->NH2+NH2+M, N2H3+H2->N2H4+H, NH2+NH2->N2H3+H, N2+OH->NH+NO, N+OH->NH+O, N2+H2O->NH2+NO, NH+OH->NH2+O, NH2+OH->NH3+O(1D), NH2+H2O->NH3+OH, NH3+M->NH2+H+M, N2O+H->NH+NO, NO+H->NH+O, CH4+SH->CH3+H2S, CO+SH->COS+H, CO+S2->COS+S, COS+NO->CS+NO2, COS+H->CO+SH, CO+S2->CS2+O, COS+S->CS2+O, H2O+NH->OH+NH2, NH2+N->NH+NH, NH3+N->NH2+NH, NO+M->O+N+M, NH+M->H+N+M, N2O+O->NO2+N, O2+M->O+O+M, HOCO+M->OH+CO+M, CO2+OH->HOCO+O(3P), CO2+H2O->HOCO+OH, H2O+CH2CO->HOCO+CH3, CH4+CO2->HOCO+CH3, H2O+CO->HOCO+H, H2+CO2->HOCO+H, COCOOH->HOCO+CO, H2O2+M->OH+OH+M, O(3P)+CO2->O(1D)+CO2, O(3P)+N2->O(1D)+N2, O(3P)+SO2->O(1D)+SO2, CH3+CH3->CH4+CH2, OH+H->O+H2, H2->H+H, HO2+CO2->HOCO+O2"
    !loop on cells
    do i=1,cellsNumber
        write(22,*) i, &
        krate(i,1)*nall(i,patmo_idx_O)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_M), &
        krate(i,2)*nall(i,patmo_idx_O)*nall(i,patmo_idx_O3), &
        krate(i,3)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_O3), &
        krate(i,4)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_O3), &
        krate(i,5)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_N2), &
        krate(i,6)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_O2), &
        krate(i,7)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_O3), &
        krate(i,8)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_O3), &
        krate(i,9)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_HO2), &
        krate(i,10)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_H2O), &
        krate(i,11)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_H), &
        krate(i,12)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_N2), &
        krate(i,13)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_N2O), &
        krate(i,14)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_N2O), &
        krate(i,15)*nall(i,patmo_idx_O)*nall(i,patmo_idx_NO2), &
        krate(i,16)*nall(i,patmo_idx_NO)*nall(i,patmo_idx_O3), &
        krate(i,17)*nall(i,patmo_idx_NO2)*nall(i,patmo_idx_O3), &
        krate(i,18)*nall(i,patmo_idx_NO2)*nall(i,patmo_idx_NO3)*nall(i,patmo_idx_M), &
        krate(i,19)*nall(i,patmo_idx_NO2)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_M), &
        krate(i,20)*nall(i,patmo_idx_HNO3)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_M), &
        krate(i,21)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_NO), &
        krate(i,22)*nall(i,patmo_idx_H)*nall(i,patmo_idx_O3), &
        krate(i,23)*nall(i,patmo_idx_O)*nall(i,patmo_idx_OH), &
        krate(i,24)*nall(i,patmo_idx_H)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_M), &
        krate(i,25)*nall(i,patmo_idx_O)*nall(i,patmo_idx_HO2), &
        krate(i,26)*nall(i,patmo_idx_H)*nall(i,patmo_idx_HO2), &
        krate(i,27)*nall(i,patmo_idx_H)*nall(i,patmo_idx_HO2), &
        krate(i,28)*nall(i,patmo_idx_H)*nall(i,patmo_idx_HO2), &
        krate(i,29)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_OH), &
        krate(i,30)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_M), &
        krate(i,31)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_HO2), &
        krate(i,32)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_NO), &
        krate(i,33)*nall(i,patmo_idx_CH3OOH)*nall(i,patmo_idx_OH), &
        krate(i,34)*nall(i,patmo_idx_CH3OOH)*nall(i,patmo_idx_OH), &
        krate(i,35)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_O2), &
        krate(i,36)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_OH), &
        krate(i,37)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_O2), &
        krate(i,38)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_M), &
        krate(i,39)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_M), &
        krate(i,40)*nall(i,patmo_idx_H2O2)*nall(i,patmo_idx_OH), &
        krate(i,41)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_OH), &
        krate(i,42)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_O), &
        krate(i,43)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_OH), &
        krate(i,44)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &
        krate(i,45)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O2), &
        krate(i,46)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O3), &
        krate(i,47)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_O), &
        krate(i,48)*nall(i,patmo_idx_H2S)*nall(i,patmo_idx_OH), &
        krate(i,49)*nall(i,patmo_idx_H2S)*nall(i,patmo_idx_O), &
        krate(i,50)*nall(i,patmo_idx_H2S)*nall(i,patmo_idx_H), &
        krate(i,51)*nall(i,patmo_idx_H2S)*nall(i,patmo_idx_HO2), &
        krate(i,52)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_O), &
        krate(i,53)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_O2), &
        krate(i,54)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_O3), &
        krate(i,55)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_NO2), &
        krate(i,56)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_O3), &
        krate(i,57)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_O2), &
        krate(i,58)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_OH), &
        krate(i,59)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_NO2), &
        krate(i,60)*nall(i,patmo_idx_S)*nall(i,patmo_idx_O2), &
        krate(i,61)*nall(i,patmo_idx_S)*nall(i,patmo_idx_O3), &
        krate(i,62)*nall(i,patmo_idx_S)*nall(i,patmo_idx_OH), &
        krate(i,63)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_HO2), &
        krate(i,64)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_NO2), &
        krate(i,65)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_O3), &
        krate(i,66)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_O2), &
        krate(i,67)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_O3), &
        krate(i,68)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_NO2), &
        krate(i,69)*nall(i,patmo_idx_HSO2)*nall(i,patmo_idx_O2), &
        krate(i,70)*nall(i,patmo_idx_HSO3)*nall(i,patmo_idx_O2), &
        krate(i,71)*nall(i,patmo_idx_SO3)*nall(i,patmo_idx_H2O), &
        krate(i,72)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_O)*nall(i,patmo_idx_M), &
        krate(i,73)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_M), &
        krate(i,74)*nall(i,patmo_idx_CH3SCH3)*nall(i,patmo_idx_OH), &
        krate(i,75)*nall(i,patmo_idx_CH3SCH3)*nall(i,patmo_idx_O), &
        krate(i,76)*nall(i,patmo_idx_CH3SCH3)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_M), &
        krate(i,77)*nall(i,patmo_idx_H2SO4), &
        krate(i,78)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_CH4), &
        krate(i,79)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_CH4), &
        krate(i,80)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_CH4), &
        krate(i,81)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_CH3O2), &
        krate(i,82)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_CH3O2), &
        krate(i,83)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_M), &
        krate(i,84)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_M), &
        krate(i,85)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CHO), &
        krate(i,86)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_CHO), &
        krate(i,87)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CHO), &
        krate(i,88)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CHO), &
        krate(i,89)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CHO), &
        krate(i,90)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CH2O), &
        krate(i,91)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH2O), &
        krate(i,92)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_H2), &
        krate(i,93)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_H2), &
        krate(i,94)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_HO2), &
        krate(i,95)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_M), &
        krate(i,96)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_S2O2), &
        krate(i,97)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_SO), &
        krate(i,98)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_SO3), &
        krate(i,99)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_SH), &
        krate(i,100)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_H), &
        krate(i,101)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_CH2O), &
        krate(i,102)*nall(i,patmo_idx_S)*nall(i,patmo_idx_S)*nall(i,patmo_idx_M), &
        krate(i,103)*nall(i,patmo_idx_S)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_M), &
        krate(i,104)*nall(i,patmo_idx_S)*nall(i,patmo_idx_S3)*nall(i,patmo_idx_M), &
        krate(i,105)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_M), &
        krate(i,106)*nall(i,patmo_idx_S4)*nall(i,patmo_idx_S4)*nall(i,patmo_idx_M), &
        krate(i,107)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_M), &
        krate(i,108)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_O), &
        krate(i,109)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH3), &
        krate(i,110)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH3), &
        krate(i,111)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_M), &
        krate(i,112)*nall(i,patmo_idx_O3)*nall(i,patmo_idx_CH3), &
        krate(i,113)*nall(i,patmo_idx_H2O2)*nall(i,patmo_idx_CH3), &
        krate(i,114)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CH3), &
        krate(i,115)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CH3), &
        krate(i,116)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_M), &
        krate(i,117)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_CH3), &
        krate(i,118)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_CH3), &
        krate(i,119)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_CH3), &
        krate(i,120)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_CH3), &
        krate(i,121)*nall(i,patmo_idx_CH3), &
        krate(i,122)*nall(i,patmo_idx_CH3), &
        krate(i,123)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_CH3), &
        krate(i,124)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_CH3), &
        krate(i,125)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_CH3), &
        krate(i,126)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CH3), &
        krate(i,127)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH2), &
        krate(i,128)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH2), &
        krate(i,129)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH2), &
        krate(i,130)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CH2), &
        krate(i,131)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_CH2), &
        krate(i,132)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_CH2), &
        krate(i,133)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_CH2), &
        krate(i,134)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_CH2), &
        krate(i,135)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CH2), &
        krate(i,136)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_CH2), &
        krate(i,137)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_CH2), &
        krate(i,138)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_CH2), &
        krate(i,139)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH), &
        krate(i,140)*nall(i,patmo_idx_CH)*nall(i,patmo_idx_NO2), &
        krate(i,141)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_CH), &
        krate(i,142)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_CH), &
        krate(i,143)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_CH), &
        krate(i,144)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CH), &
        krate(i,145)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CH), &
        krate(i,146)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_CH2), &
        krate(i,147)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_CH2), &
        krate(i,148)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_O), &
        krate(i,149)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_O), &
        krate(i,150)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_H), &
        krate(i,151)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_H), &
        krate(i,152)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_H), &
        krate(i,153)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_OH), &
        krate(i,154)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_OH), &
        krate(i,155)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_OH), &
        krate(i,156)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_CH3), &
        krate(i,157)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_CH3), &
        krate(i,158)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_CH2), &
        krate(i,159)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_O), &
        krate(i,160)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_H), &
        krate(i,161)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_H), &
        krate(i,162)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_H), &
        krate(i,163)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_H2O2), &
        krate(i,164)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_OH), &
        krate(i,165)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_HO2), &
        krate(i,166)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_CHO), &
        krate(i,167)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_CHO), &
        krate(i,168)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_CH2OH), &
        krate(i,169)*nall(i,patmo_idx_N)*nall(i,patmo_idx_O2), &
        krate(i,170)*nall(i,patmo_idx_N)*nall(i,patmo_idx_NO), &
        krate(i,171)*nall(i,patmo_idx_H)*nall(i,patmo_idx_NO2), &
        krate(i,172)*nall(i,patmo_idx_O)*nall(i,patmo_idx_NO3), &
        krate(i,173)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_M), &
        krate(i,174)*nall(i,patmo_idx_N2H4)*nall(i,patmo_idx_H), &
        krate(i,175)*nall(i,patmo_idx_N2H3)*nall(i,patmo_idx_H), &
        krate(i,176)*nall(i,patmo_idx_NH)*nall(i,patmo_idx_NO), &
        krate(i,177)*nall(i,patmo_idx_NH)*nall(i,patmo_idx_O), &
        krate(i,178)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_NO), &
        krate(i,179)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_O), &
        krate(i,180)*nall(i,patmo_idx_NH3)*nall(i,patmo_idx_O_1D), &
        krate(i,181)*nall(i,patmo_idx_NH3)*nall(i,patmo_idx_OH), &
        krate(i,182)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_H)*nall(i,patmo_idx_M), &
        krate(i,183)*nall(i,patmo_idx_NH)*nall(i,patmo_idx_NO), &
        krate(i,184)*nall(i,patmo_idx_NH)*nall(i,patmo_idx_O), &
        krate(i,185)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_H2S), &
        krate(i,186)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_H), &
        krate(i,187)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_S), &
        krate(i,188)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_NO2), &
        krate(i,189)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_SH), &
        krate(i,190)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &
        krate(i,191)*nall(i,patmo_idx_CS2)*nall(i,patmo_idx_O), &
        krate(i,192)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_NH2), &
        krate(i,193)*nall(i,patmo_idx_NH)*nall(i,patmo_idx_NH), &
        krate(i,194)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_NH), &
        krate(i,195)*nall(i,patmo_idx_O)*nall(i,patmo_idx_N)*nall(i,patmo_idx_M), &
        krate(i,196)*nall(i,patmo_idx_H)*nall(i,patmo_idx_N)*nall(i,patmo_idx_M), &
        krate(i,197)*nall(i,patmo_idx_NO2)*nall(i,patmo_idx_N), &
        krate(i,198)*nall(i,patmo_idx_O)*nall(i,patmo_idx_O)*nall(i,patmo_idx_M), &
        krate(i,199)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_M), &
        krate(i,200)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_O_3P), &
        krate(i,201)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_OH), &
        krate(i,202)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_CH3), &
        krate(i,203)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_CH3), &
        krate(i,204)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_H), &
        krate(i,205)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_H), &
        krate(i,206)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_CO), &
        krate(i,207)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_M), &
        krate(i,208)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_CO2), &
        krate(i,209)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_N2), &
        krate(i,210)*nall(i,patmo_idx_O_1D)*nall(i,patmo_idx_SO2), &
        krate(i,211)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_CH2), &
        krate(i,212)*nall(i,patmo_idx_O)*nall(i,patmo_idx_H2), &
        krate(i,213)*nall(i,patmo_idx_H)*nall(i,patmo_idx_H), &
        krate(i,214)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_O2), &
        krate(i,215)*nall(i,patmo_idx_O2), &
        krate(i,216)*nall(i,patmo_idx_O3), &
        krate(i,217)*nall(i,patmo_idx_O3), &
        krate(i,218)*nall(i,patmo_idx_OH), &
        krate(i,219)*nall(i,patmo_idx_OH), &
        krate(i,220)*nall(i,patmo_idx_HO2), &
        krate(i,221)*nall(i,patmo_idx_H2O), &
        krate(i,222)*nall(i,patmo_idx_H2O), &
        krate(i,223)*nall(i,patmo_idx_H2), &
        krate(i,224)*nall(i,patmo_idx_N2O), &
        krate(i,225)*nall(i,patmo_idx_NO2), &
        krate(i,226)*nall(i,patmo_idx_NO3), &
        krate(i,227)*nall(i,patmo_idx_NO3), &
        krate(i,228)*nall(i,patmo_idx_N2O5), &
        krate(i,229)*nall(i,patmo_idx_N2O5), &
        krate(i,230)*nall(i,patmo_idx_HNO3), &
        krate(i,231)*nall(i,patmo_idx_HNO3), &
        krate(i,232)*nall(i,patmo_idx_CH4), &
        krate(i,233)*nall(i,patmo_idx_CH3OOH), &
        krate(i,234)*nall(i,patmo_idx_CH2O), &
        krate(i,235)*nall(i,patmo_idx_CH2O), &
        krate(i,236)*nall(i,patmo_idx_CHO), &
        krate(i,237)*nall(i,patmo_idx_CO2), &
        krate(i,238)*nall(i,patmo_idx_H2O2), &
        krate(i,239)*nall(i,patmo_idx_H2O2), &
        krate(i,240)*nall(i,patmo_idx_COS), &
        krate(i,241)*nall(i,patmo_idx_SO), &
        krate(i,242)*nall(i,patmo_idx_CS2), &
        krate(i,243)*nall(i,patmo_idx_H2S), &
        krate(i,244)*nall(i,patmo_idx_SO2), &
        krate(i,245)*nall(i,patmo_idx_SO3), &
        krate(i,246)*nall(i,patmo_idx_H2SO4), &
        krate(i,247)*nall(i,patmo_idx_CH3OH), &
        krate(i,248)*nall(i,patmo_idx_CH3OH), &
        krate(i,249)*nall(i,patmo_idx_S2O2), &
        krate(i,250)*nall(i,patmo_idx_S2O), &
        krate(i,251)*nall(i,patmo_idx_N2H4), &
        krate(i,252)*nall(i,patmo_idx_NH3), &
        krate(i,253)*nall(i,patmo_idx_NH3), &
        krate(i,254)*nall(i,patmo_idx_O3)*nall(i,patmo_idx_M), &
        krate(i,255)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2), &
        krate(i,256)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2), &
        krate(i,257)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O)*nall(i,patmo_idx_O), &
        krate(i,258)*nall(i,patmo_idx_O)*nall(i,patmo_idx_N2), &
        krate(i,259)*nall(i,patmo_idx_O)*nall(i,patmo_idx_O2), &
        krate(i,260)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_O2), &
        krate(i,261)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2), &
        krate(i,262)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_O2), &
        krate(i,263)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_OH), &
        krate(i,264)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_H2), &
        krate(i,265)*nall(i,patmo_idx_N2O), &
        krate(i,266)*nall(i,patmo_idx_N2)*nall(i,patmo_idx_O2), &
        krate(i,267)*nall(i,patmo_idx_NO)*nall(i,patmo_idx_NO), &
        krate(i,268)*nall(i,patmo_idx_NO)*nall(i,patmo_idx_O2), &
        krate(i,269)*nall(i,patmo_idx_NO2)*nall(i,patmo_idx_O2), &
        krate(i,270)*nall(i,patmo_idx_NO3)*nall(i,patmo_idx_O2), &
        krate(i,271)*nall(i,patmo_idx_N2O5)*nall(i,patmo_idx_M), &
        krate(i,272)*nall(i,patmo_idx_HNO3)*nall(i,patmo_idx_M), &
        krate(i,273)*nall(i,patmo_idx_NO3)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_M), &
        krate(i,274)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_NO2), &
        krate(i,275)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_O2), &
        krate(i,276)*nall(i,patmo_idx_H)*nall(i,patmo_idx_O2), &
        krate(i,277)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_M), &
        krate(i,278)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_O2), &
        krate(i,279)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_OH), &
        krate(i,280)*nall(i,patmo_idx_O)*nall(i,patmo_idx_H2O), &
        krate(i,281)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_O2), &
        krate(i,282)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_H2O), &
        krate(i,283)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_M), &
        krate(i,284)*nall(i,patmo_idx_CH3OOH)*nall(i,patmo_idx_O2), &
        krate(i,285)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_NO2), &
        krate(i,286)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_H2O), &
        krate(i,287)*nall(i,patmo_idx_CH3O2)*nall(i,patmo_idx_H2O), &
        krate(i,288)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_HO2), &
        krate(i,289)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_H2O), &
        krate(i,290)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_HO2), &
        krate(i,291)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_H)*nall(i,patmo_idx_M), &
        krate(i,292)*nall(i,patmo_idx_H2O2)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_M), &
        krate(i,293)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_H2O), &
        krate(i,294)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_SH), &
        krate(i,295)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_SO), &
        krate(i,296)*nall(i,patmo_idx_SH)*nall(i,patmo_idx_COS), &
        krate(i,297)*nall(i,patmo_idx_CS)*nall(i,patmo_idx_SO), &
        krate(i,298)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_O), &
        krate(i,299)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_O2), &
        krate(i,300)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_S), &
        krate(i,301)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_SH), &
        krate(i,302)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_SH), &
        krate(i,303)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_SH), &
        krate(i,304)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_HSO), &
        krate(i,305)*nall(i,patmo_idx_H)*nall(i,patmo_idx_SO), &
        krate(i,306)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_SO), &
        krate(i,307)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_O2), &
        krate(i,308)*nall(i,patmo_idx_HSO)*nall(i,patmo_idx_NO), &
        krate(i,309)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_O2), &
        krate(i,310)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_O), &
        krate(i,311)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_H), &
        krate(i,312)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_NO), &
        krate(i,313)*nall(i,patmo_idx_SO)*nall(i,patmo_idx_O), &
        krate(i,314)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_SO), &
        krate(i,315)*nall(i,patmo_idx_H)*nall(i,patmo_idx_SO), &
        krate(i,316)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_SO3), &
        krate(i,317)*nall(i,patmo_idx_SO3)*nall(i,patmo_idx_NO), &
        krate(i,318)*nall(i,patmo_idx_SO3)*nall(i,patmo_idx_O2), &
        krate(i,319)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_OH), &
        krate(i,320)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_SH), &
        krate(i,321)*nall(i,patmo_idx_NO)*nall(i,patmo_idx_HSO2), &
        krate(i,322)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_SO2), &
        krate(i,323)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_SO3), &
        krate(i,324)*nall(i,patmo_idx_H2SO4), &
        krate(i,325)*nall(i,patmo_idx_SO3)*nall(i,patmo_idx_M), &
        krate(i,326)*nall(i,patmo_idx_HSO3)*nall(i,patmo_idx_M), &
        krate(i,327)*nall(i,patmo_idx_SO2), &
        krate(i,328)*nall(i,patmo_idx_SO2), &
        krate(i,329)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_CH4O3S)*nall(i,patmo_idx_M), &
        krate(i,330)*nall(i,patmo_idx_SO4), &
        krate(i,331)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_OH), &
        krate(i,332)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_H), &
        krate(i,333)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_H2), &
        krate(i,334)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_O2), &
        krate(i,335)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_O2), &
        krate(i,336)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_M), &
        krate(i,337)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_M), &
        krate(i,338)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CO), &
        krate(i,339)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CO), &
        krate(i,340)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_CO), &
        krate(i,341)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CO2), &
        krate(i,342)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CO), &
        krate(i,343)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CHO), &
        krate(i,344)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CHO), &
        krate(i,345)*nall(i,patmo_idx_H)*nall(i,patmo_idx_OH), &
        krate(i,346)*nall(i,patmo_idx_H)*nall(i,patmo_idx_H2O), &
        krate(i,347)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_OH), &
        krate(i,348)*nall(i,patmo_idx_S2O2)*nall(i,patmo_idx_M), &
        krate(i,349)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_S2O), &
        krate(i,350)*nall(i,patmo_idx_S)*nall(i,patmo_idx_SO2), &
        krate(i,351)*nall(i,patmo_idx_SO2)*nall(i,patmo_idx_SO2), &
        krate(i,352)*nall(i,patmo_idx_S)*nall(i,patmo_idx_H2S), &
        krate(i,353)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_S), &
        krate(i,354)*nall(i,patmo_idx_H2S)*nall(i,patmo_idx_CHO), &
        krate(i,355)*nall(i,patmo_idx_S2)*nall(i,patmo_idx_M), &
        krate(i,356)*nall(i,patmo_idx_S3)*nall(i,patmo_idx_M), &
        krate(i,357)*nall(i,patmo_idx_S4)*nall(i,patmo_idx_M), &
        krate(i,358)*nall(i,patmo_idx_S4)*nall(i,patmo_idx_M), &
        krate(i,359)*nall(i,patmo_idx_S8)*nall(i,patmo_idx_M), &
        krate(i,360)*nall(i,patmo_idx_S)*nall(i,patmo_idx_S)*nall(i,patmo_idx_M), &
        krate(i,361)*nall(i,patmo_idx_S)*nall(i,patmo_idx_SO), &
        krate(i,362)*nall(i,patmo_idx_CH3O), &
        krate(i,363)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_H), &
        krate(i,364)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_M), &
        krate(i,365)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_O2), &
        krate(i,366)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_HO2), &
        krate(i,367)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_H), &
        krate(i,368)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_O), &
        krate(i,369)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_M), &
        krate(i,370)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_OH), &
        krate(i,371)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_O2), &
        krate(i,372)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_CO), &
        krate(i,373)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_CH2), &
        krate(i,374)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CH2), &
        krate(i,375)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CH), &
        krate(i,376)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CH4), &
        krate(i,377)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CH4), &
        krate(i,378)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_CH3O), &
        krate(i,379)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_H), &
        krate(i,380)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_H), &
        krate(i,381)*nall(i,patmo_idx_H)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CO), &
        krate(i,382)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CO), &
        krate(i,383)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CH), &
        krate(i,384)*nall(i,patmo_idx_H)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CO2), &
        krate(i,385)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CO2), &
        krate(i,386)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_H2O), &
        krate(i,387)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CH2O), &
        krate(i,388)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CH2O), &
        krate(i,389)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_CH3), &
        krate(i,390)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CH3O), &
        krate(i,391)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CO), &
        krate(i,392)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CO), &
        krate(i,393)*nall(i,patmo_idx_CHO)*nall(i,patmo_idx_NO), &
        krate(i,394)*nall(i,patmo_idx_O)*nall(i,patmo_idx_CHO), &
        krate(i,395)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_CO), &
        krate(i,396)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CH2O), &
        krate(i,397)*nall(i,patmo_idx_H)*nall(i,patmo_idx_CH2), &
        krate(i,398)*nall(i,patmo_idx_CH3), &
        krate(i,399)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_CH3), &
        krate(i,400)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_CH3), &
        krate(i,401)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_OH), &
        krate(i,402)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_OH), &
        krate(i,403)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_H2O), &
        krate(i,404)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_H2), &
        krate(i,405)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_H2), &
        krate(i,406)*nall(i,patmo_idx_CH3O)*nall(i,patmo_idx_H2O), &
        krate(i,407)*nall(i,patmo_idx_CH2OH)*nall(i,patmo_idx_H2O), &
        krate(i,408)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_H), &
        krate(i,409)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_CH3O), &
        krate(i,410)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_CH2OH), &
        krate(i,411)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CH3), &
        krate(i,412)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_OH), &
        krate(i,413)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_OH), &
        krate(i,414)*nall(i,patmo_idx_CH3OH), &
        krate(i,415)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_H2), &
        krate(i,416)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_HO2), &
        krate(i,417)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_H2O), &
        krate(i,418)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_H2O2), &
        krate(i,419)*nall(i,patmo_idx_CH3OH)*nall(i,patmo_idx_CO), &
        krate(i,420)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CH2O), &
        krate(i,421)*nall(i,patmo_idx_CH2O)*nall(i,patmo_idx_CH3OH), &
        krate(i,422)*nall(i,patmo_idx_O)*nall(i,patmo_idx_NO), &
        krate(i,423)*nall(i,patmo_idx_N2)*nall(i,patmo_idx_O), &
        krate(i,424)*nall(i,patmo_idx_NO)*nall(i,patmo_idx_OH), &
        krate(i,425)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_NO2), &
        krate(i,426)*nall(i,patmo_idx_N2H4)*nall(i,patmo_idx_M), &
        krate(i,427)*nall(i,patmo_idx_N2H3)*nall(i,patmo_idx_H2), &
        krate(i,428)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_NH2), &
        krate(i,429)*nall(i,patmo_idx_N2)*nall(i,patmo_idx_OH), &
        krate(i,430)*nall(i,patmo_idx_N)*nall(i,patmo_idx_OH), &
        krate(i,431)*nall(i,patmo_idx_N2)*nall(i,patmo_idx_H2O), &
        krate(i,432)*nall(i,patmo_idx_NH)*nall(i,patmo_idx_OH), &
        krate(i,433)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_OH), &
        krate(i,434)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_H2O), &
        krate(i,435)*nall(i,patmo_idx_NH3)*nall(i,patmo_idx_M), &
        krate(i,436)*nall(i,patmo_idx_N2O)*nall(i,patmo_idx_H), &
        krate(i,437)*nall(i,patmo_idx_NO)*nall(i,patmo_idx_H), &
        krate(i,438)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_SH), &
        krate(i,439)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_SH), &
        krate(i,440)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_S2), &
        krate(i,441)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_NO), &
        krate(i,442)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_H), &
        krate(i,443)*nall(i,patmo_idx_CO)*nall(i,patmo_idx_S2), &
        krate(i,444)*nall(i,patmo_idx_COS)*nall(i,patmo_idx_S), &
        krate(i,445)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_NH), &
        krate(i,446)*nall(i,patmo_idx_NH2)*nall(i,patmo_idx_N), &
        krate(i,447)*nall(i,patmo_idx_NH3)*nall(i,patmo_idx_N), &
        krate(i,448)*nall(i,patmo_idx_NO)*nall(i,patmo_idx_M), &
        krate(i,449)*nall(i,patmo_idx_NH)*nall(i,patmo_idx_M), &
        krate(i,450)*nall(i,patmo_idx_N2O)*nall(i,patmo_idx_O), &
        krate(i,451)*nall(i,patmo_idx_O2)*nall(i,patmo_idx_M), &
        krate(i,452)*nall(i,patmo_idx_HOCO)*nall(i,patmo_idx_M), &
        krate(i,453)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_OH), &
        krate(i,454)*nall(i,patmo_idx_CO2)*nall(i,patmo_idx_H2O), &
        krate(i,455)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_CH2CO), &
        krate(i,456)*nall(i,patmo_idx_CH4)*nall(i,patmo_idx_CO2), &
        krate(i,457)*nall(i,patmo_idx_H2O)*nall(i,patmo_idx_CO), &
        krate(i,458)*nall(i,patmo_idx_H2)*nall(i,patmo_idx_CO2), &
        krate(i,459)*nall(i,patmo_idx_COCOOH), &
        krate(i,460)*nall(i,patmo_idx_H2O2)*nall(i,patmo_idx_M), &
        krate(i,461)*nall(i,patmo_idx_O_3P)*nall(i,patmo_idx_CO2), &
        krate(i,462)*nall(i,patmo_idx_O_3P)*nall(i,patmo_idx_N2), &
        krate(i,463)*nall(i,patmo_idx_O_3P)*nall(i,patmo_idx_SO2), &
        krate(i,464)*nall(i,patmo_idx_CH3)*nall(i,patmo_idx_CH3), &
        krate(i,465)*nall(i,patmo_idx_OH)*nall(i,patmo_idx_H), &
        krate(i,466)*nall(i,patmo_idx_H2), &
        krate(i,467)*nall(i,patmo_idx_HO2)*nall(i,patmo_idx_CO2)
    end do
    write(22,*)
    close(22)

 end subroutine patmo_dumpAllRates

 subroutine patmo_dumpAllNumberDensityDifference(ifile,nb,na)
  use patmo_commons
  use patmo_parameters
  implicit none
  integer,intent(in)::ifile
  real*8,intent(in)::nb(neqAll), na(cellsNumber, speciesNumber)
  real*8::deltaNAll(cellsNumber, speciesNumber)
  integer::i,j

  !compute deference
  do i = 1, speciesNumber
      deltaNAll(:, i) = nb((i - 1) * cellsNumber + 1 : (i * cellsNumber)) - na(:, i)
  end do

  ! do j = 1, cellsNumber
  !     do i = 1, speciesNumber
  !         write(ifile, *) j, i, deltaNAll(j, i)
  !     end do
  ! end do
  do i = 1, speciesNumber
      write (ifile, *) i, deltaNAll(:, i)
  end do
  write(ifile,*)

 end subroutine patmo_dumpAllNumberDensityDifference

 subroutine computeHescape()
   use patmo_constants
   use patmo_parameters
   implicit none
    
   Hesc = 2.5d8 * &
            (nAll(cellsNumber, patmo_idx_H) &
            + 2d0 * nAll(cellsNumber, patmo_idx_H2) &
            + 2d0 * nAll(cellsNumber, patmo_idx_H2O) &
            + 4d0 * nAll(cellsNumber, patmo_idx_CH4)) &
            / sum(nAll(cellsNumber,1:chemSpeciesNumber))
   H2esc = 2.5d8 * &
            (nAll(cellsNumber, patmo_idx_H2) &
            + nAll(cellsNumber, patmo_idx_H2O)&
            + 2d0 * nAll(cellsNumber, patmo_idx_CH4))&
            / sum(nAll(cellsNumber,1:chemSpeciesNumber))
 end subroutine computeHescape

 subroutine patmo_dumpHescape(ifile, time)
   use patmo_parameters
   implicit none
   integer,intent(in)::ifile
   real*8, intent(in)::time
   write (ifile, *) time, Hesc, H2esc
 end subroutine patmo_dumpHescape

end module patmo

####sulfur cycle model V8

###How to make core
1. Make folders for each network in tests folder
Necessary files
(1)NTW 
Ex. #@var:T=Tgas
     @format:idx,R,R,P,P,P,rate
    1,COS,OH,CO2,SH,,1.1d-13*exp(-1200/T)		
    2,COS,O,CO,SO,,2.1d-11*exp(-2200/T)		
    3,CS2,OH,SH,COS,,2.0d-15
(2)profile
(3)test:main program 
(4)options 
Ex.cellsNumber 	= 60          
   photoBinsNumber = 4440             #0.05nm=(400-180)/4440
   network 	= tests/sulfur/sulfur.ntw  
   energyMin 	= 3.099000            #400nm
   energyMax	= 6.888010            #180nm
   usePhotochemistry = T             
(5)copylist
(6)plot

2. Make files of xsecs in data folder
Ex.
 0 Branching ratio for O2       Total        1 branches    
 Lambda          Total           O/O
 1800.00         5.62E-20        5.62E-20
 1800.50         5.51E-20        5.51E-20

3. Make core
./patmo -test=(name of network folder)


###Change code

<patmo.f90:243> 
subroutine patmo_setFluxBB(starTbb,starRadius,starDistance)
□Read the Solar Flux directly from the file.

<patmo_photoRates.f90:23> 
function integrateXsec(index,tau)
□Change calculation formula of photodissociation rate constant.

<patmo_ode.f90>
□Set dn of fixed chemical species to 0.
□Add emission.
□Add dry deposition.
□Add wet deposition.
□Add aerosol formation.
□Add gravity settling of aerosol.

<patmo_reverseRates.f90>
□If put the reaction of the dummy, the reverse rate of the reaction with krate = 0d0.
□Reverse Intrinsic Functions (min and max) of Tgas(i)


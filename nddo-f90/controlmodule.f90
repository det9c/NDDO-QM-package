module control
! this globally accessed module holds all the available keywords
double precision::scftol,OPTTOL,mult,sys_charge,del1,ste1,cutoff,cutoff_sq,sdtocg
double precision::lindep_tol
character(20),dimension(:),allocatable::keyword
logical::RESTRICTED,DEORTHO,GRADIENT,TRIAL,keep,XYZ,OPTIMIZE,parameters,CENTERS &
,newton,DFTHF,LAPACK,SPARKLES,densityin,densityout,FREQUENCY,SAVE_TREE,TS,DIIS,MODES &
,CALCALL,IRC,EULER,RK4,FORWARD,SAM,FREQREAD,abinitio,debug,CGSCF,pbc
character(20)::core,basis_set
character(8)::method
integer::UPDATE,IRCSTPS,IRCSTPSZ,ERIBUFF

!parallel directives
logical master
integer::myrank,nprocs


end module control

 module tables

integer,dimension(:),allocatable::zeff,ntype,species,zstore
integer,dimension(:,:),allocatable::opt,ref,intadd,internal
character(5),dimension(:),allocatable::inttype
double precision,dimension(:,:),allocatable::g1,g2,g3,q,ac,bc,bmatrix
double precision,dimension(:,:),allocatable::bmat4int
double precision,dimension(:),allocatable::alpha,g,alpha3,gint,transvec,transvec2

double precision ,dimension(:), allocatable ::x,y,z,uss,upp,betas,betap,zs,zp,gss,gsp,gpp,    &
 gp2,hsp,d1,d2,p0,p1,p2,d3,d4,d5,d0,p3,p4,p5,p6,p7,zetad,udd,betad,gdd,hpp,zsone,zpone,zdone,alps,alpp,alpd &
,bsom1,bpom1,bdom1
! dipole moment factors
double precision ,dimension(:),allocatable::hyfsp,hyfpd
double precision,dimension(:),allocatable ::twoe
!vectorize make this a vector and not matrix
double precision ,dimension(:), allocatable ::S,H

 integer ,dimension(:), allocatable ::nbas,pairs

! d orbital stuff
double precision ,dimension(:), allocatable :: F0DD,F2DD,F4DD,F0SD,G2SD, &
              F0PD,F2PD,G1PD,G3PD, &
    IF0SD,IG2SD
double precision,dimension(:,:),allocatable ::REPD
double precision,dimension(30,30)::B
double precision,dimension(30)::F


! list the effective core charges
double precision:: eff_core(103)=(/ &
1.0,2.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.,2.,3.,3.,3.,3.,3.,3.,3.,&
                 3.,3.,3.,3.,1.,4.,4.,&
                 4.,4.,6.,5.,4.,4.,6.,5.,6.,4., &  ! start demp here (last entry on row) for 24 atoms
                 4.,6.,5.,4.,4.,1.,& !end of original list. put H* and O* at elements 90 and 91
                 1.,1.,6.,5.0,6.0,4.0,6.0,4.,1.,1.,1.,1.,1.,1.,1.,6.,1./) ! starts at 87 h=90 c=91 o=92


! list the quantum numbers for spd. Note that Hydrogen and Helium use 1s 2p for polarization functions
!integer:: nqs(18)=(/1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3/)
!integer:: nqp(18)=(/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3/)
!integer:: nqd(18)=(/3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/)
integer,dimension(103)::nqs,nqp,nqd
data nqs /2*1,8*2,8*3,18*4,18*5,13*6,1,2,2,2,2,2,3,  2,2,2,3,2,2,2,2,2,2,2,1,1,1, 2,3,2,2,2,2,1,1,1,1,1,1,1,2,1/
data nqp /2*2,8*2,8*3,18*4,18*5,13*6,2,2,2,2,2,2,3,  2,2,2,3,2,2,2,2,2,2,2,2,2,2, 2,3,2,2,2,2,2,2,2,2,2,2,2,2,2/
data nqd /30*3,18*4,31*5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/
character(3),dimension(103)::periodic=(/'H  ','HE ',& 
                'LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ',& 
                 'NA ','MG ','AL ','SI ','P  ','S  ','CL ','AR ',& 
                 'K  ','CA ',& 
                 'SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',& 
                           'GA ','GE ','AS ','SE ','BR ','KR ',& 
                 'RB ','SR ',& 
                 'Y  ','ZR ','NB ','MO ','TC ','RU ','RH ','PD ','AG ','CD ',& 
                           'IN ','SN ','SB ','TE ','I  ','XE ',& 
                 'CS ','BA ',& 
                 'LA ','CE ','PR ','ND ','PM ','SM ','EU ','GD ',& 
                 'TB ','DY ','HO ','H1 ','C2 ','C3 ','C4 ' ,&
                      'C5 ','O6 ','P7 ','C1 ','C2 ','O3','P4 ','O5 ','C6 ','C7 ',& 
                           'O8 ','N9 ','C10','C11','H12','H13',& 
                 'H14','O1 ',& 
                 'P2 ','O3 ' ,'C4 ','O5 ','C6 ','H7 ','H8 ','H9 ',& 
                 'H10','H11','H12','H13','O14','H15'             /)

double precision,dimension(103)::atmass=(/ 1.007825D+00 , 4.00260D+00 , 7.01600D+00 , &
          9.01218D+00  , 11.00931D+00 , 12.00000D+00 , &
          14.00307D+00 , 15.99491D+00 , 18.99840D+00 , &
          19.99244D+00 , 22.98980D+00 , 23.98504D+00 , &
          26.98153D+00 , 27.97693D+00 , 30.97376D+00 , &
          31.97207D+00 , 34.96885D+00 , 39.96238D+00  , &
          38.96371D+00 , 39.96259D+00 , 44.95591D+00  , &
          47.94795D+00 , 50.94396D+00 , 51.94051D+00  , &
          54.93805D+00 , 55.93494D+00 , 58.93320D+00  , &
          57.93535D+00 , 62.93960D+00 , 63.92915D+00  , &
          68.92558D+00 , 73.92118D+00 , 74.92159D+00  , &
          79.91650D+00 , 78.91830D+00 , 83.91150D+00  , &
          84.91180D+00 , 87.90560D+00 , 88.90580D+00  , &
          89.90470D+00 , 92.90640D+00 , 97.90540D+00  , &
          00.00000D+00 , 101.9043D+00 , 102.9055D+00  , &
          105.9032D+00 , 106.9050D+00 , 113.9036D+00  , &
          114.9041D+00 , 117.9034D+00 , 120.9038D+00  , &
          129.9067D+00 , 126.9044D+00 , 131.9042D+00  , &
          132.9054D+00 , 137.9052D+00 , 138.9063D+00  , &
          139.9054D+00 , 140.9076D+00 , 141.9077D+00  , &
          144.9127D+00 , 151.9197D+00 , 152.9212D+00  , &
          1.007825D+00 , 1.007825D+00 , 1.007825D+00  , & 
          157.9241D+00 , 158.9253D+00 , 163.9292D+00  , &
          164.9303D+00 ,   1.007825D+00 , 12.00000D+00   , &
          12.00000D+00  , 12.00000D+00  , 12.00000D+00   , &
          12.00000D+00 , 15.99491D+00 ,   30.97376D+00 , &
          15.99491D+00 , 12.00000D+00 , 12.00000D+00 , &
          15.99491D+00 , 14.00307D+00 , 12.00000D+00  , &
          12.00000D+00 , 1.007825D+00,  1.007825D+00  , &
          1.007825D+00 , 15.99491D+00 , 30.97376D+00  , & ! starts at 88
          15.99491D+00 , 12.00000D+00 , 15.99491D+00  , &
          12.00000D+00 , 1.007825D+00,  1.007825D+00 ,1.007825D+00 ,&
          1.007825D+00 ,1.007825D+00 ,1.007825D+00, 1.007825D+00 ,15.99491D+00,1.007825D+00/)
! sparkle atoms
integer,dimension(4)::spcharge(4)=(/2.0,1.0,-2.0,-1.0/)
double precision,dimension(4)::p0sparkle(4)=(/1.,1.,1.,1./)
double precision ,dimension(:), allocatable ::xsparkle,ysparkle,zsparkle
integer,dimension(:),allocatable::zeffsp
double precision,dimension(4)::alphasp(4)=(/1.5,1.5,1.5,1.5/)


! sam specific stuff
double precision,dimension(:,:),allocatable::rpxyz,rpdxyz,rpd2xyz
integer::ircpt



! parallel NDDO
integer,dimension(:,:),allocatable::neighbors,hpair,tpair,map_pairs
integer,dimension(:),allocatable::iatom_pairs_cpu,ifirst_pair_on_cpu,ilast_pair_on_cpu,ifirstbf,ilastbf,ifirstbf2,ilastbf2,atoms_on_cpu,ifirst_atom_on_cpu,ilast_atom_on_cpu,list_length
integer,dimension(:),allocatable::ifirsthc1c,ilasthc1c,ifirsthc2c,ilasthc2c,num_neighbors,denslist_map2,ifirst_neighbor,denslist_map1
integer::hdim_max,itotalh,itotal,denslist1_cpu,ineighbors_total,denslist2_cpu,icols_cpu,nvecs_global,icols_cpu_max
double precision,dimension(:),allocatable::pair_dens1,pair_dens2,occnum_on_cpu
integer,dimension(:),allocatable::occupied_on_cpu,ifirst_occ_on_cpu,ilast_occ_on_cpu
double precision,dimension(:,:),allocatable::vecstt

!PBCs
double precision,dimension(:,:),allocatable::fracs
double precision,dimension(3,3)::cell


 end module tables

module gaussian_basis
double precision,dimension(:),allocatable::twoeints
integer,dimension(:,:),allocatable::bfuncs,bfuncs_diff
double precision,dimension(1296)::diff_shell
integer::itotal,jmaxints
logical::debug_ab2

double precision,dimension(200)::denom
double precision,dimension(4000)::table
double precision::dkabfactor=5.914967173d0,kabkcd,lindep_tol
integer::num_s_shell,num_p_shell,num_d_shell,ishells_total,natoms,inuc,iveclength,num_gauss,nelectrons
integer,dimension(:),allocatable::iang,icol,iprim,isorb,iporb,idorb,ioffset,ioffset2
integer,dimension(:,:),allocatable::iss_pairs,ips_pairs,ipp_pairs,ids_pairs,idp_pairs,idd_pairs
integer::iss_pairs_total,ips_pairs_total,ipp_pairs_total,ids_pairs_total,idp_pairs_total,idd_pairs_total

double precision,allocatable,dimension(:,:)::shell_exp,shell_coefs,old_shell_exp,old_shell_coefs,shell_coor,old_shell_coor
integer,allocatable,dimension(:)::num_contr,old_num_contr,iangm,old_iangm
integer,allocatable,dimension(:,:)::index_basis
double precision,dimension(3)::psssm_global,psssints_mp1_global,psints_global,dkpsints_global,psintsm_global,psintsmp1_global
double precision,allocatable,dimension(:,:)::coor,density_global,density_q
double precision,allocatable,dimension(:)::zcore,glocal


double precision::ssglobal,zeta_global,dkss_global,ssm_global,ssmp1_global,zeta_plus_eta,dkab_global,dkcd_global!,zabcd 
double precision::ssss_m_global,ssss_mp1_global,repulsion_nuc
double precision,dimension(3,3)::ppssints_global,dsints_global,dkdsints_global,ppints_global,dsintsm_global,dsintsmp1_global&
,dsssm_global,dsssmp1_global,ppssintsmp1_global
double precision,dimension(3,3,3)::dpssm_global,dpssmp1_global,dspsintsmp1_global,dspsintsm_global
double precision,dimension(3,3,3,3)::ddssm_global,ddssmp1_global,dsdsmp1_global
double precision,dimension(3,3,3,3,3)::ddpsintsm_global,ddpsintsmp1_global
double precision,dimension(3,3,3,3,3,3)::dddsintsmp1_global,dddsintsm_global

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


integer,dimension(3,3)::isquish=reshape( source =(/ &
1, 2, 3, 4, 5, 6, 7, 8, 9 /), shape=(/3,3/) )
!integer,dimension(3,3)::isquish=reshape( source =(/ &
!1, 1, 1, 1, 1, 8, 3, 6, 9 /), shape=(/3,3/) )



! parallel stuff
integer,dimension(:),allocatable::iss_pairs_cpu,ips_pairs_cpu,ipp_pairs_cpu,ids_pairs_cpu,idp_pairs_cpu,idd_pairs_cpu,ifirst_on_cpu_ss,ilast_on_cpu_ss,ifirst_on_cpu_ps,ilast_on_cpu_ps &
,ifirst_on_cpu_pp,ilast_on_cpu_pp,iquads_cpu,ifirst_quad_cpu,ilast_quad_cpu,iall_on_cpu,ifirst_on_cpu_ds,ilast_on_cpu_ds,ifirst_on_cpu_dp,ilast_on_cpu_dp,ifirst_on_cpu_dd,ilast_on_cpu_dd &
,ibfcenter
integer::myrank,nprocs


end module gaussian_basis

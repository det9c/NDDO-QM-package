subroutine hfdftuhf
use constants
use tables
use scratch_array
implicit double precision (a-h,o-z)

interface
subroutine mkduhf(rhor,rhora,rhorb,x1,y1,z1)
!double precision,dimension(:,:),intent(in):: density
double precision,intent(inout)::rhor,rhora,rhorb
double precision,intent(in)::x1,y1,z1
end subroutine mkduhf

subroutine hfexcuhf(acoulomb,bcoulomb,abcoulomb,ahfx,bhfx,ekinetic)
double precision,intent(inout)::acoulomb,bcoulomb,abcoulomb,ahfx,bhfx,ekinetic
end subroutine hfexcuhf
end interface

parameter(nquad=48)
dimension weight(nquad),point(nquad),wn(numat),rij(numat,numat) &
,arad(18),wleb(194),pxleb(194),pyleb(194),pzleb(194),thleb(194),pleb(194)
parameter(constlda=-.738558766d0,fourthird=4.0d0/3.0d0)
print*,'-----------------------------------------------'
print*,'|                                             |'
print*,'|            UHF   -   DFT                    |'
print*,'-----------------------------------------------'
print*,''
write(*,*)'Density functionals will be evaluated using unrestricted HF density'
call cpusec(time1)
!s3=zero
!open(unit=1,file='SCFDENS')
!read(1,*)s3
!close(1)
! points and weights for 48 point Gauss Legendre quadrature

      weight(1)=.064737696812683922503d0
      weight(2)=.064466164435950082207d0
      weight(3)=.063924238584648186624d0
      weight(4)=.063114192286254025657d0
      weight(5)=.062039423159892663904d0
      weight(6)=.060704439165893880053d0
      weight(7)=.059114839698395635746d0
      weight(8)=.057277292100403215705d0
      weight(9)=.055199503699984162868d0
      weight(10)=.052890189485193667096d0
      weight(11)=.050359035553854474958d0
      weight(12)=.047616658492490474826d0
      weight(13)=.044674560856694280419d0
      weight(14)=.041545082943464749214d0
      weight(15)=.038241351065830706317d0
      weight(16)=.034777222564770438893d0
      weight(17)=.031167227832798088902d0
      weight(18)=.027426509708356948200d0
      weight(19)=.023570760839324379141d0
      weight(20)=.019616160457355527814d0
      weight(21)=.015579315722943848728d0
      weight(22)=.011477234579234539490d0
      weight(23)=.007327553901276262102d0
      weight(24)=.003153346052305838633d0
      weight(25)=.064737696812683922503d0
      weight(26)=.064466164435950082207d0
      weight(27)=.063924238584648186624d0
      weight(28)=.063114192286254025657d0
      weight(29)=.062039423159892663904d0
      weight(30)=.060704439165893880053d0
      weight(31)=.059114839698395635746d0
      weight(32)=.057277292100403215705d0
      weight(33)=.055199503699984162868d0
      weight(34)=.052890189485193667096d0
      weight(35)=.050359035553854474958d0
      weight(36)=.047616658492490474826d0
      weight(37)=.044674560856694280419d0
      weight(38)=.041545082943464749214d0
      weight(39)=.038241351065830706317d0
      weight(40)=.034777222564770438893d0
      weight(41)=.031167227832798088902d0
      weight(42)=.027426509708356948200d0
      weight(43)=.023570760839324379141d0
      weight(44)=.019616160457355527814d0
      weight(45)=.015579315722943848728d0
      weight(46)=.011477234579234539490d0
      weight(47)=.007327553901276262102d0
      weight(48)=.003153346052305838633d0



      point(1)=.032380170962869362033d0
      point(2)=.097004699209462698930d0
      point(3)=.161222356068891718056d0
      point(4)=.224763790394689061225d0
      point(5)=.287362487355455576736d0
      point(6)=.348755886292160738160d0
      point(7)=.408686481990716729916d0
      point(8)=.466902904750958404545d0
      point(9)=.523160974722233033678d0
      point(10)=.577224726083972703818d0
      point(11)=.628867396776513623995d0
      point(12)=.677872379632663905212d0
      point(13)=.724034130923814654674d0
      point(14)=.767159032515740339254d0
      point(15)=.807066204029442627083d0
      point(16)=.843588261624393530711d0
      point(17)=.876572020274247885906d0
      point(18)=.905879136715569672822d0
      point(19)=.931386690706554333114d0
      point(20)=.952987703160430860723d0
      point(21)=.970591592546247250461d0
      point(22)=.984124583722826857745d0
      point(23)=.993530172266350757548d0
      point(24)=.998771007252426118601d0
      point(25)=-.032380170962869362033d0
      point(26)=-.097004699209462698930d0
      point(27)=-.161222356068891718056d0
      point(28)=-.224763790394689061225d0
      point(29)=-.287362487355455576736d0
      point(30)=-.348755886292160738160d0
      point(31)=-.408686481990716729916d0
      point(32)=-.466902904750958404545d0
      point(33)=-.523160974722233033678d0
      point(34)=-.577224726083972703818d0
      point(35)=-.628867396776513623995d0
      point(36)=-.677872379632663905212d0
      point(37)=-.724034130923814654674d0
      point(38)=-.767159032515740339254d0
      point(39)=-.807066204029442627083d0
      point(40)=-.843588261624393530711d0
      point(41)=-.876572020274247885906d0
      point(42)=-.905879136715569672822d0
      point(43)=-.931386690706554333114d0
      point(44)=-.952987703160430860723d0
      point(45)=-.970591592546247250461d0
      point(46)=-.984124583722826857745d0
      point(47)=-.993530172266350757548d0
      point(48)=-.998771007252426118601d0

! end of pt and wt list for Gauss Legendre quadrature


 thleb( 1 )= 1.5707963267948965
 thleb( 2 )= 1.5707963267948965
 thleb( 3 )= 1.5707963267948965
 thleb( 4 )= 1.5707963267948965
 thleb( 5 )= 0.0E+0
 thleb( 6 )= 3.14159
 thleb( 7 )= 0.7853981633974448
 thleb( 8 )= 0.7853981633974448
 thleb( 9 )= 2.3561944901923483
 thleb( 10 )= 2.3561944901923483
 thleb( 11 )= 0.7853981633974448
 thleb( 12 )= 0.7853981633974448
 thleb( 13 )= 2.3561944901923483
 thleb( 14 )= 2.3561944901923483
 thleb( 15 )= 1.5707963267948965
 thleb( 16 )= 1.5707963267948965
 thleb( 17 )= 1.5707963267948965
 thleb( 18 )= 1.5707963267948965
 thleb( 19 )= 0.9553166181245042
 thleb( 20 )= 0.9553166181245042
 thleb( 21 )= 0.9553166181245042
 thleb( 22 )= 0.9553166181245042
 thleb( 23 )= 2.1862760354652892
 thleb( 24 )= 2.1862760354652892
 thleb( 25 )= 2.1862760354652892
 thleb( 26 )= 2.1862760354652892
 thleb( 27 )= 1.2511856323329721
 thleb( 28 )= 1.2511856323329721
 thleb( 29 )= 1.2511856323329721
 thleb( 30 )= 1.2511856323329721
 thleb( 31 )= 1.890407021256821
 thleb( 32 )= 1.890407021256821
 thleb( 33 )= 1.890407021256821
 thleb( 34 )= 1.890407021256821
 thleb( 35 )= 0.8348385661418226
 thleb( 36 )= 0.8348385661418226
 thleb( 37 )= 0.8348385661418226
 thleb( 38 )= 0.8348385661418226
 thleb( 39 )= 2.3067540874479703
 thleb( 40 )= 2.3067540874479703
 thleb( 41 )= 2.3067540874479703
 thleb( 42 )= 2.3067540874479703
 thleb( 43 )= 0.8348385661418226
 thleb( 44 )= 0.8348385661418226
 thleb( 45 )= 0.8348385661418226
 thleb( 46 )= 0.8348385661418226
 thleb( 47 )= 2.3067540874479703
 thleb( 48 )= 2.3067540874479703
 thleb( 49 )= 2.3067540874479703
 thleb( 50 )= 2.3067540874479703
 thleb( 51 )= 0.42141976340669995
 thleb( 52 )= 0.42141976340669995
 thleb( 53 )= 0.42141976340669995
 thleb( 54 )= 0.42141976340669995
 thleb( 55 )= 2.7201728901830932
 thleb( 56 )= 2.7201728901830932
 thleb( 57 )= 2.7201728901830932
 thleb( 58 )= 2.7201728901830932
 thleb( 59 )= 1.2773566640691322
 thleb( 60 )= 1.2773566640691322
 thleb( 61 )= 1.2773566640691322
 thleb( 62 )= 1.2773566640691322
 thleb( 63 )= 1.8642359895206608
 thleb( 64 )= 1.8642359895206608
 thleb( 65 )= 1.8642359895206608
 thleb( 66 )= 1.8642359895206608
 thleb( 67 )= 1.2773566640691322
 thleb( 68 )= 1.2773566640691322
 thleb( 69 )= 1.2773566640691322
 thleb( 70 )= 1.2773566640691322
 thleb( 71 )= 1.8642359895206608
 thleb( 72 )= 1.8642359895206608
 thleb( 73 )= 1.8642359895206608
 thleb( 74 )= 1.8642359895206608
 thleb( 75 )= 0.6801264219219272
 thleb( 76 )= 0.6801264219219272
 thleb( 77 )= 0.6801264219219272
 thleb( 78 )= 0.6801264219219272
 thleb( 79 )= 2.461466231667866
 thleb( 80 )= 2.461466231667866
 thleb( 81 )= 2.461466231667866
 thleb( 82 )= 2.461466231667866
 thleb( 83 )= 1.109964495414794
 thleb( 84 )= 1.109964495414794
 thleb( 85 )= 1.109964495414794
 thleb( 86 )= 1.109964495414794
 thleb( 87 )= 2.0316281581749993
 thleb( 88 )= 2.0316281581749993
 thleb( 89 )= 2.0316281581749993
 thleb( 90 )= 2.0316281581749993
 thleb( 91 )= 1.109964495414794
 thleb( 92 )= 1.109964495414794
 thleb( 93 )= 1.109964495414794
 thleb( 94 )= 1.109964495414794
 thleb( 95 )= 2.0316281581749993
 thleb( 96 )= 2.0316281581749993
 thleb( 97 )= 2.0316281581749993
 thleb( 98 )= 2.0316281581749993
 thleb( 99 )= 0.18480390510718866
 thleb( 100 )= 0.18480390510718866
 thleb( 101 )= 0.18480390510718866
 thleb( 102 )= 0.18480390510718866
 thleb( 103 )= 2.9567887484826047
 thleb( 104 )= 2.9567887484826047
 thleb( 105 )= 2.9567887484826047
 thleb( 106 )= 2.9567887484826047
 thleb( 107 )= 1.440494370798386
 thleb( 108 )= 1.440494370798386
 thleb( 109 )= 1.440494370798386
 thleb( 110 )= 1.440494370798386
 thleb( 111 )= 1.7010982827914071
 thleb( 112 )= 1.7010982827914071
 thleb( 113 )= 1.7010982827914071
 thleb( 114 )= 1.7010982827914071
 thleb( 115 )= 1.440494370798386
 thleb( 116 )= 1.440494370798386
 thleb( 117 )= 1.440494370798386
 thleb( 118 )= 1.440494370798386
 thleb( 119 )= 1.7010982827914071
 thleb( 120 )= 1.7010982827914071
 thleb( 121 )= 1.7010982827914071
 thleb( 122 )= 1.7010982827914071
 thleb( 123 )= 1.5707963267948965
 thleb( 124 )= 1.5707963267948965
 thleb( 125 )= 1.5707963267948965
 thleb( 126 )= 1.5707963267948965
 thleb( 127 )= 1.5707963267948965
 thleb( 128 )= 1.5707963267948965
 thleb( 129 )= 1.5707963267948965
 thleb( 130 )= 1.5707963267948965
 thleb( 131 )= 0.3530595115393036
 thleb( 132 )= 0.3530595115393036
 thleb( 133 )= 2.7885331420504897
 thleb( 134 )= 2.7885331420504897
 thleb( 135 )= 1.2177368152555956
 thleb( 136 )= 1.2177368152555956
 thleb( 137 )= 1.9238558383341977
 thleb( 138 )= 1.9238558383341977
 thleb( 139 )= 0.3530595115393036
 thleb( 140 )= 0.3530595115393036
 thleb( 141 )= 2.7885331420504897
 thleb( 142 )= 2.7885331420504897
 thleb( 143 )= 1.2177368152555956
 thleb( 144 )= 1.2177368152555956
 thleb( 145 )= 1.9238558383341977
 thleb( 146 )= 1.9238558383341977
 thleb( 147 )= 1.0179418913750843
 thleb( 148 )= 1.0179418913750843
 thleb( 149 )= 1.0179418913750843
 thleb( 150 )= 1.0179418913750843
 thleb( 151 )= 2.1236507622147087
 thleb( 152 )= 2.1236507622147087
 thleb( 153 )= 2.1236507622147087
 thleb( 154 )= 2.1236507622147087
 thleb( 155 )= 0.580778034089786
 thleb( 156 )= 0.580778034089786
 thleb( 157 )= 0.580778034089786
 thleb( 158 )= 0.580778034089786
 thleb( 159 )= 2.5608146195000074
 thleb( 160 )= 2.5608146195000074
 thleb( 161 )= 2.5608146195000074
 thleb( 162 )= 2.5608146195000074
 thleb( 163 )= 1.0179418913750843
 thleb( 164 )= 1.0179418913750843
 thleb( 165 )= 1.0179418913750843
 thleb( 166 )= 1.0179418913750843
 thleb( 167 )= 2.1236507622147087
 thleb( 168 )= 2.1236507622147087
 thleb( 169 )= 2.1236507622147087
 thleb( 170 )= 2.1236507622147087
 thleb( 171 )= 1.4110763938431858
 thleb( 172 )= 1.4110763938431858
 thleb( 173 )= 1.4110763938431858
 thleb( 174 )= 1.4110763938431858
 thleb( 175 )= 1.7305162597466074
 thleb( 176 )= 1.7305162597466074
 thleb( 177 )= 1.7305162597466074
 thleb( 178 )= 1.7305162597466074
 thleb( 179 )= 0.580778034089786
 thleb( 180 )= 0.580778034089786
 thleb( 181 )= 0.580778034089786
 thleb( 182 )= 0.580778034089786
 thleb( 183 )= 2.5608146195000074
 thleb( 184 )= 2.5608146195000074
 thleb( 185 )= 2.5608146195000074
 thleb( 186 )= 2.5608146195000074
 thleb( 187 )= 1.4110763938431858
 thleb( 188 )= 1.4110763938431858
 thleb( 189 )= 1.4110763938431858
 thleb( 190 )= 1.4110763938431858
 thleb( 191 )= 1.7305162597466074
 thleb( 192 )= 1.7305162597466074
 thleb( 193 )= 1.7305162597466074
 thleb( 194 )= 1.7305162597466074


 pleb( 1 )= 0.0E+0
 pleb( 2 )= 3.14159
 pleb( 3 )= 1.570795
 pleb( 4 )= 4.712384999999999
 pleb( 5 )= 0.0E+0
 pleb( 6 )= 0.0E+0
 pleb( 7 )= 1.570795
 pleb( 8 )= 4.712384999999999
 pleb( 9 )= 1.570795
 pleb( 10 )= 4.712384999999999
 pleb( 11 )= 0.0E+0
 pleb( 12 )= 3.141592653589793
 pleb( 13 )= 0.0E+0
 pleb( 14 )= 3.141592653589793
 pleb( 15 )= 0.7853981633974448
 pleb( 16 )= 2.3561944901923483
 pleb( 17 )= 5.497781836602555
 pleb( 18 )= 3.9269855098076513
 pleb( 19 )= 0.7853981633974374
 pleb( 20 )= 2.356194490192356
 pleb( 21 )= 5.497781836602562
 pleb( 22 )= 3.926985509807644
 pleb( 23 )= 0.7853981633974373
 pleb( 24 )= 2.356194490192356
 pleb( 25 )= 5.497781836602562
 pleb( 26 )= 3.926985509807644
 pleb( 27 )= 0.7853981633974507
 pleb( 28 )= 2.3561944901923426
 pleb( 29 )= 5.497781836602549
 pleb( 30 )= 3.926985509807657
 pleb( 31 )= 0.7853981633974507
 pleb( 32 )= 2.3561944901923426
 pleb( 33 )= 5.497781836602549
 pleb( 34 )= 3.926985509807657
 pleb( 35 )= 0.4377579511459213
 pleb( 36 )= 2.703834702443872
 pleb( 37 )= 5.845422048854078
 pleb( 38 )= 3.5793452975561277
 pleb( 39 )= 0.4377579511459221
 pleb( 40 )= 2.703834702443871
 pleb( 41 )= 5.845422048854077
 pleb( 42 )= 3.5793452975561286
 pleb( 43 )= 1.1330383756489856
 pleb( 44 )= 2.0085542779408074
 pleb( 45 )= 5.150141624351014
 pleb( 46 )= 4.274625722059192
 pleb( 47 )= 1.133038375648986
 pleb( 48 )= 2.0085542779408074
 pleb( 49 )= 5.150141624351014
 pleb( 50 )= 4.274625722059192
 pleb( 51 )= 0.785398163397482
 pleb( 52 )= 2.356194490192311
 pleb( 53 )= 5.497781836602518
 pleb( 54 )= 3.9269855098076886
 pleb( 55 )= 0.785398163397482
 pleb( 56 )= 2.356194490192311
 pleb( 57 )= 5.497781836602518
 pleb( 58 )= 3.9269855098076886
 pleb( 59 )= 1.2638358248914025
 pleb( 60 )= 1.8777568286983907
 pleb( 61 )= 5.019344175108597
 pleb( 62 )= 4.405423171301609
 pleb( 63 )= 1.2638358248914025
 pleb( 64 )= 1.8777568286983907
 pleb( 65 )= 5.019344175108597
 pleb( 66 )= 4.405423171301609
 pleb( 67 )= 0.3069605019035154
 pleb( 68 )= 2.834632151686278
 pleb( 69 )= 5.9762194980964844
 pleb( 70 )= 3.448547848313722
 pleb( 71 )= 0.3069605019035154
 pleb( 72 )= 2.834632151686278
 pleb( 73 )= 5.9762194980964844
 pleb( 74 )= 3.448547848313722
 pleb( 75 )= 0.785398163397451
 pleb( 76 )= 2.356194490192342
 pleb( 77 )= 5.497781836602549
 pleb( 78 )= 3.9269855098076575
 pleb( 79 )= 0.7853981633974513
 pleb( 80 )= 2.356194490192342
 pleb( 81 )= 5.497781836602549
 pleb( 82 )= 3.9269855098076575
 pleb( 83 )= 1.0512513401453431
 pleb( 84 )= 2.09034131344445
 pleb( 85 )= 5.231928659854656
 pleb( 86 )= 4.19283868655555
 pleb( 87 )= 1.0512513401453431
 pleb( 88 )= 2.0903413134444503
 pleb( 89 )= 5.231928659854656
 pleb( 90 )= 4.192838686555549
 pleb( 91 )= 0.5195449866495565
 pleb( 92 )= 2.6220476669402366
 pleb( 93 )= 5.763635013350443
 pleb( 94 )= 3.661132333059763
 pleb( 95 )= 0.5195449866495563
 pleb( 96 )= 2.622047666940237
 pleb( 97 )= 5.763635013350443
 pleb( 98 )= 3.6611323330597627
 pleb( 99 )= 0.7853981633975186
 pleb( 100 )= 2.3561944901922746
 pleb( 101 )= 5.497781836602481
 pleb( 102 )= 3.926985509807725
 pleb( 103 )= 0.7853981633975177
 pleb( 104 )= 2.3561944901922755
 pleb( 105 )= 5.497781836602482
 pleb( 106 )= 3.926985509807724
 pleb( 107 )= 1.439373887860492
 pleb( 108 )= 1.7022187657293013
 pleb( 109 )= 4.843806112139507
 pleb( 110 )= 4.580961234270698
 pleb( 111 )= 1.439373887860492
 pleb( 112 )= 1.7022187657293013
 pleb( 113 )= 4.843806112139507
 pleb( 114 )= 4.580961234270698
 pleb( 115 )= 0.13142243893442384
 pleb( 116 )= 3.0101702146553694
 pleb( 117 )= 6.151757561065576
 pleb( 118 )= 3.2730097853446302
 pleb( 119 )= 0.13142243893442384
 pleb( 120 )= 3.0101702146553694
 pleb( 121 )= 6.151757561065576
 pleb( 122 )= 3.2730097853446302
 pleb( 123 )= 1.2177368152555956
 pleb( 124 )= 1.9238558383341977
 pleb( 125 )= 5.065443184744404
 pleb( 126 )= 4.359324161665802
 pleb( 127 )= 0.3530595115393036
 pleb( 128 )= 2.7885331420504897
 pleb( 129 )= 5.930120488460696
 pleb( 130 )= 3.49464685794951
 pleb( 131 )= 6.28317988079071
 pleb( 132 )= 3.141587465619496
 pleb( 133 )= 6.28317988361821
 pleb( 134 )= 3.1415874627919957
 pleb( 135 )= 6.283179955296516
 pleb( 136 )= 3.1415873911136902
 pleb( 137 )= 6.283179955296516
 pleb( 138 )= 3.1415873911136902
 pleb( 139 )= 1.570795
 pleb( 140 )= 4.712384999999999
 pleb( 141 )= 1.570795
 pleb( 142 )= 4.712384999999999
 pleb( 143 )= 1.570795
 pleb( 144 )= 4.712384999999999
 pleb( 145 )= 1.570795
 pleb( 146 )= 4.712384999999999
 pleb( 147 )= 1.382809425929827
 pleb( 148 )= 1.7587832276599661
 pleb( 149 )= 4.900370574070173
 pleb( 150 )= 4.524396772340033
 pleb( 151 )= 1.382809425929827
 pleb( 152 )= 1.7587832276599661
 pleb( 153 )= 4.900370574070173
 pleb( 154 )= 4.524396772340033
 pleb( 155 )= 1.276710248794553
 pleb( 156 )= 1.8648824047952403
 pleb( 157 )= 5.006469751205447
 pleb( 158 )= 4.418297595204759
 pleb( 159 )= 1.2767102487945527
 pleb( 160 )= 1.8648824047952403
 pleb( 161 )= 5.006469751205447
 pleb( 162 )= 4.418297595204759
 pleb( 163 )= 0.18798690086507344
 pleb( 164 )= 2.95360575272472
 pleb( 165 )= 6.095193099134926
 pleb( 166 )= 3.32957424727528
 pleb( 167 )= 0.18798690086507402
 pleb( 168 )= 2.953605752724719
 pleb( 169 )= 6.095193099134926
 pleb( 170 )= 3.3295742472752807
 pleb( 171 )= 0.5608291557018653
 pleb( 172 )= 2.580763497887928
 pleb( 173 )= 5.722350844298134
 pleb( 174 )= 3.702416502112072
 pleb( 175 )= 0.5608291557018653
 pleb( 176 )= 2.580763497887928
 pleb( 177 )= 5.722350844298134
 pleb( 178 )= 3.702416502112072
 pleb( 179 )= 0.29408607800035097
 pleb( 180 )= 2.8475065755894424
 pleb( 181 )= 5.989093921999649
 pleb( 182 )= 3.4356734244105573
 pleb( 183 )= 0.29408607800034985
 pleb( 184 )= 2.8475065755894433
 pleb( 185 )= 5.9890939219996495
 pleb( 186 )= 3.4356734244105564
 pleb( 187 )= 1.0099671710930326
 pleb( 188 )= 2.1316254824967604
 pleb( 189 )= 5.273212828906967
 pleb( 190 )= 4.151554517503239
 pleb( 191 )= 1.0099671710930326
 pleb( 192 )= 2.1316254824967604
 pleb( 193 )= 5.273212828906967
 pleb( 194 )= 4.151554517503239


! enter bragg-slater radii

arad(1)=.25d0
arad(2)=.25d0
arad(3)=1.45d0
arad(4)=1.05d0
arad(5)=.85d0
arad(6)=.70d0
arad(7)=.65d0
arad(8)=.60d0
arad(9)=.50d0
arad(10)=0.0000000
arad(11)=1.8d0
arad(12)=1.5d0
arad(13)=1.25d0
arad(14)=1.10d0
arad(15)=1.00d0
arad(16)=1.00d0
arad(17)=1.00d0
arad(18)=1.00d0
a=a/two/autoang


!
!   The points are mapped to real space via r=a(1+x)/(1-x)
!   where x in the point and a is the Bragg Slater Radius
!
 

      open(unit=1,file='w')
      read(1,*)wleb
      close(1)

 
! remove this redundant computation after debugging     
x=x/autoang
y=y/autoang
z=z/autoang
      do i=1,numat
      do j=1,numat
      r=(x(j)-x(i))**2+(y(j)-y(i))**2+(z(j)-z(i))**2
      rij(j,i)=dsqrt(r)
      end do
      end do


      
fac=180.0d0/pi
sumx=0.0d0
xalp=0.0d0
clyp=0.0d0
b3lyp=0.0d0
vwn=0.0d0
pbe=0.0d0
b3=0.0d0
b88=0.0d0
maxphi=15
phimax=float(maxphi)
sumy=0.0d0

delta=1.0d-6

!     loop over atoms and integrate over each center
do iatom=1,numat
  do irad=1,nquad
    a=arad(zeff(iatom))
    r=a*(one+point(irad))/(one-point(irad))
  do jleb=1,194
     theta=thleb(jleb)
     phi=pleb(jleb)
     xleb=r*dsin(theta)*dcos(phi)
     yleb=r*dsin(theta)*dsin(phi)
     zleb=r*dcos(theta)
     wn=zero
         xleb=xleb+x(iatom)
         yleb=yleb+y(iatom)
         zleb=zleb+z(iatom)
        
          call fuzzy2(x,y,z,r,xleb,yleb,zleb,wn,numat,rij)
         rhor=0.0d0
! now do one center integrations (i guess)
         call mkduhf(rhor,rhora,rhorb,xleb,yleb,zleb)
       if(rhor.lt.1D-10)goto 33
!       call mkd(density,nbasis,nl,nq,xlda,e,r,x,y,z,coorb)
       sumx=sumx+(weight(irad)*r*r*two*a*rhor/(one-point(irad))**2) &
     *wleb(jleb)*4.0d0*pi*wn(iatom)

       xlda=rhora**(fourthird)+rhorb**(fourthird)
       xalp=xalp+(weight(irad)*r*r*two*a*xlda/(one-point(irad))**2) &
            *wleb(jleb)*4.0d0*pi*wn(iatom)
! density gradient at point
       call mkduhf(rhorp,rhorpa,rhorpb,xleb+delta,yleb,zleb)
       call mkduhf(rhorm,rhorma,rhormb,xleb-delta,yleb,zleb)
       dxa=rhorpa-rhorma
       dxb=rhorpb-rhormb
       call mkduhf(rhorp,rhorpa,rhorpb,xleb,yleb+delta,zleb)
       call mkduhf(rhorm,rhorma,rhormb,xleb,yleb-delta,zleb)
       dya=rhorpa-rhorma
       dyb=rhorpb-rhormb
       call mkduhf(rhorp,rhorpa,rhorpb,xleb,yleb,zleb+delta)
       call mkduhf(rhorm,rhorma,rhormb,xleb,yleb,zleb-delta)
       dza=rhorpa-rhorma
       dzb=rhorpb-rhormb

       dxa=dxa/2.0d0/delta
       dya=dya/2.0d0/delta
       dza=dza/2.0d0/delta
       dxb=dxb/2.0d0/delta
       dyb=dyb/2.0d0/delta
       dzb=dzb/2.0d0/delta

       gaa=dxa**2+dya**2+dza**2
       gbb=dxb**2+dyb**2+dzb**2
       gab=dxa*dxb+dya*dyb+dza*dzb

         call culyp(rhora,rhorb,gaa,gbb,gab,fout,dfdra,dfdrb, &
                          dfdgaa,dfdgbb,dfdgab)

         clyp=clyp+(weight(irad)*r*r*two*a*fout/(one-point(irad))**2) &
             *wleb(jleb)*4.0d0*pi*wn(iatom)
        call xcub3lyp(rhora,rhorb,gaa,gbb,gab,fout,dfdra,dfdrb, &
                             dfdgaa,dfdgbb,dfdgab)
        b3lyp=b3lyp+(weight(irad)*r*r*two*a*fout/(one-point(irad))**2) &
             *wleb(jleb)*4.0d0*pi*wn(iatom)
       
        call c_uks_vwn5(rhora,rhorb,fout,dfdra,dfdrb)

         vwn=vwn+(weight(irad)*r*r*two*a*fout/(one-point(irad))**2) &
              *wleb(jleb)*4.0d0*pi*wn(iatom)
       
         call x_uks_pbe(rhora,rhorb,gaa,gbb,fout,dfdra,dfdrb,dfdgaa,dfdgbb)
       pbe=pbe+(weight(irad)*r*r*two*a*fout/(one-point(irad))**2) &
             *wleb(jleb)*4.0d0*pi*wn(iatom)


        call x_uks_b88(rhora,rhorb,gaa,gbb,fout,dfdra,dfdrb,dfdgaa,dfdgbb)
        b88=b88+(weight(irad)*r*r*two*a*fout/(one-point(irad))**2) &
            *wleb(jleb)*4.0d0*pi*wn(iatom)
   




       
      



33  end do ! end loops over r, theta, phi quadrature points
end do 
end do

call cpusec(time2)
write(*,*)'DFT integrations required ',time2-time1,' seconds'

! get exact exchange and coulomb energies for alpha and beta spin manifolds
call hfexcuhf(acoulomb,bcoulomb,abcoulomb,ahfx,bhfx,ekinetic)
print*,ahfx,bhfx
       PRINT*,'*******************************************'
       print*,'*Density integrates to ',sumx,' electrons *'
       PRINT*,'*******************************************'

       PRINT*,'SUMMARY OF EXCHANGE AND/OR CORRELATION ENERGIES:'
       PRINT*,'X = EXCHANGE ONLY    C = CORRELATION ONLY  XC = COMBINATION'
       print*,'-----------------------------------'
       print*,'SCF Total Energy =',ekinetic+acoulomb+bcoulomb &
+abcoulomb+ahfx+bhfx+enuc
       print*,'SCF One Particle Energy =',ekinetic
       print*,'SCF Coulomb  E(J) = ',acoulomb+bcoulomb+abcoulomb
       print*,'SCF Exchange E(X) = ',(ahfx+bhfx)/autoev     
       print*,'LDA (Xalpha) E(X) = ',xalp*constlda*2.0d0**(1.0d0/3.0d0)
       print*,'PBE          E(X) = ',pbe
       print*,'BECKE88      E(X) = ',b88
       print*,'LYP          E(C) = ',clyp
       print*,'VWN          E(C) = ',vwn
       print*,'B3LYP        E(XC)= ',b3lyp+0.2d0*(ahfx+bhfx)/autoev
print*,'b88+lyp',ekinetic+acoulomb+bcoulomb+abcoulomb &
+b88*autoev+enuc+clyp*autoev

print*,'b3lyp',ekinetic+acoulomb+bcoulomb+abcoulomb+b3lyp*autoev+0.2d0*(ahfx+bhfx)
print*,'lyp',ekinetic+acoulomb+bcoulomb &
+abcoulomb+ahfx+bhfx+enuc+clyp*autoev
      
      end subroutine hfdftuhf

 
[MFDn]
Version = 15
ndiags   =        3
MPIranks  =        6
OMPthreads =      256
 
[PARAMETERS]
 
[Basis]
Nprotons  =        4
Nneutrons =        5
Nshell  =       13
TwoMj  =        1
parity =       -1
Nmin   =        0
Nmax   =        8
DeltaN =        2
WTmax  =    13.1000
 
[Many-body matrix]
dimension  =           63003395
numnonzero =        45472165505
 
[Interaction]
# using HO TBME-J Vpp, Vnn, Vpn in H2 format
Hrank  =        2
hbomeg =    20.0000
fmass  =   938.9200
LamHcm =     3.0000
 
Trelfile = TBME_TrelHO.bin
Hrelfile = TBME_HrelHO.bin
Vnnfile  = VNN_Daejeon16_hw20.0.bin
 
[Observables]
numTBops   =        2
# TBME file for relative R2 operator
TBMEfile(1) = TBME_R2.bin
 
# TBME files for additional operators
TBMEfile(2) = TBME_HCM.bin
 
 
[RESULTS]
# following Condon-Shortley phase conventions
# for conventions of reduced matrix elements and
# for Wigner-Eckart theorem, see Suhonen (2.27)
 
[Energies]
# Seq       J    n      T        Eabs        Eexc        Error    J-full
    1      1.5   1   0.500      -57.109      0.000     0.77E-04    1.5000
    2      2.5   1   0.500      -54.193      2.916     0.78E-04    2.5000
    3      0.5   1   0.500      -51.486      5.623     0.67E-04    0.5000
    4      3.5   1   0.500      -49.488      7.622     0.50E-04    3.5000
    5      1.5   2   0.500      -49.271      7.838     0.89E-04    1.5000
    6      2.5   2   0.500      -46.834     10.276     0.12E-03    2.5000
    7      3.5   2   0.500      -45.213     11.896     0.67E-03    3.5000
    8      0.5   2   0.500      -44.234     12.875     0.18E-03    0.5000
 
[Oscillator quanta]
# Seq    J    n      T      Amp(N)^2 
    1   1.5   1   0.500      0.5739      0.2582      0.1175      0.4094E-01  0.9514E-02
    2   2.5   1   0.500      0.5525      0.2686      0.1243      0.4426E-01  0.1037E-01
    3   0.5   1   0.500      0.4187      0.3250      0.1720      0.6798E-01  0.1633E-01
    4   3.5   1   0.500      0.5680      0.2580      0.1209      0.4301E-01  0.1014E-01
    5   1.5   2   0.500      0.3826      0.3413      0.1838      0.7426E-01  0.1801E-01
    6   2.5   2   0.500      0.4554      0.3067      0.1604      0.6247E-01  0.1507E-01
    7   3.5   2   0.500      0.5832      0.2500      0.1157      0.4141E-01  0.9682E-02
    8   0.5   2   0.500      0.5746      0.2459      0.1249      0.4416E-01  0.1038E-01
 
[M1 moments]
# for M1 moment conventions, see Suhonen (6.52)
# Seq    J    n      T        mu       Dl(p)      Dl(n)      Ds(p)      Ds(n)    sum(=J)
    1   1.5   1   0.500    -1.0195     0.3026     0.8268     0.0102     0.3604     1.5000
    2   2.5   1   0.500    -0.2347     0.9977     1.0969     0.0338     0.3715     2.5000
    3   0.5   1   0.500     0.8542     0.1793     0.4691     0.0114    -0.1597     0.5000
    4   3.5   1   0.500     1.2165     1.6548     1.3620     0.1498     0.3334     3.5000
    5   1.5   2   0.500     1.3480     0.7703     0.8253     0.0225    -0.1181     1.5000
    6   2.5   2   0.500     2.4444     1.4379     1.1869     0.0562    -0.1810     2.5000
    7   3.5   2   0.500     3.7548     1.0634     1.2437     0.7709     0.4220     3.5000
    8   0.5   2   0.500     2.9477     0.1207    -0.2711     0.5648     0.0857     0.5000
 
[E2 moments]
# for Q moment conventions, see Suhonen (6.53)
# Seq    J    n      T       Q(p)       Q(n)
    1   1.5   1   0.500     4.1104     3.5812
    2   2.5   1   0.500    -1.9722     0.0211
    4   3.5   1   0.500    -4.1072    -3.1023
    5   1.5   2   0.500    -4.3336    -3.6260
    6   2.5   2   0.500    -5.6994    -6.3171
    7   3.5   2   0.500     2.9903     5.0729
 
[Angular momenta]
# Seq    J    n      T      <L^2>      <S^2>      <Lp^2>     <Sp^2>     <Ln^2>     <Sn^2>        <J^2>
    1   1.5   1   0.500     3.0327     1.1359     1.2192     0.1831     0.5380     0.8810       3.7500
    2   2.5   1   0.500     7.0935     1.1808     2.4761     0.2151     0.9875     0.8884       8.7500
    3   0.5   1   0.500     2.3053     1.1103     1.1339     0.1790     0.3671     0.8563       0.7500
    4   3.5   1   0.500    12.9229     1.5217     3.8524     0.4071     0.8619     0.9462      15.7500
    5   1.5   2   0.500     5.3484     1.1202     1.7661     0.1935     0.9951     0.8754       3.7500
    6   2.5   2   0.500    10.7760     1.1524     3.5628     0.2311     1.1006     0.8556       8.7500
    7   3.5   2   0.500     8.4494     3.4353     1.4526     1.6709     1.1301     0.9520      15.7500
    8   0.5   2   0.500     1.9649     3.1663     1.1172     1.8543     0.8690     0.8941       0.7500
 
[Relative radii]
# Seq    J    n      T      r(p)       r(n)    r(matter)
    1   1.5   1   0.500     2.2131     2.3066     2.2655
    2   2.5   1   0.500     2.2206     2.3180     2.2752
    3   0.5   1   0.500     2.2635     2.4478     2.3676
    4   3.5   1   0.500     2.2117     2.3003     2.2613
    5   1.5   2   0.500     2.2736     2.4617     2.3800
    6   2.5   2   0.500     2.2534     2.4145     2.3443
    7   3.5   2   0.500     2.2392     2.2904     2.2678
    8   0.5   2   0.500     2.2661     2.3147     2.2932
 
[Other 2-body observables]
# Seq    J    n      T     see header for TBME file names of observables
    1   1.5   1   0.500     0.117890E-04
    2   2.5   1   0.500     0.120529E-04
    3   0.5   1   0.500     0.122331E-04
    4   3.5   1   0.500     0.128017E-04
    5   1.5   2   0.500     0.120699E-04
    6   2.5   2   0.500     0.129290E-04
    7   3.5   2   0.500     0.145594E-04
    8   0.5   2   0.500     0.142744E-04
 
[Occupation probabilities]
# orb      n_rad   l   2*j
    1   pro    0    0    1     1.850        1.846        1.816        1.858        1.808        1.834        1.873        1.885    
    2   pro    0    1    1    0.4004       0.3916       0.4155       0.3279       0.3865       0.4476       0.3495       0.3370    
    3   pro    0    1    3     1.496        1.494        1.461        1.549        1.477        1.417        1.504        1.485    
    4   pro    1    0    1    0.3812E-01   0.3987E-01   0.5445E-01   0.3661E-01   0.5646E-01   0.5080E-01   0.3972E-01   0.4079E-01
    5   pro    0    2    3    0.5676E-01   0.5768E-01   0.6568E-01   0.5310E-01   0.6864E-01   0.5976E-01   0.4319E-01   0.4962E-01
    6   pro    0    2    5    0.5462E-01   0.5589E-01   0.6251E-01   0.5301E-01   0.6645E-01   0.5866E-01   0.5234E-01   0.5485E-01
    7   pro    1    1    1    0.2153E-01   0.2324E-01   0.2747E-01   0.2141E-01   0.2800E-01   0.3135E-01   0.3135E-01   0.2560E-01
    8   pro    1    1    3    0.4078E-01   0.4152E-01   0.5307E-01   0.4089E-01   0.5633E-01   0.4785E-01   0.5177E-01   0.7121E-01
    9   pro    0    3    5    0.1274E-01   0.1558E-01   0.1491E-01   0.1915E-01   0.1784E-01   0.1697E-01   0.1344E-01   0.1032E-01
   10   pro    0    3    7    0.1842E-01   0.2222E-01   0.1713E-01   0.2488E-01   0.2042E-01   0.2030E-01   0.1746E-01   0.8054E-02
   11   pro    2    0    1    0.1308E-02   0.1545E-02   0.1532E-02   0.1994E-02   0.1711E-02   0.1960E-02   0.3403E-02   0.4775E-02
   12   pro    1    2    3    0.9795E-03   0.1098E-02   0.1325E-02   0.1299E-02   0.1510E-02   0.1591E-02   0.2113E-02   0.3309E-02
   13   pro    1    2    5    0.1894E-02   0.2170E-02   0.2166E-02   0.2709E-02   0.2436E-02   0.2641E-02   0.3665E-02   0.4914E-02
   14   pro    0    4    7    0.1917E-02   0.2157E-02   0.2083E-02   0.2124E-02   0.2406E-02   0.2140E-02   0.1863E-02   0.2298E-02
   15   pro    0    4    9    0.6790E-03   0.8037E-03   0.1238E-02   0.8138E-03   0.1311E-02   0.1215E-02   0.7016E-03   0.4436E-03
   16   pro    2    1    1    0.3166E-03   0.3781E-03   0.4027E-03   0.4815E-03   0.4033E-03   0.7457E-03   0.1962E-02   0.2269E-02
   17   pro    2    1    3    0.1952E-02   0.2472E-02   0.1921E-02   0.3876E-02   0.2157E-02   0.3136E-02   0.6799E-02   0.1021E-01
   18   pro    1    3    5    0.2735E-03   0.2902E-03   0.3035E-03   0.3924E-03   0.3921E-03   0.3843E-03   0.5021E-03   0.4993E-03
   19   pro    1    3    7    0.3069E-03   0.3523E-03   0.3567E-03   0.4660E-03   0.4039E-03   0.4400E-03   0.4810E-03   0.4053E-03
   20   pro    0    5    9    0.3149E-03   0.4159E-03   0.2696E-03   0.5342E-03   0.3785E-03   0.3497E-03   0.5922E-03   0.2852E-03
   21   pro    0    5   11    0.5809E-04   0.9408E-04   0.5368E-04   0.1152E-03   0.7891E-04   0.8655E-04   0.8595E-04   0.3506E-04
   22   pro    3    0    1    0.1333E-03   0.1599E-03   0.1195E-03   0.2302E-03   0.1326E-03   0.1922E-03   0.4880E-03   0.6491E-03
   23   pro    2    2    3    0.7263E-04   0.8508E-04   0.8088E-04   0.9441E-04   0.8258E-04   0.1212E-03   0.2252E-03   0.2512E-03
   24   pro    2    2    5    0.1395E-03   0.1626E-03   0.1224E-03   0.2481E-03   0.1365E-03   0.1822E-03   0.4731E-03   0.6881E-03
   25   pro    1    4    7    0.5874E-04   0.6105E-04   0.4620E-04   0.6160E-04   0.5220E-04   0.5738E-04   0.7521E-04   0.9694E-04
   26   pro    1    4    9    0.4126E-04   0.4385E-04   0.3801E-04   0.5266E-04   0.3992E-04   0.4518E-04   0.5728E-04   0.4112E-04
   27   pro    0    6   11    0.6081E-04   0.6426E-04   0.5492E-04   0.8273E-04   0.5934E-04   0.5973E-04   0.8558E-04   0.4732E-04
   28   pro    0    6   13    0.2653E-05   0.3855E-05   0.3402E-05   0.5110E-05   0.4142E-05   0.5651E-05   0.3979E-05   0.1459E-05
   29   pro    3    1    1    0.4972E-04   0.7050E-04   0.4374E-04   0.8005E-04   0.4591E-04   0.1141E-03   0.2900E-03   0.2340E-03
   30   pro    3    1    3    0.2246E-03   0.2789E-03   0.1857E-03   0.4658E-03   0.2120E-03   0.3114E-03   0.9270E-03   0.1338E-02
   31   pro    2    3    5    0.1021E-04   0.9289E-05   0.1072E-04   0.9469E-05   0.1058E-04   0.1212E-04   0.1208E-04   0.1840E-04
   32   pro    2    3    7    0.1131E-04   0.1065E-04   0.8237E-05   0.1323E-04   0.9383E-05   0.9040E-05   0.1973E-04   0.2157E-04
   33   pro    1    5    9    0.4684E-05   0.5461E-05   0.3488E-05   0.6447E-05   0.4291E-05   0.5112E-05   0.7319E-05   0.7770E-05
   34   pro    1    5   11    0.3742E-05   0.4143E-05   0.2599E-05   0.5701E-05   0.2866E-05   0.3435E-05   0.4455E-05   0.2297E-05
   35   pro    0    7   13    0.2508E-05   0.2953E-05   0.2046E-05   0.4722E-05   0.2531E-05   0.2256E-05   0.3936E-05   0.1294E-05
   36   pro    0    7   15    0.3531E-07   0.9653E-07   0.1352E-07   0.1800E-06   0.5079E-07   0.6727E-07   0.9812E-07   0.1088E-07
   37   pro    4    0    1    0.1450E-04   0.1657E-04   0.1163E-04   0.2284E-04   0.1189E-04   0.1762E-04   0.4706E-04   0.6302E-04
   38   pro    3    2    3    0.6124E-05   0.7623E-05   0.5338E-05   0.9722E-05   0.5741E-05   0.1003E-04   0.2357E-04   0.2225E-04
   39   pro    3    2    5    0.1593E-04   0.1808E-04   0.1150E-04   0.2733E-04   0.1198E-04   0.1680E-04   0.4744E-04   0.7091E-04
   40   pro    2    4    7    0.1136E-05   0.1276E-05   0.8813E-06   0.1764E-05   0.1011E-05   0.1116E-05   0.1009E-05   0.9123E-06
   41   pro    2    4    9    0.1631E-05   0.1876E-05   0.9917E-06   0.2616E-05   0.1171E-05   0.1431E-05   0.1751E-05   0.6928E-06
   42   pro    1    6   11    0.3789E-06   0.4445E-06   0.2343E-06   0.7639E-06   0.3371E-06   0.3030E-06   0.8380E-06   0.2037E-06
   43   pro    1    6   13    0.4315E-07   0.1104E-06   0.1196E-07   0.1685E-06   0.5096E-07   0.6649E-07   0.1219E-06   0.5867E-08
   44   pro    0    8   15    0.6520E-07   0.2272E-06   0.2972E-18   0.3020E-06   0.1583E-06   0.1317E-06   0.2213E-06   0.3706E-15
   45   pro    0    8   17    0.1472E-20   0.1563E-12   0.4353E-20   0.1296E-08   0.2065E-18   0.1781E-08   0.6293E-09   0.2653E-17
   46   pro    4    1    1    0.5738E-05   0.7428E-05   0.4934E-05   0.8125E-05   0.5198E-05   0.1061E-04   0.2646E-04   0.1941E-04
   47   pro    4    1    3    0.2340E-04   0.2614E-04   0.1763E-04   0.3897E-04   0.1781E-04   0.2476E-04   0.7330E-04   0.1142E-03
   48   pro    3    3    5    0.2499E-06   0.4763E-06   0.2878E-06   0.1395E-05   0.5547E-06   0.5425E-06   0.5784E-06   0.1004E-06
   49   pro    3    3    7    0.1237E-05   0.1844E-05   0.6172E-06   0.3263E-05   0.9329E-06   0.1597E-05   0.2313E-05   0.5985E-06
   50   pro    2    5    9    0.1416E-07   0.5076E-07   0.1924E-09   0.1209E-06   0.3651E-07   0.4082E-07   0.5224E-07   0.3803E-08
   51   pro    2    5   11    0.1047E-07   0.7479E-07   0.1535E-18   0.6717E-07   0.3863E-07   0.3149E-07   0.4066E-07   0.8761E-16
   52   pro    1    7   13    0.9894E-20   0.3085E-10   0.1576E-19   0.1330E-08   0.2607E-18   0.6029E-12   0.7049E-08   0.1342E-17
   53   pro    1    7   15    0.5407E-21   0.4315E-20   0.1294E-20   0.3527E-10   0.1017E-19   0.1062E-19   0.7009E-12   0.1287E-18
   54   pro    0    9   17    0.3761E-21   0.7877E-21   0.3986E-21   0.7299E-21   0.7056E-21   0.2199E-20   0.5168E-18   0.4003E-19
   55   pro    0    9   19    0.8750E-23   0.4805E-22   0.6556E-22   0.5451E-22   0.1279E-21   0.3939E-21   0.1823E-19   0.1013E-20
   56   pro    5    0    1     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   57   pro    4    2    3     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   58   pro    4    2    5     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   59   pro    3    4    7     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   60   pro    3    4    9     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   61   pro    2    6   11     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   62   pro    2    6   13     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   63   pro    1    8   15     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   64   pro    1    8   17     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   65   pro    0   10   19     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   66   pro    0   10   21     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   67   pro    5    1    1     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   68   pro    5    1    3     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   69   pro    4    3    5     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   70   pro    4    3    7     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   71   pro    3    5    9     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   72   pro    3    5   11     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   73   pro    2    7   13     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   74   pro    2    7   15     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   75   pro    1    9   17     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   76   pro    1    9   19     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   77   pro    0   11   21     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   78   pro    0   11   23     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   79   pro    6    0    1     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   80   pro    5    2    3     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   81   pro    5    2    5     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   82   pro    4    4    7     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   83   pro    4    4    9     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   84   pro    3    6   11     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   85   pro    3    6   13     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   86   pro    2    8   15     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   87   pro    2    8   17     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   88   pro    1   10   19     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   89   pro    1   10   21     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   90   pro    0   12   23     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   91   pro    0   12   25     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
# sum of proton occ.prob.     4.000        4.000        4.000        4.000        4.000        4.000        4.000        4.000    
 
    1   neu    0    0    1     1.846        1.842        1.805        1.857        1.794        1.829        1.882        1.899    
    2   neu    0    1    1    0.4012       0.4677       0.9354       0.3850        1.011       0.7510       0.3545       0.2926    
    3   neu    0    1    3     2.373        2.291        1.665        2.381        1.555        1.882        2.425        2.455    
    4   neu    1    0    1    0.4614E-01   0.4921E-01   0.8662E-01   0.4440E-01   0.9316E-01   0.7529E-01   0.4472E-01   0.4436E-01
    5   neu    0    2    3    0.6596E-01   0.6632E-01   0.8789E-01   0.6052E-01   0.9009E-01   0.8090E-01   0.4874E-01   0.5479E-01
    6   neu    0    2    5    0.7260E-01   0.7581E-01   0.8707E-01   0.7085E-01   0.9783E-01   0.7538E-01   0.5763E-01   0.5964E-01
    7   neu    1    1    1    0.2052E-01   0.2927E-01   0.1111       0.2387E-01   0.1234       0.8948E-01   0.2764E-01   0.2035E-01
    8   neu    1    1    3    0.8318E-01   0.8241E-01   0.1123       0.7896E-01   0.1126       0.1087       0.7831E-01   0.1018    
    9   neu    0    3    5    0.2696E-01   0.2565E-01   0.2082E-01   0.2666E-01   0.2331E-01   0.2527E-01   0.2108E-01   0.1765E-01
   10   neu    0    3    7    0.3401E-01   0.3793E-01   0.2604E-01   0.3767E-01   0.3274E-01   0.2372E-01   0.2586E-01   0.1113E-01
   11   neu    2    0    1    0.4205E-02   0.4481E-02   0.8540E-02   0.4609E-02   0.9061E-02   0.7964E-02   0.4856E-02   0.6132E-02
   12   neu    1    2    3    0.1735E-02   0.2146E-02   0.8373E-02   0.2118E-02   0.8245E-02   0.7871E-02   0.2631E-02   0.4034E-02
   13   neu    1    2    5    0.4888E-02   0.5002E-02   0.4036E-02   0.5374E-02   0.4910E-02   0.4165E-02   0.5119E-02   0.6068E-02
   14   neu    0    4    7    0.2139E-02   0.2419E-02   0.2204E-02   0.2567E-02   0.2423E-02   0.2633E-02   0.2847E-02   0.3300E-02
   15   neu    0    4    9    0.7582E-03   0.9766E-03   0.1509E-02   0.1011E-02   0.1622E-02   0.1448E-02   0.9318E-03   0.7708E-03
   16   neu    2    1    1    0.6200E-03   0.1600E-02   0.1992E-01   0.1022E-02   0.1927E-01   0.1725E-01   0.1775E-02   0.2613E-02
   17   neu    2    1    3    0.1088E-01   0.1070E-01   0.6932E-02   0.1119E-01   0.9046E-02   0.6885E-02   0.1029E-01   0.1341E-01
   18   neu    1    3    5    0.3598E-03   0.4596E-03   0.1066E-02   0.5721E-03   0.1019E-02   0.1353E-02   0.6273E-03   0.1649E-02
   19   neu    1    3    7    0.6107E-03   0.7169E-03   0.7256E-03   0.9204E-03   0.1062E-02   0.7056E-03   0.7749E-03   0.5148E-03
   20   neu    0    5    9    0.3667E-03   0.5376E-03   0.2790E-03   0.7774E-03   0.3675E-03   0.5358E-03   0.1057E-02   0.6060E-03
   21   neu    0    5   11    0.9807E-04   0.1629E-03   0.9363E-04   0.1980E-03   0.1636E-03   0.1118E-03   0.1248E-03   0.4408E-04
   22   neu    3    0    1    0.6706E-03   0.7166E-03   0.1679E-02   0.6890E-03   0.1751E-02   0.1493E-02   0.6570E-03   0.8256E-03
   23   neu    2    2    3    0.1564E-03   0.2414E-03   0.1774E-02   0.1824E-03   0.1633E-02   0.1590E-02   0.2404E-03   0.3440E-03
   24   neu    2    2    5    0.7654E-03   0.7494E-03   0.5191E-03   0.7921E-03   0.7672E-03   0.4714E-03   0.6967E-03   0.8689E-03
   25   neu    1    4    7    0.5402E-04   0.5991E-04   0.5948E-04   0.6296E-04   0.6302E-04   0.8287E-04   0.7435E-04   0.1182E-03
   26   neu    1    4    9    0.3911E-04   0.4516E-04   0.3299E-04   0.6466E-04   0.4151E-04   0.4185E-04   0.7264E-04   0.5130E-04
   27   neu    0    6   11    0.5915E-04   0.6206E-04   0.4883E-04   0.7962E-04   0.5330E-04   0.5466E-04   0.8312E-04   0.4822E-04
   28   neu    0    6   13    0.2508E-05   0.4165E-05   0.7174E-05   0.5133E-05   0.6124E-05   0.9993E-05   0.3849E-05   0.1944E-05
   29   neu    3    1    1    0.8619E-04   0.2601E-03   0.3361E-02   0.1336E-03   0.3030E-02   0.2961E-02   0.2570E-03   0.2880E-03
   30   neu    3    1    3    0.1634E-02   0.1585E-02   0.1119E-02   0.1631E-02   0.1647E-02   0.9942E-03   0.1391E-02   0.1748E-02
   31   neu    2    3    5    0.1404E-04   0.1837E-04   0.7622E-04   0.2120E-04   0.6783E-04   0.1127E-03   0.2919E-04   0.1145E-03
   32   neu    2    3    7    0.2926E-04   0.3648E-04   0.1969E-04   0.6756E-04   0.4279E-04   0.1761E-04   0.5638E-04   0.2613E-04
   33   neu    1    5    9    0.4755E-05   0.5810E-05   0.5759E-05   0.6970E-05   0.5833E-05   0.8344E-05   0.8278E-05   0.9879E-05
   34   neu    1    5   11    0.4081E-05   0.4447E-05   0.3008E-05   0.6188E-05   0.3604E-05   0.3818E-05   0.5285E-05   0.2808E-05
   35   neu    0    7   13    0.2351E-05   0.2987E-05   0.1947E-05   0.4868E-05   0.2405E-05   0.2396E-05   0.4262E-05   0.1723E-05
   36   neu    0    7   15    0.3394E-07   0.1152E-06   0.7600E-07   0.2615E-06   0.1016E-06   0.1792E-06   0.1611E-06   0.1651E-07
   37   neu    4    0    1    0.6464E-04   0.6904E-04   0.1827E-03   0.6656E-04   0.1931E-03   0.1612E-03   0.6127E-04   0.7980E-04
   38   neu    3    2    3    0.1549E-04   0.2511E-04   0.2073E-03   0.1899E-04   0.1892E-03   0.1883E-03   0.2440E-04   0.3149E-04
   39   neu    3    2    5    0.8283E-04   0.8082E-04   0.6048E-04   0.8746E-04   0.9572E-04   0.5182E-04   0.7101E-04   0.9116E-04
   40   neu    2    4    7    0.1789E-05   0.1830E-05   0.1707E-05   0.2584E-05   0.1761E-05   0.3775E-05   0.2276E-05   0.4169E-05
   41   neu    2    4    9    0.2595E-05   0.3292E-05   0.1418E-05   0.5141E-05   0.2409E-05   0.1658E-05   0.3335E-05   0.8624E-06
   42   neu    1    6   11    0.3660E-06   0.4218E-06   0.3073E-06   0.7356E-06   0.3813E-06   0.3740E-06   0.9407E-06   0.3056E-06
   43   neu    1    6   13    0.2803E-07   0.9659E-07   0.1818E-07   0.1311E-06   0.4869E-07   0.6219E-07   0.1018E-06   0.3895E-08
   44   neu    0    8   15    0.5746E-07   0.2263E-06   0.3077E-18   0.2762E-06   0.1388E-06   0.9418E-07   0.1749E-06   0.3838E-15
   45   neu    0    8   17    0.2216E-20   0.2861E-11   0.7335E-20   0.1287E-08   0.2007E-18   0.4923E-08   0.5197E-09   0.3979E-17
   46   neu    4    1    1    0.7439E-05   0.2227E-04   0.2984E-03   0.1110E-04   0.2655E-03   0.2669E-03   0.2269E-04   0.2339E-04
   47   neu    4    1    3    0.1313E-03   0.1260E-03   0.9757E-04   0.1313E-03   0.1532E-03   0.8162E-04   0.1079E-03   0.1470E-03
   48   neu    3    3    5    0.1436E-05   0.1558E-05   0.2010E-05   0.2910E-05   0.2098E-05   0.6401E-05   0.2113E-05   0.6867E-05
   49   neu    3    3    7    0.3931E-05   0.5474E-05   0.1399E-05   0.9462E-05   0.3559E-05   0.2001E-05   0.5942E-05   0.7734E-06
   50   neu    2    5    9    0.8805E-07   0.6136E-07   0.2263E-07   0.1487E-06   0.5109E-07   0.6292E-07   0.2122E-06   0.8747E-07
   51   neu    2    5   11    0.2083E-07   0.9373E-07   0.2155E-18   0.6143E-07   0.5500E-07   0.3509E-07   0.2801E-07   0.7116E-16
   52   neu    1    7   13    0.1816E-19   0.2082E-08   0.3193E-19   0.3660E-08   0.4877E-18   0.3793E-09   0.8876E-08   0.1115E-16
   53   neu    1    7   15    0.1253E-20   0.4915E-20   0.2213E-20   0.1326E-09   0.2429E-19   0.1822E-19   0.5271E-10   0.2462E-18
   54   neu    0    9   17    0.4808E-21   0.1162E-20   0.6282E-21   0.1207E-20   0.1040E-20   0.8250E-20   0.2710E-18   0.8234E-19
   55   neu    0    9   19    0.2669E-22   0.6201E-22   0.1298E-21   0.1007E-21   0.1307E-21   0.6315E-21   0.3651E-19   0.1888E-20
   56   neu    5    0    1     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   57   neu    4    2    3     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   58   neu    4    2    5     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   59   neu    3    4    7     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   60   neu    3    4    9     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   61   neu    2    6   11     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   62   neu    2    6   13     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   63   neu    1    8   15     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   64   neu    1    8   17     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   65   neu    0   10   19     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   66   neu    0   10   21     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   67   neu    5    1    1     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   68   neu    5    1    3     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   69   neu    4    3    5     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   70   neu    4    3    7     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   71   neu    3    5    9     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   72   neu    3    5   11     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   73   neu    2    7   13     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   74   neu    2    7   15     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   75   neu    1    9   17     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   76   neu    1    9   19     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   77   neu    0   11   21     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   78   neu    0   11   23     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   79   neu    6    0    1     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   80   neu    5    2    3     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   81   neu    5    2    5     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   82   neu    4    4    7     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   83   neu    4    4    9     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   84   neu    3    6   11     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   85   neu    3    6   13     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   86   neu    2    8   15     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   87   neu    2    8   17     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   88   neu    1   10   19     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   89   neu    1   10   21     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   90   neu    0   12   23     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
   91   neu    0   12   25     0.000        0.000        0.000        0.000        0.000        0.000        0.000        0.000    
# sum of neutron occ.prob.    5.000        5.000        5.000        5.000        5.000        5.000        5.000        5.000    
 

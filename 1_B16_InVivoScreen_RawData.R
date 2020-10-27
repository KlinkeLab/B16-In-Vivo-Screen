
##NOTE: This script is the first in a series of 2 to complete a Markov Chain analysis on tumor volume growth curves.
#This script organizes the raw data so that it can be easily fed to the Markov Chain function in the next script.


#Basic Clean-up
graphics.off()  #closes all open figures
rm(list=ls())   #clears the variables
cat("\014")     #clears the console

require(stats)
library(EBImage)
library(coda)

#DEFINING RAW GROWTH DATA
#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
#Wild Type data for comparison with DNMT3A and PTPM11 KO data

#Defining the X axis as Days 7 - 19 with data points every 2 days
DP_Days <- seq(7,19,2)

#5 NSG mice injected with B16F0 WT cells as the control group for both DNMT3A and PTPM11 KO groups
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 1, 2, 6, 7, and 11

DPWT_NSG1 <- matrix(c(DP_Days, log10(c(45.5, 270.8, 479.5, 811.9, 982.5, 1596.6, 3339.8))), nrow = 7, byrow = FALSE)
DPWT_NSG2 <- matrix(c(DP_Days, log10(c(41.8, 144.1, 303.9, 592.0, 701.9, 1494.9, 2431.5))), nrow = 7, byrow = FALSE)
DPWT_NSG3 <- matrix(c(DP_Days[1:6], log10(c(23.1, 75.9,  239.3, 367.9, 961.2, 1604.3))), nrow = 6, byrow = FALSE)
DPWT_NSG4 <- matrix(c(DP_Days, log10(c(39.8, 58.8,  201.6, 322.9, 704.5, 1886.9, 3215.2))), nrow = 7, byrow = FALSE)
DPWT_NSG5 <- matrix(c(DP_Days, log10(c(46.8, 83.2,  171.8, 611.6, 1260.5, 1918.0, 2993.9))), nrow = 7, byrow = FALSE)

DPWT_NSG <- abind(DPWT_NSG1, DPWT_NSG2, DPWT_NSG3, DPWT_NSG4, DPWT_NSG5, along = 1)
rm(list = 'DPWT_NSG1', 'DPWT_NSG2', 'DPWT_NSG3', 'DPWT_NSG4', 'DPWT_NSG5')

#Defining X axis as Days 8 - 22 with data points every 2 days
C57_DP_Days  <- seq(8,22,2)

#5 C57BL/6 mice injected with B16FO as the control group for both DNMT3A and PTPN11 KO groups
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 1, 2, 8, 9, and 15
DPWT_C571 <- matrix(c(C57_DP_Days[1:3], log10(c(86.4,	185.2,	547.1))), nrow = 3, byrow = FALSE)
DPWT_C572 <- matrix(c(C57_DP_Days[1:4], log10(c(77.2,	151.5,	471.0,	965.9))), nrow = 4, byrow = FALSE)
DPWT_C573 <- matrix(c(C57_DP_Days[1:6], log10(c(76.3,	150.8,	278.1,	665.9,	544.8,	1453.3))), nrow = 6, byrow = FALSE)
DPWT_C574 <- matrix(c(C57_DP_Days,      log10(c(40.8,	167.7,	411.7,	304.0,	975.7,	883.5,	1669.2,	2417.9))), nrow = 8, byrow = FALSE)
DPWT_C575 <- matrix(c(C57_DP_Days[1:4], log10(c(68.1,	165.6,	355.2,	597.8))), nrow = 4, byrow = FALSE)

DPWT_C57 <- abind(DPWT_C571, DPWT_C572, DPWT_C573, DPWT_C574, DPWT_C575, along = 1)
rm(list = 'DPWT_C571', 'DPWT_C572', 'DPWT_C573', 'DPWT_C574', 'DPWT_C575')

#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
#DNMT3A KO data


#5 NSG mice injected with B16FO DNMT3A KO cells 
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 3, 4, 8, 9, AND 12

D_NSG1 <- matrix(c(DP_Days, log10(c(57.8, 126.7, 254.3, 717.9, 1203.0, 1824.6, 2852.3))), nrow = 7, byrow = FALSE)
D_NSG2 <- matrix(c(DP_Days, log10(c(18.3, 153.9, 313.5, 490.6, 1001.0, 1930.2, 2748.7))), nrow = 7, byrow = FALSE)
D_NSG3 <- matrix(c(DP_Days, log10(c(48.7, 120.6, 369.0, 583.9, 1123.2, 2092.3, 3431.9))), nrow = 7, byrow = FALSE)
D_NSG4 <- matrix(c(DP_Days, log10(c(46.1, 198.5, 410.6, 541.3, 1089.5, 2141.4, 2621.0))), nrow = 7, byrow = FALSE)
D_NSG5 <- matrix(c(DP_Days[1:6], log10(c(72.0, 169.4, 319.8, 567.1, 1108.2, 2621.0))), nrow = 6, byrow = FALSE)

D_NSG <- abind(D_NSG1, D_NSG2, D_NSG3, D_NSG4, D_NSG5, along = 1)
rm(list = 'D_NSG1', 'D_NSG2', 'D_NSG3', 'D_NSG4', 'D_NSG5')


#5 C57B/6 mice injected with B16FO DNMT3A KO cells 
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 3, 4, 10, 11, and 12

#Defining X axis values based on each mouse's survival between Days 8 - 22 with data points every 2 days
D_Days  <- seq(8,20,2)

D_C571 <- matrix(c(D_Days[1:6], log10(c(30.6, 77.5,	  282.8,	474.7,	1463.0,	2293.8))), nrow = 6, byrow = FALSE)
D_C572 <- matrix(c(D_Days[1:4], log10(c(23.9,	89.7,	  398.9,  660.3))), nrow = 4, byrow = FALSE)
D_C573 <- matrix(c(D_Days[1:4], log10(c(38.7,	86.8, 	294.5,	607.0))), nrow = 4, byrow = FALSE)
D_C574 <- matrix(c(D_Days[2:7], log10(c(86.9,	323.7,	472.8,	912.3,	1437.7,	2866.9))), nrow = 6, byrow = FALSE) #note mouse 4 data starts on day 10
D_C575 <- matrix(c(D_Days[1:4], log10(c(36.6,	99.1,	  284.2,	476.7))), nrow = 4, byrow = FALSE)

D_C57 <- abind(D_C571, D_C572, D_C573, D_C574, D_C575, along = 1)
rm(list = "D_C571", "D_C572", "D_C573", "D_C574", "D_C575")

#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
#PTPN11 Data


#5 NSG mice injected with B16F0 PTPN11 KO cells
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 5, 10, 13, 14, AND 15

P_NSG1 <- matrix(c(DP_Days, log10(c(11.7, 51.7, 65.1,  175.3, 374.0, 1168.3, 2161.1))), nrow = 7, byrow = FALSE)
P_NSG2 <- matrix(c(DP_Days, log10(c(32.7, 63.2, 183.3, 227.3, 382.3, 1211.8, 2230.8))), nrow = 7, byrow = FALSE)
P_NSG3 <- matrix(c(DP_Days, log10(c(22.6, 70.0, 139.1, 389.7, 558.4, 1798.8, 3056.0))), nrow = 7, byrow = FALSE)
P_NSG4 <- matrix(c(DP_Days, log10(c(18.9, 43.6, 75.3,  187.7, 463.8, 1327.2, 1502.5))), nrow = 7, byrow = FALSE)
P_NSG5 <- matrix(c(DP_Days, log10(c(14.5, 49.3, 75.9,  171.5, 582.9, 1204.1, 2948.6))), nrow = 7, byrow = FALSE)

P_NSG <- abind(P_NSG1, P_NSG2, P_NSG3, P_NSG4, P_NSG5, along = 1)
rm(list = 'P_NSG1', 'P_NSG2', 'P_NSG3', 'P_NSG4', 'P_NSG5')


#5 C57B/6 mice injected with B16FO PTPN11 KO cells 
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3

#Defining X axis values based on each mouse's survival between Days 10 - 24 with data points every 2 days
P_Days <- seq(10,24,2)

P_C571 <- matrix(c(P_Days[1:5], log10(c(80.2,  	188.5,	306.6,	772.1,	1191.1))), nrow = 5, byrow = FALSE)
P_C572 <- matrix(c(P_Days[1:4], log10(c(58.5,	  126.0,	334.9,	1465.3))), nrow = 4, byrow = FALSE)
P_C573 <- matrix(c(P_Days[1:3], log10(c(76.4, 	115.4,	322.5))), nrow = 3, byrow = FALSE)
P_C574 <- matrix(c(P_Days[2:7], log10(c(29.7,	  110.0,	282.4,	580.6,	1419.7,	2985.2))), nrow = 6, byrow = FALSE) #note mouse 4 data starts on day 12
P_C575 <- matrix(c(P_Days[3:7], log10(c(19.4,	  117.4,	322.7,	1129.6,	1947.4))), nrow = 5, byrow = FALSE) #note mouse 5 data starts on day 14

P_C57 <- abind(P_C571, P_C572, P_C573, P_C574, P_C575, along = 1)
rm(list = 'P_C571', 'P_C572', 'P_C573', 'P_C574', 'P_C575')

#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
#Wild Type data for comparison with SPARC knockout

#Defining the X axis as Days 8 - 16 with data points every 2 days
SWT_Days <- seq(8,16,2)

#5 NSG mice injected with B16F0 WT cells as the control group for both SPARC 1C6 AND SPARC 2C2 KO groups
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 1, 2, 6, 10, and 13

SWT_NSG1 <- matrix(c(SWT_Days, log10(c(16.0, 64.8,  188.7, 441.5,  729.5))), nrow = 5, byrow = FALSE) 
SWT_NSG2 <- matrix(c(SWT_Days, log10(c(26.5, 170.4, 470.4, 1130.0, 1905.9))), nrow = 5, byrow = FALSE)
SWT_NSG3 <- matrix(c(SWT_Days, log10(c(19.8, 87.2,  271.2, 388.6,  583.9))), nrow = 5, byrow = FALSE)
SWT_NSG4 <- matrix(c(SWT_Days, log10(c(61.3, 150.5, 559.9, 1255.9, 1743.2))), nrow = 5, byrow = FALSE)
SWT_NSG5 <- matrix(c(SWT_Days, log10(c(10.1, 42.3,  89.9,  285.7,  378.7))), nrow = 5, byrow = FALSE)

SWT_NSG <- abind(SWT_NSG1, SWT_NSG2, SWT_NSG3, SWT_NSG4, SWT_NSG5, along = 1)
rm(list = 'SWT_NSG1', 'SWT_NSG2', 'SWT_NSG3', 'SWT_NSG4', 'SWT_NSG5')


#5 C57BL/6 mice injected with B16FO as the control group for both SPARC 1C6 AND SPARC 2C2 KO groups
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 1, 2, 6, 7, and 11

#Defining the X axis as Days 7 - 21 with data points every 2 days
C57_SWT_Days <- seq(7,21,2)

SWT_C571 <- matrix(c(C57_SWT_Days,      log10(c(16.8, 42.2, 68.7,  73.1,  105.9, 302.4,  335.5,  699.6))), nrow = 8, byrow = FALSE)
SWT_C572 <- matrix(c(C57_SWT_Days[1:7], log10(c(16.0, 30.1, 49.7,  59.8,  120.6, 172.6,  341.7))), nrow = 7, byrow = FALSE)
SWT_C573 <- matrix(c(C57_SWT_Days,      log10(c(15.9, 49.7, 67.6,  132.8, 169.5, 273.9,  532.1,  1093.3 ))), nrow = 8, byrow = FALSE)
SWT_C574 <- matrix(c(C57_SWT_Days[1:7], log10(c(23.0, 50.9, 141.5, 210.4, 354.7, 635.4,  989.1))), nrow = 7, byrow = FALSE) 
SWT_C575 <- matrix(c(C57_SWT_Days,      log10(c(41.6, 52.1, 111.0, 274.5, 466.3, 1410.7, 2232.8, 2086.2))), nrow = 8, byrow = FALSE) 

SWT_C57 <- abind(SWT_C571, SWT_C572, SWT_C573, SWT_C574, SWT_C575, along = 1)
rm(list = 'SWT_C571', 'SWT_C572', 'SWT_C573', 'SWT_C574', 'SWT_C575')

#___________________________________________________________________________________________________________
#SPARC 1C6 KO data


#5 NSG mice injected with B16F0 SPARC 1C6 KO cells
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 3, 4, 7, 11, and 14

S1C6_NSG1 <- matrix(c(SWT_Days, log10(c(11.7, 46.2, 171.9, 467.6, 659.2))), nrow = 5, byrow = FALSE) 
S1C6_NSG2 <- matrix(c(SWT_Days, log10(c(13.2, 40.2, 166.5, 402.7, 593.4))), nrow = 5, byrow = FALSE) 
S1C6_NSG3 <- matrix(c(SWT_Days, log10(c(18.1, 73.6, 249.7, 435.0, 694.4))), nrow = 5, byrow = FALSE) 
S1C6_NSG4 <- matrix(c(SWT_Days[2:5], log10(c(34.8, 179.2, 297.6, 479.5))), nrow = 4, byrow = FALSE) 
S1C6_NSG5 <- matrix(c(SWT_Days, log10(c(9.9,  46.1, 106.9, 304.1, 412.2))), nrow = 5, byrow = FALSE) 

S1C6_NSG <- abind(S1C6_NSG1, S1C6_NSG2, S1C6_NSG3, S1C6_NSG4, S1C6_NSG5, along = 1 )
rm(list = 'S1C6_NSG1', 'S1C6_NSG2', 'S1C6_NSG3', 'S1C6_NSG4', 'S1C6_NSG5' )


#5 C57B/6 mice injected with B16F0 SPARC 1C6 KO cells
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 3, 4, 8, 9, and 12

S1C6_C571 <- matrix(c(C57_SWT_Days,      log10(c(29.5, 42.5, 57.6, 71.0,  145.5, 338.0, 525.3,  1669.2))), nrow = 8, byrow = FALSE)
S1C6_C572 <- matrix(c(C57_SWT_Days[1:6], log10(c(21.2, 38.2, 87.2, 110.6, 229.4, 627.3))), nrow = 6, byrow = FALSE)
S1C6_C573 <- matrix(c(C57_SWT_Days,      log10(c(8.6,  24.2, 14.0, 64.1,  97.1,  144.6, 191.6,  457.4))), nrow = 8, byrow = FALSE)
S1C6_C574 <- matrix(c(C57_SWT_Days,      log10(c(10.8, 34.0, 26.5, 107.7, 169.4, 279.4, 353.8,  730.6))), nrow = 8, byrow = FALSE)
S1C6_C575 <- matrix(c(C57_SWT_Days,      log10(c(20.5, 34.8, 69.4, 89.0,  422.1, 838.2, 1473.6, 2968.4))), nrow = 8, byrow = FALSE)

S1C6_C57 <- abind(S1C6_C571, S1C6_C572, S1C6_C573, S1C6_C574, S1C6_C575, along = 1 )
rm(list = 'S1C6_C571', 'S1C6_C572', 'S1C6_C573', 'S1C6_C574', 'S1C6_C575' )

#___________________________________________________________________________________________________________
#SPARC 2C2 KO data

#5 NSG mice injected with B16F0 SPARC 2C2 KO cells
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
S2C2_NSG1 <- matrix(c(SWT_Days,      log10(c(16.5,  82.4,  195.2, 423.7, 650.3))), nrow = 5, byrow = FALSE) 
S2C2_NSG2 <- matrix(c(SWT_Days[2:5], log10(c(46.9,  187.0, 359.1, 665.9))), nrow = 4, byrow = FALSE) 
S2C2_NSG3 <- matrix(c(SWT_Days,      log10(c(15.0,  73.5,  255.1, 559.3, 729.5))), nrow = 5, byrow = FALSE) 
S2C2_NSG4 <- matrix(c(SWT_Days,      log10(c(31.8,  124.8, 284.6, 486.6, 793.9))), nrow = 5, byrow = FALSE) 
S2C2_NSG5 <- matrix(c(SWT_Days,      log10(c(11.1,  26.9,  75.4,  287.7, 457.0))), nrow = 5, byrow = FALSE) 

S2C2_NSG <- abind(S2C2_NSG1, S2C2_NSG2, S2C2_NSG3, S2C2_NSG4, S2C2_NSG5, along = 1 )
rm(list = 'S2C2_NSG1', 'S2C2_NSG2', 'S2C2_NSG3', 'S2C2_NSG4', 'S2C2_NSG5' )


#5 C57B/6 mice injected with B16F0 SPARC 1C6 KO cells
#Data vector is PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 3, 4, 8, 9, and 12

S2C2_C571 <- matrix(c(C57_SWT_Days[2:8], log10(c(12.7,  27.1,  53.2,  155.0, 402.1, 641.8,  1386.7))), nrow = 7, byrow = FALSE)
S2C2_C572 <- matrix(c(C57_SWT_Days[2:8], log10(c(16.5,  32.3,  145.5, 228.0, 682.4, 1296.8, 1665.6))), nrow = 7, byrow = FALSE)
S2C2_C573 <- matrix(c(C57_SWT_Days[2:8], log10(c(12.6,  172.6, 254.6, 289.9, 595.4, 861.6,  2078.3))), nrow = 7, byrow = FALSE)
S2C2_C574 <- matrix(c(C57_SWT_Days[4:8], log10(c(27.4,  57.3,  74.8,  227.6, 303.4))), nrow = 5, byrow = FALSE)
S2C2_C575 <- matrix(c(C57_SWT_Days[3:7], log10(c(188.8, 194.7, 293.0, 484.4, 636.6))), nrow = 5, byrow = FALSE)

S2C2_C57 <- abind(S2C2_C571, S2C2_C572, S2C2_C573, S2C2_C574, S2C2_C575, along = 1 )
rm(list = 'S2C2_C571', 'S2C2_C572', 'S2C2_C573', 'S2C2_C574', 'S2C2_C575' )

#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
#Wild Type data for comparison with WISP1 knockout

#Establishing X values for WT data to compare to WISP1 KO
WWT_Days <- seq(7,19,2)
C57_WWT_Days <- c(7, 10, 15, 16, 17, 18, 19, 20, 21)
C57_WWT_Days2 <- c(7, 9, 14, 16, 20, 21)
#5 NSG mice injected with B16F0 WT cells as the control group for the WISP1 KO group
#Data rows are PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 1, 3, 4, 7, and 8

WWT_NSG1 <- matrix(c(WWT_Days[1:5], log10(c(38.7, 84.4,  701.6, 1400.9, 2851.0))), nrow = 5, byrow = FALSE)
WWT_NSG2 <- matrix(c(WWT_Days,      log10(c(38.5, 96.9,  251.0, 373.0,  508.9,  1789.2, 3202.2))), nrow = 7, byrow = FALSE)
WWT_NSG3 <- matrix(c(WWT_Days[1:6], log10(c(49.9, 142.2, 443.9, 835.4,  1565.0, 2561.2))), nrow = 6, byrow = FALSE)
WWT_NSG4 <- matrix(c(WWT_Days,      log10(c(45.2, 103.8, 177.1, 356.3,  530.1,  1424.1, 3208.0))), nrow = 7, byrow = FALSE)
WWT_NSG5 <- matrix(c(WWT_Days[1:6], log10(c(75.3, 144.4, 233.4, 510.9,  1382.6, 2021.4))), nrow = 6, byrow = FALSE)

WWT_NSG <- abind(WWT_NSG1, WWT_NSG2, WWT_NSG3, WWT_NSG4, WWT_NSG5, along = 1 )
rm(list = 'WWT_NSG1', 'WWT_NSG2', 'WWT_NSG3', 'WWT_NSG4', 'WWT_NSG5' )
                     
#5 C57B/6 mice injected with B16F0 WT cells as the control group for the WISP1 KO group
#Data rows are PER MOUSE, each point is an individual mouse tumor volume in mm3
#Data from mice 1, 2, 3, 6, AND 7

WWT_C571 <- matrix(c(C57_WWT_Days[1:4],        log10(c(99.4,  151.3, 312.5,  579.7))), nrow = 4, byrow = FALSE)
WWT_C572 <- matrix(c(C57_WWT_Days[1:3],        log10(c(129.4, 178.6, 1231.7))), nrow = 3, byrow = FALSE)
WWT_C573 <- matrix(c(C57_WWT_Days[c(1:3,5:9)], log10(c(84.9,  218.2, 284.3,  549.8,  684.0, 1224.5, 1537.7, 1889.2))), nrow = 8, byrow = FALSE)
WWT_C574 <- matrix(c(C57_WWT_Days[1:3],        log10(c(63.6,  201.6, 321.1))), nrow = 3, byrow = FALSE)
WWT_C575 <- matrix(c(C57_WWT_Days[c(1:3,5:9)], log10(c(103.0, 119.1, 329.1,  559.9,  675.4, 1172.1, 1352.2, 1168.7))), nrow = 8, byrow = FALSE)
WWT_C576 <- matrix(c(C57_WWT_Days2[1:5],       log10(c(30.22, 34.38, 248.33, 352.83, 564.5))), nrow = 5, byrow = FALSE)
WWT_C577 <- matrix(c(C57_WWT_Days2[c(1:4,6)],  log10(c(39.02, 44.01, 283.27, 627.31, 2029.89))), nrow = 5, byrow = FALSE)
WWT_C578 <- matrix(c(C57_WWT_Days2[1:3],       log10(c(36.56, 88.8,  466.3))), nrow = 3 , byrow = FALSE)

WWT_C57 <- abind(WWT_C571, WWT_C572, WWT_C573, WWT_C574, WWT_C575, WWT_C576, WWT_C577, WWT_C578, along = 1)
rm(list = "WWT_C571", "WWT_C572", "WWT_C573", "WWT_C574", "WWT_C575", "WWT_C576", "WWT_C577", 'WWT_C578')
#___________________________________________________________________________________________________________
#WISP1 KO data

#5 NSG mice injected with B16F0 WISP1 KO cells
#Data rows are PER DAY, each point is an individual mouse tumor volume in mm3
#Data from mice 2, 5, 6, 9, and 10

W_NSG1 <- matrix(c(WWT_Days[1:5], log10(c(57.0, 123.5, 608.0, 1094.1, 2926.1))), nrow = 5, byrow = FALSE)
W_NSG2 <- matrix(c(WWT_Days[1:6], log10(c(62.7, 196.1, 449.2, 838.9,  1415.5, 2377.5))), nrow = 6, byrow = FALSE)
W_NSG3 <- matrix(c(WWT_Days[1:6], log10(c(52.7, 140.3, 551.4, 903.9,  1905.9, 2940.7))), nrow = 6, byrow = FALSE)
W_NSG4 <- matrix(c(WWT_Days[1:6], log10(c(70.9, 107.6, 286.7, 794.5,  2289.8, 3409.6))), nrow = 6, byrow = FALSE)
W_NSG5 <- matrix(c(WWT_Days[1:6], log10(c(52.3, 119.1, 266.0, 714.2,  1956.3, 2993.9))), nrow = 6, byrow = FALSE)

W_NSG <- abind(W_NSG1, W_NSG2, W_NSG3, W_NSG4, W_NSG5, along = 1)
rm(list = "W_NSG1", "W_NSG2", "W_NSG3", "W_NSG4", "W_NSG5")

#Establishing X values for WISP1 KO
C57_W_Days <- c(7, 10, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25)
C57_W_Days2 <- c(7, 9, 14, 16, 21, 23, 25, 26, 30)

#5 C57B/6 mice injected with B16F0 WISP1 KO cells
#Data rows are PER DAY, each point is an individual mouse tumor volume in mm3
#Data from mice 4, 5, 8, 9, and 10

W_C571 <- matrix(c(C57_W_Days[1:7],  log10(c(58.2, 205.9, 318.1,  729.2,  1093.6, 1032.2, 1139.6))), nrow = 7, byrow = FALSE)
W_C572 <- matrix(c(C57_W_Days,       log10(c(38.0, 64.7,  35.1,   267.6,  355.2,  362.9,  483.5,   773.5, 761.3, 838.6, 1096.6, 1240.0))), nrow = 12, byrow = FALSE)
W_C573 <- matrix(c(C57_W_Days[1:8],  log10(c(70.9, 77.6,  366.2,  586.0,  773.3,  1246.9, 902.0,   816.8))), nrow = 8, byrow = FALSE)
W_C574 <- matrix(c(C57_W_Days[1:4],  log10(c(31.2, 70.0,  350.8,  514.1))), nrow = 4, byrow = FALSE)
W_C575 <- matrix(c(C57_W_Days,       log10(c(29.2, 46.2,  98.5,   236.9,  535.2,  465.5,  602.3,   603.4, 733.2, 966.3, 1525.6, 1510.8))), nrow = 12, byrow = FALSE)
W_C576 <- matrix(c(C57_W_Days2[1:8], log10(c(17.6, 4.19,  113,    130.06, 261.04, 365.99, 1592.79, 1592.79))), nrow = 8, byrow = FALSE)
W_C577 <- matrix(c(C57_W_Days2[1:6], log10(c(4.0,  96.4,  598.41, 608.43, 680.68, 2852.17))), nrow = 6, byrow = FALSE)
W_C578 <- matrix(c(C57_W_Days2,      log10(c(16.8, 14.26, 144.2,  89.01,  184.47, 163.03, 335.1,   335.1, 2052.5))), nrow = 9, byrow = FALSE)

W_C57 <- abind(W_C571, W_C572, W_C573, W_C574, W_C575, W_C576, W_C577, W_C578, along = 1)
rm(list = "W_C571", "W_C572", "W_C573", "W_C574", "W_C575", "W_C576", "W_C577", "W_C578")

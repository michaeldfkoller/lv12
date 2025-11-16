import markovlv
from markovlv import MARKOVLV
from pylab import *
import numpy 
from scipy import *
import scipy.io

def read_array(filename, mylen):
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
    """
    return (numpy.loadtxt(filename))

disp('Calculating LTC')
# ============================================================
iVon  = input("Start Alter (40) ") 
iNach = input("End Alter (100) ")
iProduceAudit = 0
#symTex = open("test.tex","w")
# ============================================================
pij = read_array('pij.dat', 4)
post = read_array('post.dat', 4)
pre = read_array('pre.dat', 4)
PrePrem1 = read_array('PrePrem1.dat',4)
zins = read_array("zins.dat", 2)
# ============================================================
fid = open('MResult.dat','w');
pHandlePrem  = MARKOVLV(122,10,2)
pHandleBenefits  = MARKOVLV(122,10,2)
pHandlePrem.vSetNrStates(9)
pHandleBenefits.vSetNrStates(9)
# ============================================================
pHandlePrem.vSetStartTime(120)
pHandlePrem.vSetStopTime(20)
pHandleBenefits.vSetStartTime(120)
pHandleBenefits.vSetStopTime(20)
# ============================================================

disp('Init pij')
for i in  range(pij.shape[0]):
    dTemp= pHandlePrem.dSetPij(long(pij[i,0]), long(pij[i,1]),long(pij[i,2]), pij[i,3])
#    print "Pij  %d %d %d %10.4f" %(long(0.005+pij[i,0]), long(0.005+pij[i,1]-1),long(0.005+pij[i,2]-1), dTemp)
    dTemp= pHandleBenefits.dSetPij(long(pij[i,0]), long(pij[i,1]),long(pij[i,2]), pij[i,3])
#    print "Pij2 %d %d %d %10.4f" % (long(0.005+pij[i,0]), long(0.005+pij[i,1]),long(0.005+pij[i,3]), pij[i,3])

disp('Init post')

#for i in  range(post.shape[0]):
dTemp= pHandlePrem.dSetPost(0,0,0,0)
dTemp= pHandleBenefits.dSetPost(0,0,0,0)

disp('Init pre')

for i in  range(101):
  for j in range(4,8):
    dTemp= pHandleBenefits.dSetPre(i, j, j, 12000.)
#    print "PRE %d %d  %d %10.4f" %(long(0.005+pre[i,0]), long(0.005+pre[i,1]-1),long(0.005+pre[i,2]-1), dTemp)

for i in  range(PrePrem1.shape[0]):
    dTemp= pHandlePrem.dSetPre(long(PrePrem1[i,0]), long(PrePrem1[i,1]),long(PrePrem1[i,2]), PrePrem1[i,3])
#    print "PRE2 %d %10.4f" %(long(0.005+pre[i,0]), dTemp)

for i in  range(zins.shape[0]):
    for j in range(9):
      pHandlePrem.dSetDisc(long(zins[i,0]), j,j, zins[i,1])
      pHandleBenefits.dSetDisc(long(zins[i,0]), j,j, zins[i,1])
#      print long(0.005+zins[i,0]), j,j, zins[i,1]
# ============================================================
values = zeros((100,4), dtype='f')
j= 0
for i in range(iVon, iNach+1):
    P_1 = pHandleBenefits.dGetDK(i,1,1) /  pHandlePrem.dGetDK(i,1,1) + 0.000001
#    print  i, pHandleBenefits.dGetDK(i,1,1), pHandlePrem.dGetDK(i,1,1), pHandleBenefits.dGetDK(i,4,1)
    P_2 = pHandleBenefits.dGetDK(i,2,1) /  pHandlePrem.dGetDK(i,2,1)
    strTemp = 'x  %5d   P1  %10.5f    P2  %10.5f     Verh %5.2f %%'%(i,P_1, P_2, 100 * P_2 / P_1)
    disp(strTemp)
    fid.write(strTemp+"\n")
    values[j,0]=  i
    values[j,1]= P_1
    values[j,2]= P_2
    values[j,3]= P_2/P_1
    j = j + 1
#pHandleBenefits.vPrintTeX(symTex,True,"Test", True)
#symTex.close()

# ============================================================
valuesP = zeros((100,8), dtype='f') 
valuesB = zeros((100,8), dtype='f')
j2= 0
for i in range (iVon, 101):
  j2 = j2+1
  for k in range(8):
      valuesP[j2,k] = pHandlePrem.dGetDK(i,k,1)
      valuesB[j2,k] = pHandleBenefits.dGetDK(i,k,1)

# ============================================================
figure(1)
subplot(2,1,1)
plot(values[0:j-1,0],values[0:j-1,1],'-b', values[0:j-1,0],values[0:j-1,2],'-.ro')
#Title('Absolute Premium for 12k / 6k')
subplot(2,1,2)
plot(values[0:j-1,0],values[0:j-1,3],'-dg')
#Title('Relative Premium P_2 / P_1')
fid.close()
# ============================================================
show()
# ============================================================



import os
import sys
import math
import string
import math
import time
from scipy import *
import numpy
import weave

class StopLoss:

    def __init__(self):
        self.QxFile = "qx.dat"
        print 'QxFile:' , self.QxFile
        print("Plese enter new value:")
        strInput=sys.stdin.readline()
        strInput=strInput.strip() # Newline (\n) am Ende der Eingabe entfernen
        if strInput == "":
            print("Using above mentioned Parameter")
        else:
            self.QxFile = strInput

        self.InforceFile = "pk.dat"
        print 'Inforce:' , self.InforceFile
        print("Plese enter new value:")
        strInput=sys.stdin.readline()
        strInput=strInput.strip() # Newline (\n) am Ende der Eingabe entfernen
        if strInput == "":
            print("Using above mentioned Parameter")
        else:
            self.InforceFile = strInput

        self.ResultFile = "result.dat"
        print 'Result:' , self.ResultFile
        print("Plese enter new value:")
        strInput=sys.stdin.readline()
        strInput=strInput.strip() # Newline (\n) am Ende der Eingabe entfernen
        if strInput == "":
            print("Using above mentioned Parameter")
        else:
            self.ResultFile = strInput

        self.delta = 2500.0
        self.nrclasses = 2**20
        print ' Parameters \n ---------- \n'#            print alter, sex
        strTemp = " QxFile: %20s \n InForce %20s \n Delta : %20.0f \n Classes %20d \n Max Loss %20.4e \n \n" % (self.QxFile, self.InforceFile, self.delta, self.nrclasses,self.delta*self.nrclasses)
        print strTemp
        print 'Loading Data'#            print alter, sex
        self.qx = scipy.io.read_array(self.QxFile)
        self.pk = scipy.io.read_array(self.InforceFile)
#        print self.qx
#        print self.pk
        self.Now = time.clock()

    def calclambda(self):
        # Funktion berechnet Paprameter Lambda for compoud poisson
        # use [lambda] = calclambda(pk,qx)
        actlam = 0
        for i in range(self.pk.shape[0]):
            alter = self.pk[i,0]
            sex   = self.pk[i,1]
            qxakt = self.qx[alter,sex]
            actlam = actlam + qxakt
        return(actlam) 

    def calcpois(self, ev, lam):
    #
    # berechnet gesamtschadenverteilung 
    # use [gv]= calcpois(ev, lambda, self.nrclasses)
    #
        gvfft = fft(ev)
        gvfft = exp(lam *(gvfft - 1))
        gv = real(ifft(gvfft))
        return(gv)

    def makeeinzelvert(self):
    #
    # Berechnet Einzzelschadenverteilung fuer pk
    # use  [ev] = makeeinzelvert(pk,self.qx,self.delta,self.nrclasses)
    #
        temp   = zeros([self.nrclasses], dtype='f')
        lam = 0.0
        for i in range(self.pk.shape[0]):
            alter  = self.pk[i,0]
            sex    = self.pk[i,1]
            qxakt  = self.qx[alter,sex]
            lam +=  qxakt
            leist  = self.pk[i,2]
            index  = round(leist / self.delta )
            temp[index] += qxakt
        return(temp / lam)

    def dojob(self):
        lam = self.calclambda()
        print "Lambda", lam
        ev = self.makeeinzelvert()
#        fid2 = open("temp.dat","w")
#        scipy.io.write_array(fid2,  ev)
#        fid2.close()        
        gv = self.calcpois(ev, lam)
#        fid2 = open("temp.dat","w")
#        scipy.io.write_array(fid2,  gv)
#        fid2.close()
        print 'Elased Time (FFT+Reading):', time.clock()-self.Now


        fid = open(self.ResultFile,'w')

        temp = 0.

        for i in range(self.pk.shape[0]):
            alter = self.pk[i,0]
            sex   = self.pk[i,1]
            qxakt = self.qx[alter,sex]
            leist = self.pk[i,2]
            temp += leist * qxakt
        ExpLoss = temp
        strOut = ""
        strLine = "\n Pensionskasse: N " + str(self.pk.shape[0])
        print strLine
        strOut += strLine
        strLine = "\n Erwartungswert: " + str(temp)
        print strLine
        strOut += strLine

        a= zeros([4,100], dtype='f')

        cumgv = gv.copy()
#         SLP = gv.copy()
        for i in range(1,gv.shape[0]):
             cumgv[i] = cumgv[i-1] + gv[i]
#         MyLoss = self.delta * ( gv.shape[0] - 1.)   
#         SLP[i] = cumgv[i-1] + gv[i]
#         for i in range(gv.shape[0],1,-1):


        
        for i in range(80):
            a[0,i] = i * 0.05
        for i in range(80):    
            index  = int(round(a[0,i] * temp / self.delta ))
            p = 1. -cumgv[index]
            e = 0
#            print gv[range(index,self.nrclasses)].shape[0],gv[range(index,self.nrclasses)].shape[1]
            expression ="zz = array(range(index, self.nrclasses)).transpose()"
            weave.blitz(expression, compiler_name = 'mingw32')
            zz = matrix(zz)
#            print zz.shape[0], zz.shape[1]
#            print zz * gv[range(index,self.nrclasses)]
            expression = "e = self.delta * zz * gv[range(index,self.nrclasses)]"
            weave.blitz(expression, compiler_name = 'mingw32')
            a[1,i] = index
            a[2,i] = p
            a[3,i]= e
            strLine = "\n " + str(a[0,i]) + " ES  Index "+ str(a[1,i]) +" Prob " + str(a[2,i]) +" E((S-X)+) " + str(a[3,i])
#            print strLine
            strOut += strLine
        fid.write(strOut)
        fid.close()
        index2 =  int(round(79 * 0.05 * ExpLoss / self.delta))
        aa= array(range(1,index2)) * self.delta
        myfigure =  figure(2)
        subplot(211)
        plot(aa,gv[range(1,index2)],'r',label="Composite Distribution")
        plot(aa,ev[range(1,index2)]/100,'b:', label ="Individual Distribution")
        title('Distribution')
        legend()
        ab= a[0,range(80)] * ExpLoss
        ac= 1 - a[2,range(80)]
        subplot(212)
        plot(ab, ac,'g', label="Probability", linewidth = 6)
        grid(b=1,color="r")
        [XMin, XMax, YMin, YMax] = axis()
        [XMin, XMax, YMin, YMax] = axis([XMin, XMax, 0, 1 ])
        title('CdF')
        for i in range(7):
            LocIndex = 10 * (i+1)
            xPos = 0.5*(XMax - XMin) + XMin
            yPos = ((6-i)+2)/10.*(YMax-YMin) + YMin
            if i == 6:
                MyText =  "1.0 x ES = %10.0f"  % (ExpLoss)
            else:
                MyText = str(a[0,LocIndex]) + " x ES --> Prob: %6.4f %20s %10.0f"  % (1-round(1e4*a[2,LocIndex])/1e4 ,'E[(S-X)+]:', round(a[3,LocIndex]))  
            text(xPos, yPos, MyText)
        print("Plese Adjust Graphic")
        strInput=sys.stdin.readline()
            
        savefig("slp.ps", dpi = 300, orientation ="landscape", papertype = "a4")
        print 'Elased Time:', time.clock()-self.Now

def main():
#############################################################################################################    
        slp = StopLoss()
        slp.dojob()

if __name__ == "__main__":
    main()
    


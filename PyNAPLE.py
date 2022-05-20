#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 17:25:21 2022

@author: daniel
"""

#import glob
#import tk_tools

dependencyList = []

try:
    import os
except ImportError:
    dependencyList.append("os")

try:
    import re
except ImportError:
    dependencyList.append("re")

try:
    import sys
except ImportError:
    dependencyList.append("sys")

try:
    import time
except ImportError:
    dependencyList.append("time")

try:
    import platform
except ImportError:
    dependencyList.append("platform")
#import numpy as np
#from pyds9 import *

try:
    import tkinter as tk
    from tkinter import filedialog as fd
    from tkinter import simpledialog as sd    
    from tkinter import scrolledtext as st
except ImportError:
    dependencyList.append("tkinter")

try:
    import subprocess as sp
except ImportError:
    dependencyList.append("subprocess")

try:
    import mechanize as mech
except ImportError:
    dependencyList.append("mechanize")
#from astropy.io import fits
#import multiprocessing as mp

try:
    from bs4 import BeautifulSoup
except ImportError:
    dependencyList.append("bs4")

try:
    import matplotlib.pyplot as plt
except ImportError:
    dependencyList.append("matplotlib")
#from matplotlib.image import imread

try:
    from tkcalendar import DateEntry
except ImportError:
    dependencyList.append("tkcalendar")

try:
    from shapely.geometry import MultiPoint, Polygon
except ImportError:
    dependencyList.append("shapely.geometry")

try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

if platform.system() == "Darwin":   ### if its a Mac
    from tkmacosx import Button
else:
    from tkinter import Button
   
if len(dependencyList):
    depend = open("DEPENDENCIES", "w")
    for i in dependencyList:
        depend.write(i)
    depend.close()
    print("Imports not found")
    print("Writing missing dependencies to DEPENDENCIES")
    quit()

def main():
    
    if len(sys.argv) == 1:
        root = tk.Tk()
        root.resizable(width=False, height=False)
        root.protocol("WM_DELETE_WINDOW", lambda arg=root: on_closing(arg))
        app = PyNAPLEGUI(master=root)
        app.mainloop()
    elif "-h".casefold() in sys.argv or "--help".casefold() in sys.argv:
        runHelp()
        return
    else:
        dummy = fakeGUI()
        retcode = dummy.handleArgs([a.casefold() for a in sys.argv[1:]])
        if retcode:
            return

class PyNAPLE:
    
    def __init__(self, app):
        
        
        def dirDecode(LAT, LONG):
            """
            Takes Latitude and Longitude and generates a string for the location.
            I.E. Lat: 4.5, Long: -12 would become N4.5W12.
            """
            while LONG > 180:
              LONG -= 360
            if (LAT >= 0):
                latDir = "N"+str(LAT)
            else:
                latDir = "S"+str(abs(LAT))
            if (LONG >= 0):
                longDir = "E"+str(LONG)
            else:
                longDir = "W"+str(abs(LONG))
        
            return latDir+longDir
        
        self.app = app
        self.directory = app.directory.get()
        epochStr = str(app.epoch).replace('-', '')
        localStr = dirDecode(self.app.latitude.get(), self.app.longitude.get())
        self.eventDir = str(app.directory.get())+"/processed/"+epochStr+"_"+localStr
        self.pynaple()

    def pynaple(self):
        print("Now Running...")
        sp.check_call(["mkdir -p "+ self.eventDir], shell=True, executable="/bin/bash", 
                      cwd=self.directory, stdout=sys.stdout, stderr=sys.stderr)
        self.makePlot()
        self.processedPairs = []
        while self.app.searchSpace.get() >= self.app.currentSearchSpace.get():
            print()
            preNames, postNames = self.findImages()
            preImages, postImages = self.getData(preNames, postNames)
            if len(preImages) and len(postImages):
                pairlist = self.compare(preImages, postImages)
                self.app.pairsFound.set(len(pairlist))
                for pair in pairlist:
                    if self.pairIsNew(pair, self.processedPairs):
                        success = self.process(*pair)
                        self.processedPairs.append(pair)
                        self.app.pairsProcessed.set(len(self.processedPairs))
                        if success:
                            pass
                            #TODO add change evaluation
                
            if self.app.currentSearchSpace.get() == self.app.searchSpace.get():
                print()
                print("Process Complete.")
                break
            print("Expanding search space...")
            self.app.currentSearchSpace.set(self.app.currentSearchSpace.get()+self.app.searchIncrement.get())
            
    def pairIsNew(self, pair, pairlist):
        for item in pairlist:
             if pair[0].name == item[0].name and pair[1].name == item[1].name:
                 return False
        return True
    
    def makePlot(self):
        
        self.fig, self.ax = plt.subplots(1,1)
        plt.xlim(self.app.longitude.get() - (0.5 + (self.app.searchSpace.get()/2.)), self.app.longitude.get() + (0.5 + (self.app.searchSpace.get()/2.)))
        plt.ylim(self.app.latitude.get() - (0.5 + (self.app.searchSpace.get()/2.)), self.app.latitude.get() + (0.5 + (self.app.searchSpace.get()/2.)))
        self.ax.scatter(self.app.longitude.get(), self.app.latitude.get(), marker='x', label='Center')
        plt.xlabel("Selenographic Longitude")
        plt.ylabel("Selenographic Latitude")
        searchBox = Polygon(MultiPoint([[self.app.longitude.get()+self.app.searchSpace.get()/2., self.app.latitude.get()-self.app.searchSpace.get()/2.],
                                        [self.app.longitude.get()-self.app.searchSpace.get()/2., self.app.latitude.get()-self.app.searchSpace.get()/2.],
                                        [self.app.longitude.get()-self.app.searchSpace.get()/2., self.app.latitude.get()+self.app.searchSpace.get()/2.],
                                        [self.app.longitude.get()+self.app.searchSpace.get()/2., self.app.latitude.get()+self.app.searchSpace.get()/2.],
                                        [self.app.longitude.get()+self.app.searchSpace.get()/2., self.app.latitude.get()-self.app.searchSpace.get()/2.]]))
        
        self.ax.plot(*searchBox.exterior.xy, linestyle='dotted', label="Search Space")
        plt.legend(loc='upper left')
        self.fig.savefig(self.eventDir+"/location.png")
    
    def createDirectory(self, preimg, postimg):
        
        dirString = self.eventDir+"/"+preimg.name+"_"+postimg.name
        sp.check_call(["mkdir -p "+ dirString], shell=True, 
                      executable="/bin/bash", cwd=self.directory, 
                      stdout=sys.stdout, stderr=sys.stderr)
        return dirString
        
    def process(self, preimg, postimg):
        
        self.app.currentPreImage.set(preimg.name)
        self.app.currentPostImage.set(postimg.name)
        
        print("Temporal pair: "+preimg.name+" & "+postimg.name)
        
        
#        preThread = mp.Process(target=imageDownloadHelper, args=[preimg])
#        #postThread = mp.Process(target=imageDownloadHelper, args=[postimg])
#        preThread.start()
#        #postThread.start()
#        imageDownloadHelper(postimg)
#        
#        preThread.join()
#        #postThread.join()
        
        preimg.downloadImage()
        postimg.downloadImage()
        
        self.plotPair(preimg, postimg)
        pairDir = self.createDirectory(preimg, postimg)
        preImage = preimg.name
        postImage = postimg.name
        proceed = False
        finished = False
        fileLoc = pairDir+"/deffile"
        sp.check_call(["touch "+fileLoc], shell=True, executable="/bin/bash", 
                      cwd=self.directory, stdout=sys.stdout, stderr=sys.stderr)
        with open(fileLoc, 'r') as infile:
            checkpoint = infile.readline()
        if checkpoint == '':
            checkpoint = 0
        checkpoint = int(checkpoint)
        cpReached = 0
        if checkpoint == -1:
            print("Pair already processed")
            return
        
        print("Image processing will now begin...")
        if checkpoint > 0:
            print("Resuming from last checkpoint")
    
    
        isis_operations = ["cam2map from="+self.directory+"/img/"+preImage+".cal.cub to="+preImage+".cal.map.cub",
                           "cam2map from="+self.directory+"/img/"+postImage+".cal.cub map="+preImage+".cal.map.cub to="+postImage+".cal.map.cub matchmap=yes",
                           "tonematch from="+postImage+".cal.map.cub match="+preImage+".cal.map.cub to="+postImage+".cal.map.eq.cub",
                           "autoregtemplate algorithm="+self.app.algorithm.get()+" tolerance="+str(self.app.tolerance.get())+" psamp=15 pline=15 ssamp=215 sline=215 topvl=autoRegTemp_"+preImage+"_"+postImage+".def",
                           "coreg from="+postImage+".cal.map.eq.cub match="+preImage+".cal.map.cub to="+postImage+".cal.map.cub deffile=autoRegTemp_"+preImage+"_"+postImage+".def onet=controlNetwork_"+preImage+"_"+postImage+".net",
                           "cnetedit cnet=controlNetwork_"+preImage+"_"+postImage+".net onet=controlNetwork_"+preImage+"_"+postImage+".net",
                           #"warp from="+postImage+".cal.map.eq.cub to="+postImage+".cal.map.warp.cub cnet=controlNetwork_"+preImage+"_"+postImage+".net cube="+preImage+".cal.map.cub weighted=no degree=2",
                           #'printf "'+preImage+'.cal.map.cub\n'+postImage+'.cal.map.warp.cub\n" > imList.lis',
                           'printf "'+preImage+'.cal.map.cub\n'+postImage+'.cal.map.cub\n" > imList.lis',
                           "cubeit fromlist=imList.lis to=duoCube.cub",
                           "bandtrim from=duoCube.cub to=duoCube.trim.cub",
                           "cropspecial from=duoCube.trim.cub to=duoCube.crop.trim.cub",
                           "explode from=duoCube.crop.trim.cub to=duoCube",
                           "autoregtemplate algorithm="+self.app.algorithm.get()+" tolerance=0.85 psamp=15 pline=15 ssamp=55 sline=55 topvl=autoRegTemp_"+preImage+"_"+postImage+".def",
                           "coreg from=duoCube.band0002.cub match=duoCube.band0001.cub deffile=autoRegTemp_"+preImage+"_"+postImage+".def onet=controlNetwork_"+preImage+"_"+postImage+"_V2.net rows=1000 columns=500", #was 500 150
                           "cnetedit cnet=controlNetwork_"+preImage+"_"+postImage+"_V2.net onet=controlNetwork_"+preImage+"_"+postImage+"_V2.net",
                           "warp from=duoCube.band0002.cub to="+postImage+".warp.crop.cal.map.cub cnet=controlNetwork_"+preImage+"_"+postImage+"_V2.net cube=duoCube.band0001.cub weighted=no degree=2",
                           "tonematch from="+postImage+".warp.crop.cal.map.cub match=duoCube.band0001.cub to="+postImage+".eq.warp.crop.cal.map.cub",
                           "ratio numerator="+postImage+".eq.warp.crop.cal.map.cub denominator=duoCube.band0001.cub to="+postImage+"_over_"+preImage+".cub",
                           "ratio numerator=duoCube.band0001.cub denominator="+postImage+".eq.warp.crop.cal.map.cub to="+preImage+"_over_"+postImage+".cub",
                           "mv duoCube.band0001.cub "+preImage+".crop.cal.map.cub"]
    #                       "isis2std from="+preImage+".cal.map.cub to="+preImage+".png",
    #                       "isis2std from="+postImage+".cal.map.warp_V2.cub to="+postImage+".png",
    #                       "isis2std from="+postImage+"_over_"+preImage+".cub to="+postImage+"_over_"+preImage+".png
    
        for operation in isis_operations:
            if not proceed:
                if cpReached >= checkpoint:
                    print(time.ctime().split()[3]+": Now running " + operation.split()[0])
                    proceed = sp.check_call([". ./../../../runIsis.sh && "+operation],
                                            shell=True, executable="/bin/bash",
                                            cwd=pairDir, stdout=DEVNULL,
                                            stderr=sys.stdout)
    
                cpFile = open(pairDir+"/deffile", "w+")
                cpReached = cpReached + 1
                cpFile.write(str(cpReached))
                cpFile.close()
    
        if not proceed:
            if cpReached >= checkpoint:
                removal = [postImage+".cal.map.cub",
                           postImage+".cal.map.eq.cub",
                           postImage+".cal.map.warp.cub",
                           postImage+".cal.map.warp_V2.cub",
                           preImage+".cal.map.cub",
                           postImage+".crop.cal.map.cub",
                           postImage+".warp.crop.cal.map.cub",
                           "duoCube.band0002.cub",
                           "duoCube.cub",
                           "duoCube.trim.cub",
                           "duoCube.crop.trim.cub",
                           "*_stats.log",
                           "autoRegTemp_*",
                           "imList*.lis"]
    
                print("Cleaning up...")
                for trash in removal:
                    proceed = sp.check_call(["rm -f "+trash], shell=True,
                                            executable="/bin/bash", cwd=pairDir,
                                            stdout=sys.stdout)
    
                finished = True
            else:
                finished = False
            cpFile = open(pairDir+"/deffile", "w+")
            cpReached = cpReached + 1
            cpFile.write(str(cpReached)+"\n")
            cpFile.close()
    
        cpFile = open(pairDir+"/deffile", "w+")
        if finished:
            cpFile.write("-1\n")
            cpFile.close()
            print("Processing finished")
            return
        else:
            cpFile.write(str(cpReached)+"/n")
            cpFile.close()
            print("Error processing pair")
            return

        
    def getData(self, preNames, postNames):
        preImages = []
        postImages = []
        for preImg in preNames:
            preImages.append(LunarImage(preImg, self.directory))
        for postImg in postNames:
            postImages.append(LunarImage(postImg, self.directory))
        return preImages, postImages
        
    def compare(self, inprelist, inpostlist):
        prelist = inprelist.copy()
        postlist = inpostlist.copy()
        pairs = []
        for postimg in postlist:
            for preimg in prelist:
                overlap = self.overlap(preimg, postimg)
                if (overlap.area > 0.001) and self.shadowsSame(preimg, postimg):
                    temporalPair = [preimg, postimg]    
                    pairs.append(temporalPair)
        
        print(str(len(pairs))+" temporal pairs found")
        return pairs
    
    def shadowsSame(self, preimg, postimg):
        sameDirection = ((preimg.getRelSolarLong() * postimg.getRelSolarLong()) > 0.)
        angleTolerance = (abs(preimg.getIncidenceAngle() - postimg.getIncidenceAngle()) <= self.app.delta_theta.get())
        angleClose = (abs(preimg.getIncidenceAngle() - postimg.getIncidenceAngle()) <= 5.)
        return ((angleTolerance or angleClose) and sameDirection)
    
    def overlap(self, preimg, postimg):
        
        prePoly = Polygon(MultiPoint([preimg.getUpperLeft(),
                                      preimg.getLowerLeft(), 
                                      preimg.getLowerRight(), 
                                      preimg.getUpperRight()]))
    
        postPoly = Polygon(MultiPoint([postimg.getUpperLeft(), 
                                       postimg.getLowerLeft(), 
                                       postimg.getLowerRight(), 
                                       postimg.getUpperRight()]))
    
        #return preimg.getBoundary().intersection(postimg.getBoundary())
        return prePoly.intersection(postPoly)
        
    def plotPair(self, preimg, postimg):
        
        for item in plt.gcf().axes:
            item.get_lines()[0].set_color("black")
        self.ax.fill(*self.overlap(preimg, postimg).exterior.xy, 
                 hatch='//', fill=False, color='red')
        self.fig.savefig(self.eventDir+"/location.png")
    
    def findImages(self):
        """
        Searches the Lunar Orbital Data Explorer web repository for images
        which meet the specified constraints.
    
        """
        print("Searching for Images...")
        
        LRO_START = "2009-01-01T00:00:00.000"
        LRO_END   = "2100-12-31T23:59:59.999"
        
        def longRange(long):
            while long < 0.:
                long += 360
            return long
        
        def genList(start, end):
            success = False
            while not success:
                try:
                    success = True
                    br = mech.Browser()
                    br.open("http://ode.rsl.wustl.edu/moon/productsearch.aspx")
                    br.select_form(name="Form1")
                    br.toggle_single(name="cb__LRO__LROC__EDRNAC__")
                    br["txtMaxLatitude"]         = str(self.app.latitude.get() + self.app.currentSearchSpace.get()/2.)
                    br["txtWestLongitude"]       = str(longRange(self.app.longitude.get() - self.app.currentSearchSpace.get()/2.))
                    br["txtEastLongitude"]       = str(longRange(self.app.longitude.get() + self.app.currentSearchSpace.get()/2.))
                    br["txtMinLatitude"]         = str(self.app.latitude.get() - self.app.currentSearchSpace.get()/2.)
                    br["txtStart"]               = start
                    br["txtEnd"]                 = end
                    br["txtIncidenceAngleStart"] = str(0)
                    br["txtIncidenceAngleEnd"]   = str(self.app.theta_i.get())
        
                    searchSubmit = br.find_control(name="btnSearchb")
                    searchSubmit.disabled = False
        
                    br.submit(searchSubmit.id)
                    br.open("http://ode.rsl.wustl.edu/moon/searchResultsList.aspx")
                    imgList = []
                    for link in br.links():
                        if (None != re.match("[M]\d+[RL][E]",link.text)):
                            imgList.append(str(link.text))
        
                # If the link isnt found try again. Wait 15s to avoid becoming an unintentional DDoS
                except mech.URLError:
                    print("Error connecting to server, trying again in 15s...")
                    time.sleep(15)
                    success = False
            
            success = False
            return imgList
        
        date = str(self.app.epoch)
        preList = genList(LRO_START, date)
        postList = genList(date, LRO_END)

        print(str(len(preList))+" before and "+str(len(postList))+" after images found")        
        return preList, postList



class LunarImage:
    
    """
    Lunar Image class is used to create a LunarImage object, which holds all
    the relevant data on the image location, incidence angle.
    """
    def __init__(self, imageName, workingDir):

        """
        In initialising the LunarImage object the relevant data is collected by
        calling the collectData method and propagating its output to the
        LunarImage object.
        """        
        self.name = str(imageName)
        self.directory = workingDir
        filepath = os.path.normpath(workingDir+"/img/IMAGE_DATA")
        if os.path.isfile(filepath):
            datafound = self.checkFile(filepath)
        else:
            datafound = False
        
        if not datafound:
            datafound = self.downloadData(filepath)    
            
            
    def checkFile(self, filepath):
        with open(filepath) as infile:
            for line in infile:
                if self.name in line:
                    dataList = line.split(" ")
                    self.upperLeftLatitude   = float(dataList[1])
                    self.upperLeftLongitude  = float(dataList[2])
                    self.lowerLeftLatitude   = float(dataList[3])
                    self.lowerLeftLongitude  = float(dataList[4])
                    self.lowerRightLatitude  = float(dataList[5])
                    self.lowerRightLongitude = float(dataList[6])
                    self.upperRightLatitude  = float(dataList[7])
                    self.upperRightLongitude = float(dataList[8])
                    self.emissionAngle       = float(dataList[9])
                    self.incidenceAngle      = float(dataList[10])
                    self.phaseAngle          = float(dataList[11])
                    self.solarLatitude       = float(dataList[12])
                    self.solarLongitude      = float(dataList[13])
                    return True
        return False
    
    def downloadData(self, filepath):
        
        timeout = [1,5,10,30,60]
        count = 0
        found = False
        while not found:
            try:
                found = self.collectData()
            except Exception as e:
                if "reset" in str(e):
                    print("Connection lost.... Retrying in "+str(timeout[count])+"seconds")
                    time.sleep(timeout[count])
                else:
                    print(e)
            else:
                count += 1
        
        with open(filepath, 'a') as infile:
            infile.write(self.name+" "+str(self.upperLeftLatitude)+" "+
                         str(self.upperLeftLongitude)+" "+str(self.lowerLeftLatitude)+" "+
                         str(self.lowerLeftLongitude)+" "+str(self.lowerRightLatitude)+" "+
                         str(self.lowerRightLongitude)+" "+str(self.upperRightLatitude)+" "+
                         str(self.upperRightLongitude)+" "+str(self.emissionAngle)+" "+
                         str(self.incidenceAngle)+" "+str(self.phaseAngle)+" "+
                         str(self.solarLatitude)+" "+str(self.solarLongitude)+"\n")
        return True

    def collectData(self):
        """
        Collects data on the image. If the data can be found in the pairDef.txt
        file, it uses that data. If the data is not present in pairDef.txt it
        connects to the WMS LROC web repository and dowloads the data, also
        propagating it to pairDef.txt for use in future runs.
        """
        success = False
        while not success:
            try:
                br = mech.Browser()
                br.set_handle_robots(False)
                formFound = False
                formAttempts = 0
                while not formFound:
                    try:
                        br.open("http://wms.lroc.asu.edu/lroc/search")
                        br.select_form(nr=0)
                        formFound = True
                    except mech.FormNotFoundError:
                        print("form not found")
                        time.sleep(3)
                        formFound = False
                        if formAttempts == 10:
                            return []

                br["filter[product_id]"] = self.name
                br.submit()
                try: br.find_link(text=self.name)
                except mech.LinkNotFoundError:
                    print("Image "+self.name+" not found")
                    return []
                dataPage = br.follow_link(text=self.name)
                soup = BeautifulSoup(dataPage, features="html5lib")
                soup.prettify()
                for row in soup.findAll('tr')[1:]:
                    col = row.findAll('td')
                    if ("Upper left latitude" in str(col[0])):
                        self.upperLeftLatitude = float(str(col[1])[4:-5])
                    if ("Upper left longitude" in str(col[0])):
                        self.upperLeftLongitude = float(str(col[1])[4:-5])
                    if ("Lower left latitude" in str(col[0])):
                        self.lowerLeftLatitude = float(str(col[1])[4:-5])
                    if ("Lower left longitude" in str(col[0])):
                        self.lowerLeftLongitude = float(str(col[1])[4:-5])
                    if ("Lower right latitude" in str(col[0])):
                        self.lowerRightLatitude = float(str(col[1])[4:-5])
                    if ("Lower right longitude" in str(col[0])):
                        self.lowerRightLongitude = float(str(col[1])[4:-5])
                    if ("Upper right latitude" in str(col[0])):
                        self.upperRightLatitude = float(str(col[1])[4:-5])
                    if ("Upper right longitude" in str(col[0])):
                        self.upperRightLongitude = float(str(col[1])[4:-5])
                    if ("Emission angle" in str(col[0])):
                        self.emissionAngle = float(str(col[1])[4:-5])
                    if ("Incidence angle" in str(col[0])):
                        self.incidenceAngle = float(str(col[1])[4:-5])
                    if ("Phase angle" in str(col[0])):
                        self.phaseAngle = float(str(col[1])[4:-5])
                    if ("Sub solar latitude" in str(col[0])[4:-5]):
                        self.solarLatitude = float(str(col[1])[4:-5])
                    if ("Sub solar longitude" in str(col[0])):
                        self.solarLongitude = float(str(col[1])[4:-5])
                        
                success = True

            except mech.URLError:
                success = False

        return success

    def longRange(self, long):
        if long > 180:
            return long - 360.
        else:
            return long

    def getUpperLeft(self):
        """ Returns the coordinates of the upper left most pixel in the image. """
        return self.longRange(self.upperLeftLongitude), self.upperLeftLatitude

    def getLowerLeft(self):
        """ Returns the coordinates of the lower left most pixel in the image. """
        return self.longRange(self.lowerLeftLongitude), self.lowerLeftLatitude

    def getLowerRight(self):
        """ Returns the coordinates of the lower right most pixel in the image. """
        return self.longRange(self.lowerRightLongitude), self.lowerRightLatitude

    def getUpperRight(self):
        """ Returns the coordinates of the upper right most pixel in the image. """
        return self.longRange(self.upperRightLongitude), self.upperRightLatitude

    def getIncidenceAngle(self):
        """ Returns the incidence angle at which the image was collected. """
        return self.incidenceAngle

    def getPhaseAngle(self):
        """ Returns the phase angle at which the image was collected. """
        return self.phaseAngle
       
    def getSolarLong(self):
        """ Returns the solar longitude at the time of image collection. """
        return self.solarLongitude

    def getRelSolarLong(self):
        """ Returns the displacement of the solar longitude from the normal. """
        return ((self.upperLeftLongitude + self.upperRightLongitude + self.lowerLeftLongitude + self.lowerRightLongitude) / 4.) - self.solarLongitude    

    def getBoundary(self):
        #add longrange
        return Polygon(MultiPoint([[self.upperLeftLongitude, self.upperLeftLatitude],
                                   [self.lowerLeftLongitude, self.lowerLeftLatitude],
                                   [self.lowerRightLongitude, self.lowerRightLatitude],
                                   [self.upperRightLongitude, self.upperRightLatitude]]))
    
    def downloadImage(self):
        """
        Downloads the specified image from the WMS LROC web repository into
        the specified directory.
        """
        
        if os.path.isfile(self.directory+"/img/"+self.name+".IMG"):
            return
        success = False
        while not success:
            try:
                success = True
                br = mech.Browser()
                br.set_handle_robots(False)
                br.open("http://wms.lroc.asu.edu/lroc/search")
                br.select_form(nr=0)
                br["filter[product_id]"] = self.name
                br.submit()
                br.find_link(text="Get download URLs for this page")
                downloadLink = br.follow_link(text="Get download URLs for this page")
                print("Downloading image "+self.name)
                linkToDownload = ("wget -nv "+downloadLink.read().decode("utf-8").split("\n")[0])

                sp.check_call(linkToDownload, shell=True, executable="/bin/bash", cwd=self.directory+"/img", stderr=sys.stdout, stdout=DEVNULL)

            except mech.URLError:
                success = False
            except mech.LinkNotFoundError:
                success = False
            except Exception as e:
                #if 'no space left on drive' in e:
                #    return "ERROR"
                #else:                    
                print(e)
                success = False
        
        calibrationOperations = ["lronac2isis from="+self.name+".IMG to="+self.name+".cub",
                                 "spiceinit from="+self.name+".cub cksmithed=yes spksmithed=yes shape=user model=/$ISISROOT/../data/base/dems/LRO_LOLA_LDEM_global_128ppd_20100915_0002.cub",
                                 "spicefit from="+self.name+".cub",
                                 "lronaccal from="+self.name+".cub to="+self.name+".cal.cub",
                                 "rm "+self.name+".IMG",
                                 "rm "+self.name+".cub"]
        
        print("Calibrating image "+self.name)
        for operation in calibrationOperations:
            sp.check_call([". ./../runIsis.sh && "+operation], shell=True, executable="/bin/bash", cwd=self.directory+"/img", stdout=DEVNULL, stderr=DEVNULL)
                        
        return self.name


class fakeGUI:
    def __init__(self):
        self.gatherInfo()
        self.directory = os.path.normpath(os.getcwd())
        self.splash_screen()
        #ADD INFO GATHERING 
    
    def gatherInfo(self):
        
        correctInfo = False
        while not correctInfo:
            latCheck, lonCheck, dateCheck =  False, False, False
    
            # Prompt user input of Latitude and loops until valid value is given
            while not latCheck:
                LATITUDE = input("Enter Latitude: ")
                try:
                    LATITUDE = float(LATITUDE)
                    latCheck = True
                except ValueError:
                    print("Latitude must be a number.")
                    print("")
    
            # Prompt user input of Longitude and loops until valid value is given
            while not lonCheck:
                LONGITUDE = input("Enter Longitude: ")
                try:
                    LONGITUDE = float(LONGITUDE)
                    lonCheck = True
                except ValueError:
                    print("Longitude must be a number.")
                    print("")
    
            # Prompy user input of date and loops until valid value is given
            while not dateCheck:
                DATE = input("Enter Date: ")
                try:
                    DATE = str(DATE) #MAKE IT CHECK ITS A VALID DATE
                    dateCheck = True
                except ValueError:
                    print("Please enter a valid date.")
                    print("")
    
    
            # Prompt confirmation of values before beginning algorithm
            print("Lat: "+str(LATITUDE)+" Long: "+str(LONGITUDE)+" Date: "+str(DATE))
            isAnswered = False
            while not isAnswered:
                proceedCheck = input("Run PyNAPLE with these parameters? (Y/N): ")
                acceptableAnswer = ["y", "Y", "n", "N"]
                if proceedCheck in acceptableAnswer:
                    isAnswered = True
    
            if (proceedCheck == "y") or (proceedCheck == "Y"):
                correctInfo = True
                self.latitude = LATITUDE
                self.longitude = LONGITUDE
                self.epoch = DATE
    
    def splash_screen(self):
        print("")
        print("    \||/                                                          \||/")
        print("    \||/      ________   ___   _   ___  ______ _      _____       \||/")
        print("  .<><><>.    | ___ \ \ / / \ | | / _ \ | ___ \ |    |  ___|    .<><><>.")
        print(" .<><><><>.   | |_/ /\ V /|  \| |/ /_\ \| |_/ / |    | |__     .<><><><>.")
        print(" '<><><><>'   |  __/  \ / | . ` ||  _  ||  __/| |    |  __|    '<><><><>'")
        print(" '<><><><>'   | |     | | | |\  || | | || |   | |____| |___    '<><><><>'")
        print(" '<><><><>'   \_|     \_/ \_| \_/\_| |_/\_|   \_____/\____/    '<><><><>'")
        print("  '<><><>'                                                      '<><><>'")
        print("")
        
    def handleArgs(self, arglist):
        
        for arg in arglist:
            if arg == '-h'.casefold():
                runHelp()
        
        
        
class StdRedirector(object):
    """ A stdout redirector class to redirect the stdout to the GUI text widget. """
    def __init__(self, text_widget):
        self.text_space = text_widget

    def fileno(self):
        """ Returns 1 to indicate that its the STDOUT. """
        return 1

    def write(self, string):
        """ Writes the stdout stream to the GUI text widget and updates it. """
        self.text_space.configure(state=tk.NORMAL)
        self.text_space.insert("end", string)
        self.text_space.update_idletasks()
        self.text_space.see("end")
        self.text_space.configure(state=tk.DISABLED)

    def flush(self):
        """ Updates idle tasks to act as a flush method. """
        self.text_space.configure(state=tk.NORMAL)
        self.text_space.update_idletasks()
        self.text_space.configure(state=tk.DISABLED)


class ErrRedirector(object):
    """ A stderr redirector class to redirect the stdout to the GUI text widget. """
    def __init__(self, text_widget):
        self.text_space = text_widget

    def fileno(self):
        """ Returns 2 to indicate that its the STDERR. """
        return 2

    def write(self, string):
        """ Writes the stderr stream to the GUI text widget and updates it. """
        self.text_space.configure(state=tk.NORMAL)
        self.text_space.insert("end", string)
        self.text_space.update_idletasks()
        self.text_space.see("end")
        self.text_space.configure(state=tk.DISABLED)

    def flush(self):
        """ Updates idle tasks to act as a flush method. """
        self.text_space.configure(state=tk.NORMAL)
        self.text_space.update_idletasks()
        self.text_space.configure(state=tk.DISABLED)


class PyNAPLEGUI(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.title("PyNAPLE")
        #self.master.geometry("600x800")
        self.configRead()
        self.directory = tk.StringVar(value=os.path.normpath(os.getcwd()))
        self.latitude = tk.DoubleVar(value=0.0)
        self.longitude = tk.DoubleVar(value=0.0)
        self.checked = tk.BooleanVar(value=True)
        self.searchSpace = tk.DoubleVar(value=self.configlist[0])
        self.searchIncrement = tk.DoubleVar(value=self.configlist[1])
        self.theta_i = tk.DoubleVar(value=self.configlist[2])
        self.delta_theta = tk.DoubleVar(value=self.configlist[3])
        self.algorithm = tk.StringVar(value=self.configlist[4])
        self.tolerance = tk.DoubleVar(value=self.configlist[5])
        self.currentSearchSpace = tk.DoubleVar(value=0.0)
        self.pairsFound = tk.IntVar(value=0)
        self.pairsProcessed = tk.IntVar(value=0)
        self.currentPreImage = tk.StringVar(value='')
        self.currentPostImage = tk.StringVar(value='')
        self.create_widgets()


    def configRead(self):
        with open("CONFIG") as config:
            self.configlist = [a.split(" = ")[1] for a in config.read().splitlines()]
            

    def create_widgets(self):
       
        self.topFrame = tk.Frame(self, width=900)
        self.topFrame.grid(row=0, column=0, padx=10, pady=10)
        
        self.frame = tk.Frame(self.topFrame, width=550, height=200)
        self.frame.grid(row=0, column=0, sticky='ew')

        self.folderLab = tk.Entry(self.frame, textvariable=self.directory)
        self.folderLab.grid(row=0, column=0, columnspan=4, sticky='ew')

        self.folderBttn = Button(self.frame, command=self.reqDir, text='Choose Dir')
        self.folderBttn.grid(row=0, column=5, sticky='ew')
       
        self.latLab = tk.Label(self.frame, text='Latitude, \u03C6 (-90 to 90):')
        self.latLab.grid(row=1, column=0, columnspan=2, sticky='w')
        
        self.latEntry = tk.Entry(self.frame, textvariable=self.latitude, width=10)
        self.latEntry.grid(row=1, column=2, columnspan=1, sticky='w')
        
        self.lonLab = tk.Label(self.frame, text='Longitude, \u03BB (-180 to 180):')
        self.lonLab.grid(row=2, column=0, columnspan=2, sticky='w')
        
        self.lonEntry = tk.Entry(self.frame, textvariable=self.longitude, width=10)
        self.lonEntry.grid(row=2, column=2, columnspan=1, sticky='w')
        
        self.epochLab = tk.Label(self.frame, text='Flash epoch:')
        self.epochLab.grid(row=1, column=3, columnspan=1)
   
        self.epochEntry = DateEntry(self.frame, date_pattern='y-mm-dd')
        self.epochEntry.grid(row=1, column=5)
        
        self.defaultCheck = tk.Checkbutton(self.frame, text='Use default values', variable=self.checked, onvalue=True, command=self.disenframe)
        self.defaultCheck.grid(row=2, column=3, columnspan=3)
        
#        self.vline = tk.ttk.Separator(self.frame, orient='vertical')
#        self.vline.grid(row=0, column=6, rowspan=3, sticky='ns', padx=10)
        
        ########################
        
        self.line1 = tk.ttk.Separator(self.topFrame, orient='horizontal')
        self.line1.grid(row=1, column=0, columnspan=6, sticky='news', pady=10)
        
        ########################
        
        self.valsFrame = tk.Frame(self.topFrame, width=550, height=200)
        self.valsFrame.grid(row=0, column=1, sticky='ew')
        
        self.sspLab = tk.Label(self.valsFrame, text='Search Space (\u00B0)', state='disabled')
        self.sspLab.grid(row=0, column=0, sticky='w')
        
        self.sspEntry = tk.Entry(self.valsFrame, textvariable=self.searchSpace, width=10, state='disabled')
        self.sspEntry.grid(row=0, column=1)
        
        self.ssiLab = tk.Label(self.valsFrame, text='Increment (\u00B0)', state='disabled')
        self.ssiLab.grid(row=0, column=2, sticky='w')        

        self.ssiEntry = tk.Entry(self.valsFrame, textvariable=self.searchIncrement, width=10, state='disabled')
        self.ssiEntry.grid(row=0, column=3)

        self.thetaLab = tk.Label(self.valsFrame, text='Maximum \u03B8 (\u00B0)', state='disabled')
        self.thetaLab.grid(row=1, column=0, sticky='w')
            
        self.thetaEntry = tk.Entry(self.valsFrame, textvariable=self.theta_i, width=10, state='disabled')
        self.thetaEntry.grid(row=1, column=1)
        
        self.deltaThetaLab = tk.Label(self.valsFrame, text='Maximum \u0394\u03B8 (\u00B0)', state='disabled')
        self.deltaThetaLab.grid(row=1, column=2, sticky='w')
        
        self.deltaThetaEntry = tk.Entry(self.valsFrame, textvariable=self.delta_theta, width=10, state='disabled')
        self.deltaThetaEntry.grid(row=1, column=3)
        
        self.algorithmLab = tk.Label(self.valsFrame, text='Algorithm' , state='disabled')
        self.algorithmLab.grid(row=2, column=0, sticky='w')
        
        options = ['MaximumCorrelation', 'MinimumDifference', 'AdaptiveGruen']
        
        self.algorithmList = tk.OptionMenu(self.valsFrame, self.algorithm, *options)
        self.algorithmList.grid(row=2, column=0, columnspan=2, sticky='e')
        self.algorithmList.configure(state='disabled', width=14, font=("Helvetica", 8), anchor='w')
        
        self.toleranceLab = tk.Label(self.valsFrame, text='Tolerance (0-1)', state='disabled')
        self.toleranceLab.grid(row=2, column=2, sticky='w')
        
        self.toleranceEntry = tk.Entry(self.valsFrame, textvariable=self.tolerance, width=10, state='disabled')
        self.toleranceEntry.grid(row=2, column=3)
        
        ########################
        
        self.runFrame = tk.Frame(self.topFrame, width=550, height=100)
        self.runFrame.grid(row=3, column=0, columnspan=2, sticky='ew')
        
        self.runButton = Button(self.runFrame, text=u'\u23F5', font=("Helvetica", 15), fg='green', activeforeground='green', command=self.checkRunPynaple)
        self.runButton.grid(row=0, column=0, columnspan=2, sticky='ew')
        
        self.stopAftButton = Button(self.runFrame, text=u'\u23EF', font=("Helvetica", 15), fg='yellow', activeforeground='yellow', command=self.stopNext, state='disabled')
        self.stopAftButton.grid(row=0, column=2, sticky='ew')
        
        self.stopButton = Button(self.runFrame, text=u'\u23F9', font=("Helvetica", 15), fg='red', activeforeground='red', command=self.stop, state='disabled')
        self.stopButton.grid(row=0, column=3, sticky='ew')
        
        ########################
        
        self.line2 = tk.ttk.Separator(self.topFrame, orient='horizontal')
        self.line2.grid(row=4, column=0, columnspan=6, sticky='news', pady=10)
        
        ########################
        
        self.outputFrame = tk.Frame(self.topFrame)
        self.outputFrame.grid(row=5, column=0, columnspan=2, sticky='ew')
        
        self.outputWindow = st.ScrolledText(self.outputFrame, fg='Black', bg='White', width=90, height=25, state='disabled')
        self.outputWindow.pack(side='left')
        sys.stdout = StdRedirector(self.outputWindow)
        sys.stderr = ErrRedirector(self.outputWindow)
        
        self.monitorFrame = tk.Frame(self.outputFrame, width=125)
        self.monitorFrame.pack(side='right', fill='both')
        
        self.currentSSLab = tk.Label(self.monitorFrame, text='Current Search Area')
        self.currentSSLab.grid(row=0, column=0, sticky='w', columnspan=2)
        
        self.currentSSEntry = tk.Entry(self.monitorFrame, textvariable=self.currentSearchSpace, width=5, state='disabled')
        self.currentSSEntry.grid(row=1, column=1, sticky='e', columnspan=2)
        
        self.pairsFoundLab = tk.Label(self.monitorFrame, text='Pairs Found')
        self.pairsFoundLab.grid(row=2, column=0, sticky='w', columnspan=2)
        
        self.pairsFoundEntry = tk.Entry(self.monitorFrame, textvariable=self.pairsFound, width=5, state='disabled')
        self.pairsFoundEntry.grid(row=3, column=1, sticky='e', columnspan=2)
        
        self.pairsProcLab = tk.Label(self.monitorFrame, text='Pairs Processed')
        self.pairsProcLab.grid(row=4, column=0, sticky='w', columnspan=2)
        
        self.pairsProcEntry = tk.Entry(self.monitorFrame, textvariable=self.pairsProcessed, width=5, state='disabled')
        self.pairsProcEntry.grid(row=5, column=1, sticky='e', columnspan=2)
        
        self.currentPairLab = tk.Label(self.monitorFrame, text='Current Pair')
        self.currentPairLab.grid(row=6, column=0, sticky='w', columnspan=2)
        
        self.currentPairEntryPre = tk.Entry(self.monitorFrame, textvariable=self.currentPreImage, width=15, state='disabled')
        self.currentPairEntryPre.grid(row=7, column=1, sticky='e', columnspan=2)
        
        self.currentPairEntryPost = tk.Entry(self.monitorFrame, textvariable=self.currentPostImage, width=15, state='disabled')
        self.currentPairEntryPost.grid(row=8, column=1, sticky='e', columnspan=2)
        
        self.pack()

    
    def verifyInput(self):
        #todo check inputs
        return True
   
    def checkRunPynaple(self):
        
        if self.verifyInput():
            self.epoch = self.epochEntry.get_date()
            okText = "Run PyNAPLE on \n\u03C6 = "+str(self.latitude.get())+", \u03BB = "+str(self.longitude.get())+"\nat "+str(self.epoch)+"?"
            cont = tk.messagebox.askokcancel("Run PyNAPLE?", okText)
            
            if cont:
                self.runPyNAPLE()
            else:
                return
            
        else:
            return
        
    def runPyNAPLE(self):
        self.lockInput()
        self.pynapleInstance = PyNAPLE(self)
#        self.pynapleInstance = mp.Process(target=PyNAPLE, args=[self])
#        self.pynapleInstance.start()
    
    def lockInput(self):
        for child in self.frame.winfo_children():
            child.config(state='disabled')
        for child in self.valsFrame.winfo_children():
            child.config(state='disabled')
        self.stopAftButton.config(state='normal')
        self.stopButton.config(state='normal')

    def releaseInput(self):
        for child in self.frame.winfo_children():
            child.config(state='normal')
        for child in self.valsFrame.winfo_children():
            child.config(state='disabled' if self.checked.get() else 'normal')
        self.stopAftButton.config(state='disabled')
        self.stopButton.config(state='disabled')

    def stopNext(self):
        #TODO add guts
        self.releaseInput()
    
    def stop(self):
        #TODO add guts
        self.releaseInput()
    
    def disenframe(self):
        self.searchSpace.set(value=self.configlist[0])
        self.searchIncrement.set(value=self.configlist[1])
        self.theta_i.set(value=self.configlist[2])
        self.delta_theta.set(value=self.configlist[3])
        self.algorithm.set(value=self.configlist[4])
        self.tolerance.set(value=self.configlist[5])
        
        for child in self.valsFrame.winfo_children():
            child.config(state='disabled' if self.checked.get() else 'normal')
        

    def reqDir(self):
        direct = fd.askdirectory()
        if os.path.isdir(str(direct)):
            path = os.path.normpath(direct)
            self.directory.set(path)
        else:
            self.directory.set('')

        #path, file = os.path.split(filepath)
        #img = np.array(fits.open(filepath)[0].data)
        #self.frameSize = img.shape
        #self.pttnVar.set(patternise(file))
        self.folderLab.xview_moveto(1)
        
def runHelp():
        print("")
        print("    \||/                                                          \||/")
        print("    \||/      ________   ___   _   ___  ______ _      _____       \||/")
        print("  .<><><>.    | ___ \ \ / / \ | | / _ \ | ___ \ |    |  ___|    .<><><>.")
        print(" .<><><><>.   | |_/ /\ V /|  \| |/ /_\ \| |_/ / |    | |__     .<><><><>.")
        print(" '<><><><>'   |  __/  \ / | . ` ||  _  ||  __/| |    |  __|    '<><><><>'")
        print(" '<><><><>'   | |     | | | |\  || | | || |   | |____| |___    '<><><><>'")
        print(" '<><><><>'   \_|     \_/ \_| \_/\_| |_/\_|   \_____/\____/    '<><><><>'")
        print("  '<><><>'                                                      '<><><>'")
        print("")
        print("")
        print("NAME")
        print("        PyNAPLE - Python NAC Automated Temporal Pair Lunar Evaluator")
        print("")
        print("SYNOPSIS")
        print("        python PyNAPLE_<version> [-h] [-c] [-t] [-p img1 img2]")
        print("")
        print("DESCRIPTION")
        print("        PyNAPLE is an automated image processing pipeline for the det-")
        print("        ection of craters formed during lunar impact flash events. By ")
        print("        suppling PyNAPLE with the location and date of a lunar impact ")
        print("        flash event, the software automatically retrieves images which")
        print("        may contain the resulatant crater, and processes them in order")
        print("        to locate the crater.")
        print("")
        print("OPTIONS")
        print("        -h    Prints the help page, detailing the programs usage.")
        print("")
#        print("        -t    Starts the program in terminal mode")
#        print("")
#        print("        -c    Starts the program in terminal mode with multiprocessing")
#        print("              enabled")
#        print("")
#        print("        -p img1 img2")
#        print("              Processed img1 and img2 as a temporal pair. No checks a-")
#        print("              re made as to the images suitablility.")
#        print("")
        print("")

def pynapleHelper(guiObj):
    guiObj.instance = PyNAPLE(guiObj)
    
def imageDownloadHelper(LunarImg):
    """ Helper class to allow multiprocessing to spawn process. """
    LunarImg.downloadImage()
        
def patternise(file):
    filename, ext = os.path.splitext(file)
    return re.sub('[0-9][0-9][0-9][0-9]', 'XXXX', filename[::-1])[::-1] + ext

def on_closing(root):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        root.destroy()

main()

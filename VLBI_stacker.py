import os, re, time, datetime, sys, math, fnmatch
from os.path import join, getsize, isfile
from os import listdir
from datetime import date
from collections import deque
import Utilities
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
ti = time.time()
import numpy as np
import csv
import string
import StringIO

def findmaxb(uvdata):
	maxbaseline = 0
	antab = uvdata.table('AN',1)
	for row in antab :
		for row2 in antab :
			xsep = row.stabxyz[0] - row2.stabxyz[0]
			ysep = row.stabxyz[1] - row2.stabxyz[1]
			zsep = row.stabxyz[2] - row2.stabxyz[2]
			hypxy =  math.sqrt((xsep * xsep) + (ysep * ysep))
			hypxyz = math.sqrt((zsep * zsep) + (hypxy * hypxy))
			if hypxyz > maxbaseline :
				maxbaseline = hypxyz
	cellsize = (1.22 * (300000000 / uvdata.header.crval[2]) / maxbaseline) / 3.141592 * 180 * 3600 / 5
	print "maxbaseline = ", maxbaseline, "cellsize = ", cellsize
	return cellsize,cellsize

AIPS.userno = 1001
filestart = 'HDFC'
point_cen_RA = -170.791333333   #RA pointing centre of observations in decimal degrees
point_cen_DEC = 62.216111111
path = '~/nas/radcliff/EG078B/multi_source_s_cal/eMERGE/UV/'
#MANUALLY LOAD HDFS0001 first and rename to STACK.UV.1 (SORT this out when have more time)

with open('eMERGE_matched_VLBI_pos.csv','rb') as csvfile:
	spamreader = csv.reader(csvfile,delimiter=',')
	x = np.asarray(list(spamreader))
print x[0]
for file in os.listdir('./'):
	for i in range(len(x)):
		if file[:8] == x[i][5].upper():
			print file
			if AIPSUVData('STACK','UVDATA',1,1).exists() == False:
				fitld = AIPSTask('FITLD')
				fitld.datain = ('PWD:' + file)
				fitld.ncount = 1
				fitld.doconcat = 1
				fitld.clint = 0
				fitld.wtthresh = 0
				fitld.outname = 'STACK'
				fitld.outclass = 'UVDATA'
				fitld.outdisk = 1
				fitld.digicor = -1
				fitld.go()

				uvdata2 = AIPSUVData('STACK','UVDATA',1,1)
				
				uvfix = AIPSTask('UVFIX')
				uvfix.indata = uvdata2
				uvfix.outdisk = 1
				print 'Shifting by: %.3f %.3f' % (float(x[i][11]),float(x[i][12]))
				uvfix.shift[1:] = float(x[i][11]), float(x[i][12])
                                uvfix.fqcenter = -1
				uvfix.go()
				uvdata2.zap()
				uvdata2 = WizAIPSUVData('STACK','UVFIX',1,1)

				uvdata2.header['crval'][4] = point_cen_RA
				uvdata2.header.update()
				uvdata2.header['crval'][5] = point_cen_DEC
				uvdata2.header.update()
				uvdata2.rename(name='STACK',klass='UVDATA', seq=0)
				uvdata2 = AIPSUVData('STACK','UVDATA',1,1)

			else:
				fitld = AIPSTask('FITLD')
				fitld.datain = ('PWD:' + file)
				fitld.ncount = 1
				fitld.doconcat = 1
				fitld.clint = 0
				fitld.wtthresh = 0
				fitld.outname = 'STACK'
				fitld.outclass = 'UV'
				fitld.outdisk = 1
				fitld.digicor = -1
				fitld.go()

				uvdata = WizAIPSUVData('STACK','UV',1,1)
				uvfix = AIPSTask('UVFIX')
				uvfix.indata = uvdata
				uvfix.outdisk = 1
				print 'Shifting by: %.3f %.3f' % (float(x[i][11]),float(x[i][12]))
				uvfix.shift[1:] = float(x[i][11]), float(x[i][12])
				uvfix.fqcenter = -1
				uvfix.go()

				uvdata.zap()
				uvdata = WizAIPSUVData('STACK','UVFIX',1,1)

				uvdata.header['crval'][4] = point_cen_RA
				uvdata.header.update()
				uvdata.header['crval'][5] = point_cen_DEC
				uvdata.header.update()
				
				uvdata=AIPSUVData('STACK','UVFIX',1,1)
				dbapp = AIPSTask('DBAPP')
				dbapp.indata = uvdata
				dbapp.indisk = 1
				dbapp.inseq = 1
				dbapp.outdata = uvdata2
				dbapp.outdisk = 1
				dbapp.go()
				uvdata.zap()
				

uvsrt = AIPSTask('UVSRT')                                                      
uvsrt.sort = 'TB'                                                              
uvsrt.indata = uvdata2                                                         
                                                        
uvsrt.outdisk = 1                                                              
uvsrt.go()   
uvdata2.zap()
uvdata2 = WizAIPSUVData('STACK','UVSRT',1,1)                                                                
'''
indxr = AIPSTask('INDXR')                                                      
indxr.indata = uvdata2                                                         
indxr.cparm[1]=360                                                             
indxr.cparm[2] = 360                                                           
indxr.cparm[3] = 0.25                                                          
indxr.go()                                                                     
'''
nchan = 32
imsize = [512,512]
niter = 1
imagr = AIPSTask('IMAGR')                                                      
imagr.nchav = nchan #use imagr to get a clean model!                           
imagr.indata = uvdata2                                                         
imagr.sources[1] = ''                                                          
imagr.outname = 'STACK_IM'                                                     
imagr.cellsize[1:] = findmaxb(uvdata2)                                          
imagr.imsize[1:] = imsize                                                      
imagr.nboxes = 1                                                               
imagr.nfield = 1                                                               
imagr.outdisk = 1                                                              
imagr.uvwtfn = 'NA'                                                           
imagr.niter = niter                                                            
imagr.go() 
         
'''
for i in range(1):
	j=i+1

	if j < 10:
		uvdataname = filestart+'000'+str(j)+'_split.fits' #Load the UV data file & image
	if j < 100 and j>9:
		uvdataname = filestart+'00'+str(j)+'_split.fits'
	if j < 1000 and j>99:
		uvdataname = filestart+'0'+str(j)+'_split.fits'
	name = str(x[i][0])

	fitld = AIPSTask('FITLD')
	fitld.datain = ('PWD:' + uvdataname)
	fitld.ncount = 1
	fitld.doconcat = 1
	fitld.clint = 0
	fitld.wtthresh = 0
	fitld.outname = str(x[i][0])
	fitld.outclass = 'UV'
	fitld.outdisk = 1
	fitld.digicor = -1
	fitld.go()

	y1=x[i][1]
	y2=x[i][2]
	rasec1 = 15*float(y1[0:2])+(float(y1[3:5])/4.)+(float(y1[6:11])/240.)
	print y1[0:2], y1[3:5], y1[6:11]
	dasec1 = float(y2[0:2])+(float(y2[3:5])/60.)+(float(y2[6:11])/3600.)
	print y2[0:2],y2[3:5],y2[6:11]
	print (float(rasec1)-360.), (float(x[i][3])-360.), dasec1, x[i][4]

 	dec1rad=dasec1*math.pi/180.
	diffra = (((float(x[i][3])-rasec1)*np.cos(dec1rad))*3600.)
	diffdec = ((float(x[i][4])-dasec1)*3600.)

	uvdata = AIPSUVData(name,'UV',1,1)
	uvfix = AIPSTask('UVFIX')
	uvfix.indata = uvdata
	uvfix.outdisk = 1
	uvfix.shift[1:] = float(diffra), float(diffdec)
	uvfix.go()

	uvdata.zap()
	uvdata2 = WizAIPSUVData(name,'UVFIX',1,1)

	uvdata2.header['crval'][4] = point_cen_RA
	uvdata2.header.update()
	uvdata2.header['crval'][5] = point_cen_DEC
	uvdata2.header.update()
	uvdata2.rename(name='STACK',klass='UV', seq=0)
	uvdata2 = AIPSUVData('STACK','UV',1,1)

for i in range(1,151):
	j=i+1

	if j < 10:
		uvdataname = 'HDFS000'+str(j)+'_split.fits' #Load the UV data file & image
	if j < 100 and j>9:
		uvdataname = 'HDFS00'+str(j)+'_split.fits'
	if j < 1000 and j>99:
		uvdataname = 'HDFS0'+str(j)+'_split.fits'
	name = str(x[i][0])

	fitld = AIPSTask('FITLD')
	fitld.datain = ('PWD:' + uvdataname)
	fitld.ncount = 1
	fitld.doconcat = 1
	fitld.clint = 0
	fitld.wtthresh = 0
	fitld.outname = str(x[i][0])
	fitld.outclass = 'UV'
	fitld.outdisk = 1
	fitld.digicor = -1
	fitld.go()

	y1=x[i][1]
	y2=x[i][2]
	rasec1 = 15*float(y1[0:2])+(float(y1[3:5])/4.)+(float(y1[6:11])/240.)
	print y1[0:2], y1[3:5], y1[6:11]
	dasec1 = float(y2[0:2])+(float(y2[3:5])/60.)+(float(y2[6:11])/3600.)
	print y2[0:2],y2[3:5],y2[6:11]
	print (float(rasec1)-360.), (float(x[i][3])-360.), dasec1, x[i][4]

 	dec1rad=dasec1*math.pi/180.
	diffra = (((float(x[i][3])-rasec1)*np.cos(dec1rad))*3600.)
	diffdec = ((float(x[i][4])-dasec1)*3600.)

	uvdata = AIPSUVData(name,'UV',1,1)
	uvfix = AIPSTask('UVFIX')
	uvfix.indata = uvdata
	uvfix.outdisk = 1
	uvfix.shift[1:] = float(diffra), float(diffdec)
	uvfix.go()

	uvdata.zap()
	uvdata = WizAIPSUVData(name,'UVFIX',1,1)

	uvdata.header['crval'][4] = point_cen_RA
	uvdata.header.update()
	uvdata.header['crval'][5] = point_cen_DEC
	uvdata.header.update()

	uvdata=AIPSUVData(name,'UVFIX',1,1)
	dbapp = AIPSTask('DBAPP')
	dbapp.indata = uvdata
	dbapp.indisk = 1
	dbapp.inseq = 1
	dbapp.outdata = uvdata2
	dbapp.outdisk = 1
	dbapp.go()
	uvdata.zap()

uvsrt = AIPSTask('UVSRT')
uvsrt.sort = 'TB'
uvsrt.indata = uvdata2
uvsrt.outdata = uvdata2
uvsrt.outdisk = 1
uvsrt.go()

indxr = AIPSTask('INDXR')
indxr.indata = uvdata2
indxr.cparm[1]=360
indxr.cparm[2] = 360
indxr.cparm[3] = 0.25
indxr.go()

nchan = 32
imsize = [512,512]
niter=500

imagr = AIPSTask('IMAGR')
imagr.nchav = nchan #use imagr to get a clean model!
imagr.indata = uvdata2
imagr.sources[1] = ''
imagr.outname = 'STACK_IM'
imagr.cellsize[1:] = findmaxb(uvdata)
imagr.imsize[1:] = imsize
imagr.nboxes = 1
imagr.nfield = 1
imagr.outdisk = 1
imagr.uvwtfn = ''
imagr.niter = niter
imagr.go()

'''

'''
uvdata2 = WizAIPSUVData(name, 'UVSUB',2,1)
uvdata2.header['crval'][4] = point_cen_RA
uvdata2.header.update()
uvdata2.header['crval'][5] = point_cen_DEC
uvdata2.header.update()
'''

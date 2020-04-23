#!BPY

"""
Name: 'Corona Render (.nra2) ...'
Blender: 245
Group: 'Export'
Tooltip: 'Exports and renders nra2 format'
"""

import math
import sys
import struct, os
import Blender
from Blender import *
from Blender.BGL import *
from Blender.Mathutils import *
import bpy
import subprocess
from threading import Thread

# try:
# 	import psyco
# 	psyco.full()
# 	# print "PSYCO found - !!"
# except:
# 	# print "PSYCO not found - no compiling "
# 	pass


global materials

materials = []


################ Set path for corona executable and tmpdir #################

coronadir = "DREGGN_INSTALLDIR"
datadir = "DREGGN_DATADIR"
tmpdir = "DREGGN_TMPDIR"
blenddir = "DREGGN_BLENDDIR"
corona_export_version = "DREGGN_VERSION"

############################################################################

class renderer(Thread):
	def __init__ (self,cmdline):
		Thread.__init__(self)
		self.cmdline = cmdline
	def run(self):
		print "calling "+self.cmdline
		fp = os.popen(self.cmdline,"r")
		while 1:
			line = fp.readline()
			if not line: break

def getTextureFileName(filename):
	f = filename[filename.rfind('/')+1:-4]
	f = f[f.rfind('\\')+1:]
	f = f.replace(' ', '_')
	# as textures don't change over time, this is not necessary:
	# f += "%04d"%Blender.Get("curframe")
	return f

def writetex(image,filename):
	
	try:
		size = image.getSize()
	except:
		return
	
	ifilename = image.filename
	xsize = size[0]
	ysize = size[1]
	
	imagefilename = sys.dirname(filename) + "/" + getTextureFileName(ifilename)+".tga"

	bpp = image.getDepth()
	if bpp != 24 and bpp != 32:
		print "ERROR: not 24 or 32 bits per pixel!"
		return

	if not sys.exists(imagefilename):
		print "Writing image "+str(imagefilename)
		simage = open(imagefilename, 'wb')
		simage.write(struct.pack("B", 0))
		simage.write(struct.pack("B", 0))
		simage.write(struct.pack("B", 2))
		simage.write(struct.pack("h", 0))
		simage.write(struct.pack("h", 0))
		simage.write(struct.pack("B", 0))
		simage.write(struct.pack("h", 0))
		simage.write(struct.pack("h", 0))
		simage.write(struct.pack("h", xsize))
		simage.write(struct.pack("h", ysize))
		simage.write(struct.pack("B", bpp))
		simage.write(struct.pack("B", 0))

		for y in range(0,ysize):
			for x in range(0, xsize):
				pix = image.getPixelI(x, ysize-1-y)
				simage.write("%c%c%c" % (pix[2], pix[1], pix[0]))
				if bpp == 32: simage.write("%c" % pix[3])
		simage.close()
		
def writenormalmap(image, filename, intensity):
	
	try:
		size = image.getSize()
	except:
		return
	
	ifilename = image.filename
	ppmname = sys.dirname(filename) + "/" + getTextureFileName(ifilename)+".ppm"
	
	pgmname = sys.dirname(filename) + "/" + getTextureFileName(ifilename)+".pgm"

	bpp = image.getDepth()
	
	ftype = 'P5'
	xsize = size[0]
	ysize = size[1]
	
	if not sys.exists(ppmname):
		
		simage = open(pgmname, 'wb')
		
		simage.write("%s\n" % ftype)
		simage.write("%d %d\n" % (xsize, ysize))
		simage.write("255\n")
			
		for y in range(0,ysize):
			for x in range(0, xsize):
				pix = image.getPixelI(x, ysize-1-y)
				if bpp == 8:
					simage.write("%c" % pix[0])
				else:
					gray = int((pix[2] + pix[1] + pix[0])/3)
					simage.write("%c" % gray )

		simage.close()
			
		bump2normmap = datadir+"tools/bump2normalmap"
		
		cmd = ("%s %s %s %s\n") % (bump2normmap, intensity, pgmname, ppmname)
		
		os.system(cmd)
		
		delpgm = ("rm %s" % pgmname)
		
		os.system(delpgm)
	
##################################################

def writemat(mat, mat_str, mat_num, num_mats, filename):
	FLAGS = mat.getMode()
	
	brdf_mat = num_mats
	glass = 1
	
	ior1 = round(mat.IOR,5 )
	spect = 1-mat.specTransp
	ior2 = ior1-round(spect, 5)
	
	# if mat.enableSSS:
	# 	if FLAGS & Material.Modes['RAYTRANSP']:
	# 		mat_str += 'solidglass %s %s # %s\n' % (ior1, ior2, num_mats)
	# 	else:
	# 		hard = str(mat.getHardness()*20)
	# 		mat_str += 'milkglass %f # %s\n' % (hard, num_mats)
	# 	
	if FLAGS & Material.Modes['RAYTRANSP']:
		if mat.glossTra < 1.0:
			hard = mat.getHardness()*20
			mat_str += 'milkglass %f # %s\n' % (hard, num_mats)

		elif ior1 == 1.0:
			mat_str += 'glass %s # %s\n' % (spect, num_mats)
	
		else:
			mat_str += 'solidglass %s %s # %s\n' % (ior1, ior2, num_mats)
			glass = 0
	
	elif FLAGS & Material.Modes['RAYMIRROR']:
		mat_str += 'glass 1.0 # %s\n' % (num_mats)
	
	elif mat.spec != 0:
		hard = str(mat.getHardness()*20)
		mat_str += 'ashi %s %s # %s\n' % (hard, hard, num_mats)
	
	else:
		mat_str += 'diffuse # %s\n' %  (num_mats)
		
	num_mats += 1
	textures = mat.getTextures()
	diffFileName = 'none'
	specFileName = 'none'
	emitFileName = 'none'
	bumpFileName = 'none'
	bumpStrength = 0

	for tex in textures:
		if tex == None: continue
		if tex.tex.type != Texture.Types.IMAGE: continue
		if tex.mtCol != 0:
			diffFileName = getTextureFileName(tex.tex.image.filename)+".tga"
			writetex(tex.tex.image, filename)
		if tex.mtSpec != 0:
			specFileName = getTextureFileName(tex.tex.image.filename)+".tga"
			writetex(tex.tex.image, filename)
		if tex.mtNor != 0:
			bumpFileName = getTextureFileName(tex.tex.image.filename)+".ppm"
			bumpStrength = tex.norfac / 5.0
			writenormalmap(tex.tex.image, filename, bumpStrength)
		if tex.mtEmit!= 0:
			emitFileName = getTextureFileName(tex.tex.image.filename)+".tga"
			writetex(tex.tex.image, filename)
	
	c = mat.rgbCol
	ref = mat.getRef()
	
	if mat.name != "hair":
		mat_str += 'normals  %s # %s\n' % (brdf_mat, num_mats)
		num_mats += 1
		
	## emitter Texture ##
	if mat.emit != 0 :
		if emitFileName == 'none':
			mat_str += 'color e %s %s %s %s # %s\n' % (round(c[0]*mat.emit*10000,5), round(c[1]*mat.emit*10000,5), round(c[2]*mat.emit*10000,5), brdf_mat, num_mats)
		else:
			mat_str += 'texture e %s %s # %s\n' % (brdf_mat, emitFileName, num_mats)			
		num_mats += 1
	
	## diffuse Texture ##
	if diffFileName == 'none' :
		if FLAGS & Material.Modes['RAYTRANSP'] and not FLAGS & Material.Modes['RAYMIRROR']:
			# this is legacy stuff. absorption is now handled in medium_*
			# if not glass:
			# 	mat_str += 'color d %s %s %s %s # %s\n' % (round(1/(c[0]+0.001)-1,5), round(1/(c[1]+0.001)-1,5), round(1/(c[2]+0.001)-1,5), brdf_mat, num_mats)
			# else:
			mat_str += 'color d %s %s %s %s # %s\n' % (round(c[0],5), round(c[1],5), round(c[2],5),  brdf_mat, num_mats)
		else:
			mat_str += 'color d %s %s %s %s # %s\n' % (round(c[0]*ref,5), round(c[1]*ref,5), round(c[2]*ref,5), brdf_mat, num_mats)
		num_mats += 1
	
	else:
		mat_str += 'texture d %s %s # %s\n' % (brdf_mat, diffFileName, num_mats)
		num_mats += 1
	
	## SSS ###
	
	if mat.enableSSS:
		s=mat.sssCol
		mat_str += 'color v %s %s %s %s # %s\n' % (round(s[0],5), round(s[1],5), round(s[2],5), brdf_mat, num_mats)
		num_mats += 1
		scale = mat.sssScale*10 # account for 1.f => 1m to 1.f => 1dm translation
		# mat_str += 'medium_iso %s %s # %s\n' % (round(scale,5), brdf_mat, num_mats)
		ar = mat.sssRadiusRed
		ag = mat.sssRadiusGreen
		ab = mat.sssRadiusBlue
		mat_str += 'medium_iso_rgb %s %s %s %s # %s\n' % (round(scale*ar,5), round(scale*ag,5), round(scale*ab,5), brdf_mat, num_mats)
		num_mats += 1
		
	## specular Texture ##
	if mat.spec != 0:
		if specFileName == 'none':
			s = mat.specCol
			spec = mat.getSpec()
			s[0] = s[0]*spec;s[1] = s[1]*spec;s[2] = s[2]*spec;
			mat_str += 'color s %s %s %s %s # %s\n' % (round(s[0],5), round(s[1],5), round(s[2],5), brdf_mat, num_mats)
		else:
			mat_str += 'texture s %s %s # %s\n' % (brdf_mat, diffFileName, num_mats)
		num_mats += 1

	if bumpFileName != 'none' :
		mat_str += 'normalmap %s %s # %s\n' % (brdf_mat, bumpFileName, num_mats)
		num_mats += 1

	mat_str += 'mult %s ' % (num_mats - brdf_mat - 1)
	for i in range(brdf_mat+1, num_mats):
		mat_str += '%s ' % (i)
	mat_str += '%s # %s %s\n' % (brdf_mat, num_mats,  mat.name)
	
	mat_num.append(num_mats)
	num_mats += 1
	return mat_str, mat_num, num_mats 
	

def write(filename):
	print "\nExporting Corona Files"
	time1 = sys.time()
	if not filename.lower().endswith('.nra2'):
		filename += '.nra2'
	
	# leave edit mode
	editmode = Window.EditMode()
	if editmode: Window.EditMode(0)
		
	# get all objects
	scn = bpy.data.scenes.active
	objects = scn.objects
	
	#start filestring
	file_str = ''
	mat_str = ''
	
	# set environment
	if EnvType == 0:
		file_str += "cloudy_sky\n"
		
	elif EnvType == 1:
		sunfound = 0
		for o in objects:
			if o.type == "Lamp":
				ldata = o.getData()
				if ldata.type == 1:
					omatrix = o.matrix
					vect = Mathutils.Vector([0,0,1,0])
					svect = vect*omatrix
					svect.normalize()
					turb = Turbidity.val
					file_str += "daylight "+str(round(-svect[0],5))+" "+str(round(-svect[1],5))+" "+str(round(-svect[2],5))+" "+str(turb)+"\n"
					sunfound = 1		
		if not sunfound:
			file_str += "daylight #\n"
			print "No Sun found - using default !"
					
	elif EnvType == 2:
		file_str += "rgbeis %s %s\n" % (EnvFile.val, EnvGain.val)
	
	elif EnvType == 3:
		file_str += "black\n"

	num_mats = 0
	mat_num = []
	camn = 1
	materials = []
	
	#### print mats #####
	
	mats = bpy.data.materials
	for mat in mats:
		if mat.users <= 0: continue 
		if mat.name not in materials:
			materials.append(mat.name)
		mat_str, mat_num, num_mats = writemat(mat, mat_str, mat_num, num_mats, filename)
		
	file_str += '%s\n%s' % (num_mats, mat_str)
	numObjects = 0
#	objnum = len(objects)
#	objc = 0
	ob_str = ''
	curr_mat = 0
	shaders = []
	expmat = []
	LAYERS = scn.layers
	
	for o in objects:
		
# 		Progressbar 2x slower !!
#		objc += 1
#		onum = float(objnum-objc)/float(objnum)
#		Window.DrawProgressBar(1-onum, "Obj - %s" %o.name)
		
		if o.type == "Camera" and o.Layers & 1: # o.layers[0] in LAYERS:
			
			#size
			context = Scene.GetCurrent().getRenderingContext()
			xsize = context.imageSizeX()
			ysize = context.imageSizeY()

			#aspect ratio
			if xsize > ysize:
				ratio = float(xsize)/float(ysize)
			else:
				ratio = 1.0
			#lens
			odata = o.getData()
			lensp = odata.getLens()
			lens = lensp#*1.2439024
			#transform
			omatrix = o.getInverseMatrix()
			
			oquat = omatrix.toQuat()
			
			nquat = Quaternion()
			nquat.w = 0;nquat.x = 0;nquat.y = 1;nquat.z = 0
			
			mquat = Quaternion()
			
			mquat.w = nquat.w*oquat.w - nquat.x*oquat.x - nquat.y*oquat.y - nquat.z*oquat.z;
			mquat.x = nquat.x*oquat.w + nquat.w*oquat.x + nquat.y*oquat.z - nquat.z*oquat.y;
			mquat.y = nquat.w*oquat.y - nquat.x*oquat.z + nquat.y*oquat.w + nquat.z*oquat.x;
			mquat.z = nquat.w*oquat.z + nquat.x*oquat.y - nquat.y*oquat.x + nquat.z*oquat.w;
			
			#write cam file
			camfname = Blender.sys.splitext(filename)[0] + "0" + str(camn) + '.cam'
			camn += 1
			fstopv = Fstop.val-2
			shutterv = Shutter.val
			tdist = FocusDistance.val
			fwidth = FilmWidth.val
			chromatic_aberration = ChromAber.val*490000
			tca = TChromAber.val
			iso = FilmIso.val
			camfile = open(camfname, 'wb')
			camfile.write(struct.pack("I", 0))
			camfile.write(struct.pack("fff", o.loc[0]*10, o.loc[1]*10, o.loc[2]*10))
			camfile.write(struct.pack("ffff",-mquat.w, mquat.x, mquat.y, mquat.z))
			camfile.write(struct.pack("f", 1))
			camfile.write(struct.pack("ii", xsize, ysize))
			camfile.write(struct.pack("fff", 0,0,0))
			camfile.write(struct.pack("ff", chromatic_aberration,tca))
			camfile.write(struct.pack("f", iso))
			camfile.write(struct.pack("fff", 0,0,0))
			camfile.write(struct.pack("fffffffff", 0,0,0,0,0,0,0,0,0))
			camfile.write(struct.pack("f", tdist*10))
			camfile.write(struct.pack("f", 1))
			camfile.write(struct.pack("f", 1.0))
			camfile.write(struct.pack("ff", fwidth*ratio, 0.24))
			camfile.write(struct.pack("i", fstopv))
			camfile.write(struct.pack("f", (lens/100)))
			camfile.write(struct.pack("f", 5.0)) # chromatic aberration
			camfile.write(struct.pack("i", shutterv))
			
			camfile.close()
			
			
		if o.type == "Mesh" and o.Layers & 1: # o.layers[0] in LAYERS:
			
			if ExportGeom.val:
				try:
					tmpMesh.Get('-dummy-')
				except:
					tmpMesh = Mesh.New('-dummy-')
				
				matrix = o.getMatrix()
				tmpMesh.getFromObject(o, 0, 1)
				tmpMesh.transform(matrix, 1, 0)
		
			mats = o.data.materials + o.getMaterials()
			
			# write .ra2, .uv and .n file for each material
			for m in mats:
				if m.users<=0: continue 
		
				mfilename = m.name.replace("/", "")
				mfilename = mfilename.replace(" ", "")
				mfilename = mfilename.replace(".", "")
				mfilename += "%04d"%Blender.Get("curframe")
								
				matlen = len(materials)

				for i in range(0, matlen):
					if m.name == materials[i]:
						curr_mat = mat_num[i]
										
				basefilename = sys.dirname(filename) + "/" + mfilename
				# print basefilename+".ra2"
				
				if m.name not in expmat:
					# open new files TODO: append frame number!
					ob_str += "%s %s.ra2\n" % (curr_mat, mfilename)
					expmat.append(m.name)
					numObjects += 1
					if ExportGeom.val:
						meshfile = open(basefilename + ".ra2", "wb")
						normfile = open(basefilename + ".n", "wb")
						uvfile = open(basefilename + ".uv", "wb")
				
				else:
					# append the data to files
					if ExportGeom.val:
						meshfile = open(basefilename + ".ra2", "ab")
						normfile = open(basefilename + ".n", "ab")
						uvfile  = open(basefilename + ".uv", "ab")
				if ExportGeom.val:
					faceuv = tmpMesh.faceUV
					faces = tmpMesh.faces
					if mfilename == "hair":
						for edge in tmpMesh.edges:
							meshfile.write(struct.pack("fff", edge.v1.co[0]*10, edge.v1.co[1]*10, edge.v1.co[2]*10))
							meshfile.write(struct.pack("fff", edge.v2.co[0]*10, edge.v2.co[1]*10, edge.v2.co[2]*10))
							meshfile.write(struct.pack("fff", edge.v2.co[0]*10, edge.v2.co[1]*10, edge.v2.co[2]*10))
					
					else:
						for face in faces:
						
							# this seems to be broken in blender 247+ svn!! :(
							if mats[face.mat] != m:continue
								
							lenfaces = len(face.v)
							# if lenfaces == 2: print "line!"
							for lfn in range(2, lenfaces):
								for i in [0,lfn-1,lfn]:
									meshfile.write(struct.pack("fff", face.v[i].co[0]*10, face.v[i].co[1]*10, face.v[i].co[2]*10))
									if faceuv: uvfile.write(struct.pack("ff",face.uv[i][0], face.uv[i][1]))
									else: uvfile.write(struct.pack("ff",0,0))
																						
								if face.smooth:
									for i in [0,lfn-1,lfn]: 
										normfile.write(struct.pack("fff",face.v[i].no[0], face.v[i].no[1], face.v[i].no[2]))
								else:
									for i in [0,lfn-1,lfn]:
										normfile.write(struct.pack("fff",face.no[0], face.no[1], face.no[2]))
			
					meshfile.close()
					normfile.close()
					uvfile.close()
			if ExportGeom.val:
				tmpMesh.verts=None
						
	file_str += str(numObjects) + '\n'
	file_str += ob_str
	
	matfile = open(filename, "wb") 
	matfile.write(file_str)
	matfile.close()
	
	# restore mode
	Window.EditMode(editmode)
	time2 = sys.time()
	time = time2-time1
	print "Exported in %s seconds" % (int(time))




# =======================================================================
#  start single frame export/render
# =======================================================================

def render_single(filename):
	# run corona
	rmode="none"
	if ExecuteCorona.val:
		if RenderMode.val:
			if MLTexp.val:
				rmode = "mlt_exp"
			else:
				rmode = "mlt"
		else:
			rmode = "bdpt"
	frame = Blender.Get("curframe")
	render_anim(filename,frame,frame+1,rmode)




# =======================================================================
#  start batch mode animation render
# =======================================================================

def render_anim(filename,start_frame,end_frame,mode):
	# render info
	scale = float(ScaleSize.val)/100
	render_xsize = SizeX.val*scale
	render_ysize = SizeY.val*scale
	samples = NumSamples.val
	running = 0
	for frame in range(start_frame, end_frame):
		Blender.Set("curframe", frame)
		framefilename = Blender.sys.splitext(filename)[0] + "%04d"%frame + '.nra2'
		write(framefilename)
		if running == 1:
			renderthread.join()
			running = 0
			os.system("rm "+sys.dirname(filename) + "/*%04d.ra2"%(frame-1))
			os.system("rm "+sys.dirname(filename) + "/*%04d.n"%(frame-1))
			os.system("rm "+sys.dirname(filename) + "/*%04d.uv"%(frame-1))
			os.system("rm "+sys.dirname(filename) + "/*%04d.nra2"%(frame-1))
			os.system("rm "+sys.dirname(filename) + "/*%04d??.cam"%(frame-1))

		# start render thread on this nra2 file
		renderthread = renderer("%scorona %s %s -s %d -w %d -h %d "%(coronadir,mode,framefilename,samples,render_xsize,render_ysize))
		if mode != "none":
			running = 1
			renderthread.start()
	# if running == 1:
		# renderthread.join()
		# running = 0

def render_anim_callback(filename):
	scn = bpy.data.scenes.active
	context = scn.getRenderingContext()
	start_frame = context.startFrame()
	end_frame = context.endFrame()
	if RenderMode.val: mode = "batch"
	else : mode = "batch-bdpt"
	render_anim(filename, start_frame, end_frame+1, mode)




# =======================================================================
#   GUI
# =======================================================================

# gui raster
margin = 6
col0 = margin
col1 = 80 + margin
col2 = 160 + margin
col3 = 240 + margin
cwd  = 75
row0 = 150
row1 = 130
row2 = 110
row3 =  90
row4 =  70
row5 =  50
rhg = 18
loff = 7
	
# Assign event numbers to buttons
evtNoEvt	= 0
evtExport   = 1
evtExportAnim   = 2
openCamera = 3
openEnv = 4
openRSet = 5
evtgp = 6
evtloadimg = 7
evtchangesize = 8
evtFocusS = 9
evtFocusC = 10
evtFocusT = 11
evtgplane = 12
evtccache = 13
openMaterials = 14
evtsetMat = 15
evtrpreview = 16
evtrrpreview = 17
evtselectmat = 18

ExecuteCorona = Draw.Create(1)
DefaultExport = Draw.Create(1)
RenderMode = Draw.Create(1)
MLTexp = Draw.Create(1)

NumSamples = Draw.Create(100)

Scene.GetCurrent().getRenderingContext().imageSizeX(480)
Scene.GetCurrent().getRenderingContext().imageSizeY(320)

sceneSizeX = Scene.GetCurrent().getRenderingContext().imageSizeX()
sceneSizeY = Scene.GetCurrent().getRenderingContext().imageSizeY()

SizeX = Draw.Create(sceneSizeX)
SizeY = Draw.Create(sceneSizeY)

strScaleSize = "Scale Size %t | 200 % %x200 | 175 % %x175 | 150 % %x150 | 125 % %x125 | 100 % %x100 | 75 % %x75 | 50 % %x50 | 25 % %x25"
ScaleSize = Draw.Create(100)

#ExpotGeom
ExportGeom = Draw.Create(1)

#Cam Properties
FilmWidth = Draw.Create(36.0)
ChromAber = Draw.Create(0.0)
TChromAber = Draw.Create(0.0)
FocusDistance = Draw.Create(10.0)

EnvGain = Draw.Create(1.0)
Turbidity = Draw.Create(2.0)
GroundPlane = Draw.Create(0)

## Exposure
strExpo = "Shutter %t | 60s %x0 | 30s %x1 | 15s %x2 | 8s %x3 | 4s %x4 | 2s %x5 | 1s %x6 | 1/2s %x7 | 1/4s %x8 | 1/8s %x9 | 1/15s %x10 | 1/30s %x11 | 1/60s %x12 | 1/125s %x13 | 1/250s %x14 | 1/500s %x15 | 1/1000s %x16 | 1/2000s %x17 | 1/4000s %x18 | 1/8000 %x19"
Shutter = Draw.Create(14)
strFstop = "f-stop %t | 0.5 %x0 | 0.7 %x1 | 1.0 %x2 | 1.4 %x3 | 2.0 %x4 | 2.8 %x5 | 4.0 %x6 | 5.6 %x7 | 8 %x8 | 11 %x9 | 16 %x10 | 22 %x11 | 32 %x12 | 45 %x13 | 64 %x14 | 90 %x15 | 128 %x16"
Fstop = Draw.Create(8)
FilmIso = Draw.Create(100)

## Environment Type
EnvType = Draw.Create(3)

try:
	EnvFile = Draw.Create(EnvFile.val)
except:
	EnvFile = Draw.Create("")
	

strEnvType = "Env Type %t | Cloudy Sky %x0 | Physical Sky %x1 | HDRI Map %x2 | none %x3"

selmat = Draw.Create(0)
submat = Draw.Create(0)
mattype = Draw.Create(0)
diffuse = Draw.Create(0.0,0.0,0.0)
dmulti = Draw.Create(0.0)
specular = Draw.Create(0.0,0.0,0.0)
smulti = Draw.Create(0.0)
emulti = Draw.Create(0.0)
sexpo = Draw.Create(0.0)
ior = Draw.Create(0.0)
archi = Draw.Create(0)
sssmulti = Draw.Create(0.0)
ssscol = Draw.Create(0.0,0.0,0.0)
absorb = Draw.Create(0)
dispers = Draw.Create(0.0)
absr = Draw.Create(0.0)
absg = Draw.Create(0.0)
absb = Draw.Create(0.0)

##########################################################
# load/save other settings to blender scene
##########################################################
def readWriteValue(dict, name, draw, write, default):
	if write == 0:
		try:
			try:
				draw.val = dict[name]
			except KeyError:
				draw.val = default
		except ValueError:
			print "ValueError: on property \"%s\""%name
	else:
		dict[name] = draw.val

def readwriteSettings(write=0):
	global DefaultExport, RenderMode, MLTexp, ScaleSize
	global Fstop, FocusDistance
	global FilmIso, Shutter, FilmWidth, ChromAber, TChromAber
	global GroundPlane
	global EnvType, EnvFile, EnvGain, Turbidity,ExecuteCorona
	try:
		d = Blender.Scene.GetCurrent().properties['corona'].convert_to_pyobject()
	except KeyError:
		d = {}
	readWriteValue(d, 'DefaultExport', DefaultExport, write, 1)
	readWriteValue(d, 'RenderMode', RenderMode, write, 1)
	readWriteValue(d, 'MLTexp' , MLTexp, write, 1)
	readWriteValue(d, 'Exportgeom', ExportGeom, write, 1)
	readWriteValue(d, 'sizex', SizeX, write, 480)
	readWriteValue(d, 'sizey', SizeY, write, 320)
	readWriteValue(d, 'scale size', ScaleSize, write, 100)
	readWriteValue(d, 'film width', FilmWidth, write, 36)
	readWriteValue(d, 'chromatic aberration', ChromAber, write, 0.0)
	readWriteValue(d, 'transverse chromatic aberration', TChromAber, write, 0.0)
	readWriteValue(d, 'fstop', Fstop, write, 8)
	readWriteValue(d, 'focusdistance', FocusDistance, write, 10.0)
	readWriteValue(d, 'filmiso', FilmIso, write, 100)
	readWriteValue(d, 'shutter', Shutter, write, 14)
	readWriteValue(d, 'envtype', EnvType, write, 3)
	readWriteValue(d, 'envfile', EnvFile, write, None)
	readWriteValue(d, 'envgain', EnvGain, write, 1.0)
	readWriteValue(d, 'turbidity', Turbidity, write, 2.0)
	readWriteValue(d, 'samples', NumSamples, write, 100)
	readWriteValue(d, 'ExecuteCorona', ExecuteCorona, write, 1)
	
	if write != 0:
		Blender.Scene.GetCurrent().properties['corona'] = d
## read Registry ##
readwriteSettings()

def setFocus(target):
	loc2 = 0
	global FocusDistance
	
	scn = Scene.GetCurrent()
	cam = scn.objects.camera
	objects = scn.objects
	loc1 = cam.getLocation()
	if target == "S":
		selObj = Object.GetSelected()[0]
		try:	
			loc2 = selObj.getLocation()
		except:
			print "select an object to focus\n"
			
	if target == "C":
		loc2 = Window.GetCursorPos()
		
	if target == "T":
		try:
			for t in objects:
				if t.name == (cam.name+".target"):
					loc2 = t.loc
		except:
			print "No target found using Cursor\n"
			loc2 = Window.GetCursorPos()
			
	FocusDistance.val = ((((loc1[0]-loc2[0])**2)+((loc1[1]-loc2[1])**2)+((loc1[2]-loc2[2])**2))**0.5)


####################################################################
##### Render preview

def render_preview(mat, mode):
	filename = "%sdefault.nra2" % tmpdir
	mat_pre = ""
	mat_num_pre = []
	num_mats_pre = 19
	mat_exp = "cloudy_sky\n"
	mat_pre,mat_num_pre,num_mats_pre = writemat(mat, mat_pre, mat_num_pre, num_mats_pre, tmpdir+"previews/")
	mat_exp += str(num_mats_pre) + "\n"
	mat_exp += """diffuse # 0
normals  0 # 1
color d 0.39244 0.39244 0.39244 0 # 2
mult 2 1 2 0 # 3 back
ashi 980 980 # 4
normals  4 # 5
color e 29471.96364 29471.96364 29471.96364 4 # 6
color d 1.0 1.0 1.0 4 # 7
color s 0.5 0.5 0.5 4 # 8
mult 4 5 6 7 8 4 # 9 Emitter_li
diffuse # 10
normals  10 # 11
color e 38313.54976 32878.32076 21139.72733 10 # 12
color d 1.0 0.85814 0.55176 10 # 13
mult 3 11 12 13 10 # 14 Emitter_re
diffuse # 15
normals  15 # 16
texture d 15 greychecker.tga # 17
mult 2 16 17 15 # 18 plane
"""
	mat_exp += mat_pre
	
	mat_exp +="""4
14 Emitter_re.ra2
9 Emitter_li.ra2
%d probe.ra2
18 plane.ra2"""%(num_mats_pre-1)
	
	os.system("rm "+tmpdir+"previews/preview.txt")
	prefilename = tmpdir+"previews/preview.nra2"
	preimage = open(prefilename, 'wb')
	
	preimage.write(mat_exp)
	preimage.close()
	
	xsize=ysize=100
	samples = 50
	if mode == "preview2": samples = 400
	
	renderthread = renderer("%scorona %s %s -s %d -w %d -h %d "%(coronadir,mode,prefilename,samples,xsize,ysize))
	renderthread.start()
	renderthread.join()
	imgfname = tmpdir+"previews/"+activemat.name+".tga"
	os.system("mv "+tmpdir+"previews/preview.tga "+imgfname)
	key = activemat.name+".tga"
	try:
		bpy.data.images[key].reload()
		bpy.data.images[key].fakeUser = 1
	except:
		try:
			tmpimage = Image.Load(imgfname)
			bpy.data.images[key].fakeUser = 1
		except:
			print "could not load preview image %s!"%key
	# os.system("rm "+tmpdir+"previews/*.tga")
	os.system("rm "+imgfname)
	
####################################################################

######  Draw Camera  ###############################################
def drawCamera():
   
	global evtNoEvt, evtExport, evtExportAnim, evtFocusS, evtFocusC, openCamera, openEnv , openRSet, evtchangesize
	global SizeX, SizeY, strScaleSize, ScaleSize, Fstop, FocusDistance, GroundPlane, FilmWidth, ChromAber, TChromAber
	global Shutter, FilmIso

	drawButtons()
	
	BGL.glColor3f(0.9,0.9,0.9)
	BGL.glRectf(col0*sx-2,180*sy,(col0+cwd)*sx+1,204*sy)
	
	BGL.glRasterPos2i(int(col0*sx),int(157*sy)); stext = Draw.Text("Shutter Speed : ","small")
	Shutter = Draw.Menu(strExpo, evtNoEvt,int(col1*sx), int(row0*sy), int(cwd*sx), int(rhg*sy), Shutter.val, "Set exposure time.")
	BGL.glRasterPos2i(int(col0*sx),int(137*sy)); ftext = Draw.Text("f-stop : ","small")
	Fstop = Draw.Menu(strFstop, evtNoEvt, int(col1*sx), int(row1*sy), int(cwd*sx), int(rhg*sy), Fstop.val, "Defines the f-Stop. Larger f-stop means more depth of field.")
	BGL.glRasterPos2i(int(col0*sx),int(117*sy)); itext = Draw.Text("Film ISO : ","small")
	FilmIso = Draw.Number("", evtNoEvt,int(col1*sx),int(row2*sy),int(cwd*sx),int(rhg*sy),FilmIso.val,25,6400, "Set ISO speed.")
	BGL.glRasterPos2i(int(col0*sx),int(97*sy)); ftext = Draw.Text("Focus Distance : ","small")
	FocusDistance = Draw.Number("", evtNoEvt, int(col1*sx), int(row3*sy), int(cwd*sx), int(rhg*sy), FocusDistance.val, 0.0, 10000, "Distance from the camera at which objects will be in focus.")
	Draw.Button("S", evtFocusS, int(col2*sx), int(row3*sy), int(cwd/3*sx), int(rhg*sy), "Get the distance from the selected object.")
	Draw.Button("C", evtFocusC, int((col2+cwd/3)*sx), int(row3*sy), int(cwd/3*sx), int(rhg*sy), "Get the distance from the 3d cursor.")
	Draw.Button("T", evtFocusT, int((col2+2*cwd/3)*sx), int(row3*sy), int(cwd/3*sx), int(rhg*sy), "Get the distance from camera target.")
	BGL.glRasterPos2i(int(190*sx),int(157*sy)); fwtext = Draw.Text("Film Width: ","small")
	FilmWidth = Draw.Number("", evtNoEvt, int(col3*sx), int(row0*sy), int(cwd*sx), int(rhg*sy), FilmWidth.val, 0.0, 100.0, "Width of the film in mm.")
	BGL.glRasterPos2i(int(190*sx),int(137*sy)); fwtext = Draw.Text("C/A: ","small")
	ChromAber = Draw.Number("", evtNoEvt, int(col3*sx), int(row1*sy), int(cwd*sx), int(rhg*sy), ChromAber.val, 0.0, 1.0, "Chromatic aberration (0 - no, 1.0 - a lot).")
	BGL.glRasterPos2i(int(190*sx),int(117*sy)); fwtext = Draw.Text("T C/A: ","small")
	TChromAber = Draw.Number("", evtNoEvt, int(col3*sx), int(row2*sy), int(cwd*sx), int(rhg*sy), TChromAber.val, .0, 1.0, "Extra transverse chromatic aberration (0-1.0).")
	
	BGL.glColor3f(0.9,0.9,0.9) ; BGL.glRasterPos2i(int(col0*sx),int((row4+loff)*sy)) ; Draw.Text("Size:")
	SizeX = Draw.Number("X: ", evtchangesize, int(col1*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), SizeX.val, 1, 4096, "Width of the render.")
	SizeY = Draw.Number("Y: ", evtchangesize, int(col2*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), SizeY.val, 1, 3072, "Height of the render.")
	ScaleSize = Draw.Menu(strScaleSize, evtNoEvt, int(col3*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), ScaleSize.val, "Scale image size.")

def drawEnv():
####################################################
	global evtNoEvt, evtExport, evtExportAnim, evtFocusS, evtFocusC, openCamera, openEnv, openRSet, evtloadimg, evtgplane
	global strEnvType, EnvType, EnvFile, EnvGain, Turbidity 
	
	drawButtons()
	
	BGL.glColor3f(0.9,0.9,0.9)
	BGL.glRectf(col1*sx-2,180*sy,(col1+cwd)*sx,204*sy)
	
	EnvType = Draw.Menu(strEnvType, evtNoEvt, int(col1*sx), int(150*sy), int(cwd*sx), int(rhg*sy), EnvType.val, "Set the enviroment type.")
	if EnvType.val == 2:
		EnvFile = Draw.String("Probe: ", evtNoEvt, int(col1*sx), int(row1*sy), int((col2-col0+cwd)*sx), int(rhg*sy), EnvFile.val, 100, "File name of the hdr probe.")
		EnvGain = Draw.Number("Gain: ", evtNoEvt, int(col3*sx), int(row0*sy), int(cwd*sx), int(rhg*sy), EnvGain.val, 0.001, 1000.00, "Gain.")
		Draw.Button("Load", evtloadimg, int(col3*sx), int(row2*sy), int(cwd*sx),int(rhg*sy),"Load environment map.")
		
		if EnvFile.val:
			try:
				pimage = Image.Load(EnvFile.val)
				Draw.Image(pimage,int(col1*sx), int(row5*sy),sx*(col2-col1+cwd)/float(pimage.getSize()[0]),sy*78/float(pimage.getSize()[1]))
			except:
				BGL.glColor3f(0.9,0.2,0.2)
				BGL.glRasterPos2i(int(col1*sx),int((row2+loff)*sy)); ftext = Draw.Text("File not found! ","normal")
			
	if EnvType.val == 1:
		Turbidity = Draw.Number("Sky turbidity", evtNoEvt, int(col1*sx), int(row1*sy), int((col2-col1+cwd)*sx), int(rhg*sy), Turbidity.val, 1.5, 5.0, "Sky turbidity.")
		
	GroundPlane = Draw.Button("Ground Plane", evtgp, int(col2*sx), int(row0*sy), int(cwd*sx), int(rhg*sy), "Place large ground plane at (0,0,0).")
	
def drawRSet():
####################################################
	global RenderMode, MLTexp, NumSamples
	global evtNoEvt, evtExport, evtExportAnim, evtFocusS, evtFocusC, openCamera, openEnv, openRSet, evtloadimg, evtgplane
	
	drawButtons()
	
	BGL.glColor3f(0.9,0.9,0.9)
	BGL.glRectf(col2*sx-2,182*sy,(col2+cwd)*sx+1,204*sy)
	
	RenderMode = Draw.Toggle("MLT", evtNoEvt, int(col1*sx), int(row0*sy), int(cwd*sx), int(rhg*sy), RenderMode.val, "Enable metropolis sampling.")
	if RenderMode.val:
		MLTexp = Draw.Toggle("bidir", evtNoEvt, int(col2*sx), int(row0*sy), int(cwd*sx), int(rhg*sy), MLTexp.val, "Enable bidirectional MLT.")
	NumSamples = Draw.Number("samples", evtNoEvt, int(col1*sx), int(row1*sy), int((col2-col1+cwd)*sx), int(rhg*sy), NumSamples.val, 1, 100000000, "Nnumber of samples per pixel for batch mode.")
	
def drawMat():
####################################################-textures still missing !
	global materials, selmat, selmats, submat, activemat, diffuse, dmulti, specular, smulti, sexpo, ior, ssscol, sssmulti
	global emulti, mattype, aobject, activeobmats, omats, mmats
	global evtNoEvt, evtExport, evtExportAnim, evtFocusS, evtFocusC, openCamera, openEnv, openRSet, evtloadimg
	global evtgplane, archi, dispers, absorb, absr, absg, absb
	
	drawButtons()
	
	BGL.glColor3f(0.9,0.9,0.9)
	BGL.glRectf(col3*sx-2,180*sy,(col3+cwd)*sx,204*sy)
	
	mats = bpy.data.materials

	for mat in mats:
		if mat.users == 0:
			continue 
		if mat.name not in materials:
			materials.append(mat.name)
		
	matlen = len(materials)
	
	strmats = "Scenematerials %t |"
	for i in range(matlen):
		strmats += "| "+(materials[i]+" %x"+str(i))
	strmats += "| No Material Assigned %x"+str(matlen+1)
	
	strtypes = "Type %t | Diffuse %x0 | Ashi %x1 | Glass %x2 | Emitter %x3 | Milkglass %x4"
	
	activeobmats = []
	activeobmat = 0
		
	aobject = bpy.data.scenes.active.objects.active
	if aobject != None and aobject.type == "Mesh":
		omats = aobject.getMaterials()
		mmats = aobject.getData(mesh=0).materials
		activeobmats = (mmats + omats)
	else:
		return
		print "hmmmm"
		activemat = Material.Get(materials[selmat.val])
		amat = activemat.name

	if len(activeobmats)!=0:	
		activeobmat = activeobmats[submat.val]

	mcount = len(activeobmats)
	
	strsubs = "Submaterial %t |"
	io = 0
	im = 0
	for i in range(mcount):
		if activeobmats[i] in omats:
			strsubs += "| object %x"+str(i)
		else:
			strsubs += "| me- "+(str(im+1)+" %x"+str(i))
			im += 1
			
	if not activeobmat:
		selmat.val = matlen+1
		amat = "No Material Assigned"
		
	if activeobmat:
		for i in range(matlen):
			if activeobmat.name == materials[i]:
				selmat.val = i
				activemat = Material.Get(materials[selmat.val])
				amat = activemat.name
			
	dreggn = Draw.Button("get", evtselectmat, int((col2+cwd-10)*sx), int(row0*sy), int(10*sx), int(rhg*sy), "Edit material of selected object.")
	selmats = Draw.String("", evtsetMat, int(col2*sx), int(150*sy), int((cwd-10)*sx), int(rhg*sy), amat, 100, "Materialname")
	# selmat = Draw.Menu(strmats, evtsetMat, int((col2+cwd-10)*sx), int(row0*sy), int(10*sx), int(rhg*sy), selmat.val, "Select Material")
	
	if activeobmat:
		submat = Draw.Menu(strsubs, evtNoEvt, int(col3*sx), int(row0*sy), int((cwd-20)*sx), int(rhg*sy), submat.val, "Select Submaterial")
		Draw.Label("of  "+str(mcount), int((col3+cwd-20)*sx),int((row0)*sy),int(cwd*sx),int(rhg*sy))
		
		mattype = Draw.Menu(strtypes, evtsetMat, int(col2*sx), int(row1*sy), int(cwd*sx), int(rhg*sy), mattype.val, "Set material type.")
	
		mdcol = activemat.rgbCol
		mdmulti = activemat.ref
		mscol = activemat.specCol
		msmulti = activemat.spec
		msexpo = activemat.hard*20
		mior = activemat.IOR
		mensss = activemat.enableSSS
		mssscol = activemat.sssCol
		msssmulti = activemat.sssScale
		memulti = activemat.emit*10000
		vdisp = 1.0-activemat.specTransp
			
		ar = activemat.sssRadiusRed
		ag = activemat.sssRadiusGreen
		ab = activemat.sssRadiusBlue
		
		if mattype.val == 0: # diffuse
			
			diffuse = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row2*sy), int(cwd*sx), int(rhg*sy), (mdcol[0],mdcol[1],mdcol[2]))
			dmulti = Draw.Number("Diffuse", evtsetMat, int(col3*sx), int(110*sy), int(cwd*sx), int(rhg*sy),mdmulti, 0, 1)
			
		elif mattype.val == 1: # ashi
			
			diffuse = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row2*sy), int(cwd*sx), int(rhg*sy), (mdcol[0],mdcol[1],mdcol[2]))
			dmulti = Draw.Number("Diffuse", evtsetMat, int(col3*sx), int(row2*sy), int(cwd*sx), int(rhg*sy),mdmulti, 0, 1)
			specular = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row3*sy), int(cwd*sx), int(rhg*sy), (mscol[0],mscol[1],mscol[2]))
			smulti = Draw.Number("Specular", evtsetMat, int(col3*sx), int(row3*sy), int(cwd*sx), int(rhg*sy), msmulti, 0, 1)
			sexpo = Draw.Number("Exponent", evtsetMat, int(col3*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), msexpo, 20, 2000)
			
		elif mattype.val == 2: # glass
			
			archi = Draw.Toggle("archi", evtsetMat, int((col3+cwd/2)*sx),int(row1*sy),int(cwd/2*sx),int(rhg*sy), archi.val, "Disable refraction.")
			if archi.val:
				diffuse = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row2*sy), int(cwd*sx), int(rhg*sy), (mdcol[0],mdcol[1],mdcol[2]))
				dmulti = Draw.Number("Transmission", evtsetMat, int(col3*sx), int(row2*sy), int(cwd*sx), int(rhg*sy),mdmulti, 0, 1)
				specular = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row3*sy), int(cwd*sx), int(rhg*sy), (mscol[0],mscol[1],mscol[2]))
				smulti = Draw.Number("Specular", evtsetMat, int(col3*sx), int(row3*sy), int(cwd*sx), int(rhg*sy), msmulti, 0, 1)
				sexpo = Draw.Number("Exponent", evtsetMat, int(col3*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), msexpo, 20, 2000)
			else:
				absorb = Draw.Toggle("Medium", evtsetMat, int(col3*sx),int(row1*sy),int(cwd/2*sx),int(rhg*sy), mensss, "Enable homogeneous medium.")
				ior = Draw.Number("IOR ", evtsetMat,int(col2*sx),int(row2*sy),int(cwd*sx),int(rhg*sy), mior, 0.1, 10, "Index of refraction at 400nm.")
				dispers = Draw.Number("Dispersion ", evtsetMat,int(col3*sx),int(row2*sy),int(cwd*sx),int(rhg*sy), vdisp, 0.001, 1.0, "Difference to lower IOR at 700nm.")
				if absorb.val:
					absr = Draw.Number("r", evtsetMat, int(col2*sx), int(row3*sy), int(cwd*sx/3), int(rhg*sy), ar, 0.0, 1000.0, "Expected path length in medium (red).")
					absg = Draw.Number("g", evtsetMat, int((col2+cwd/3)*sx), int(row3*sy), int(cwd*sx/3), int(rhg*sy), ag, 0.0, 1000.0, "Expected path length in medium (green).")
					absb = Draw.Number("b", evtsetMat, int((col2+2*cwd/3)*sx), int(row3*sy), int(cwd*sx/3), int(rhg*sy), ab, 0.0, 1000.0, "Expected path length in medium (blue).")
					ssscol = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), (mssscol[0],mssscol[1],mssscol[2]), "Volume scattering color. Set to black for Beer's law only.")
					sssmulti = Draw.Number("scale", evtsetMat, int(col3*sx), int(row3*sy), int(cwd*sx), int(rhg*sy),msssmulti, 0.001, 30.0, "Expected path length multiplier.")
				
		elif mattype.val == 3:
			
			diffuse = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row2*sy), int(cwd*sx), int(rhg*sy), (mdcol[0],mdcol[1],mdcol[2]))
			emulti = Draw.Number("Emission", evtsetMat, int(col3*sx), int(row2*sy), int(cwd*sx), int(rhg*sy),memulti, 0, 20000)

		elif mattype.val == 4: # milkglass
			
			diffuse = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row2*sy), int(cwd*sx), int(rhg*sy), (mdcol[0],mdcol[1],mdcol[2]), "Multiply when entering the medium.")
			dmulti = Draw.Number("Diffuse", evtsetMat, int(col3*sx), int(row2*sy), int(cwd*sx), int(rhg*sy),mdmulti, 0, 1, "Diffuse color multiplier.")
			specular = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row3*sy), int(cwd*sx), int(rhg*sy), (mscol[0],mscol[1],mscol[2]), "Determines Fresnel term at surface.")
			smulti = Draw.Number("Specular", evtsetMat, int(col3*sx), int(row3*sy), int(cwd*sx), int(rhg*sy), msmulti, 0, 1, "Specular color multiplier.")
			sexpo = Draw.Number("Exponent", evtsetMat, int(col3*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), msexpo, 20, 2000, "Phong exponent.")
			absorb = Draw.Toggle("Medium", evtsetMat, int(col3*sx),int(row1*sy),int(cwd*sx),int(rhg*sy), mensss, "Enable homogeneous medium.")
			if absorb.val:
				absr = Draw.Number("r", evtsetMat, int(col2*sx), int(row5*sy), int(cwd*sx/3), int(rhg*sy), ar, 0.0, 1000.0, "Expected path length in medium (red).")
				absg = Draw.Number("g", evtsetMat, int((col2+cwd/3)*sx), int(row5*sy), int(cwd*sx/3), int(rhg*sy), ag, 0.0, 1000.0, "Expected path length in medium (green).")
				absb = Draw.Number("b", evtsetMat, int((col2+2*cwd/3)*sx), int(row5*sy), int(cwd*sx/3), int(rhg*sy), ab, 0.0, 1000.0, "Expected path length in medium (blue).")
				ssscol = Draw.ColorPicker(evtsetMat, int(col2*sx), int(row4*sy), int(cwd*sx), int(rhg*sy), (mssscol[0],mssscol[1],mssscol[2]), "Volume scattering color. Set to black for Beer's law only.")
				sssmulti = Draw.Number("scale", evtsetMat, int(col3*sx), int(row5*sy), int(cwd*sx), int(rhg*sy),msssmulti, 0.001, 30.0, "Expected path length multiplier.")
			
			
		#### Preview ###
		imgsx = sx*(col1-col0+cwd)/128.0
		imgsy = sy*(row0-row3+rhg)/128.0
		if imgsx < imgsy:
			scale = imgsx
		else:
			scale = imgsy
		try:
			dimg = Draw.Image(bpy.data.images[activemat.name+".tga"],col0*sx,row3*sy,scale,scale)
		except:
			pass
		
		bwd = scale*128.0
		Draw.Button("update", evtrpreview, int(col0*sx), int(row4*sy), int(bwd/2), int(rhg*sy), "Render quick preview.")
		Draw.Button("refine", evtrrpreview, int(col0*sx+bwd/2), int(row4*sy), int(bwd/2), int(rhg*sy), "Render high quality preview.")
	

def drawGUI():
######################################################
	global Screen
	
	if Screen==0:
		drawMat()
	if Screen==1:
		drawCamera()
 	if Screen==2:
 		drawEnv()
	if Screen==3:
 		drawRSet()


def event(evt, val):  # function that handles keyboard and mouse events
	# FIXME: stupid redraw event on mouse over to ensure consistency of material editor
	# this also kills tooltips.
	# if(Screen == 0): Draw.Redraw()
	if evt == Draw.ESCKEY or evt == Draw.QKEY:
		stop = Draw.PupMenu("OK?%t|Cancel export %x1")
		if stop == 1:
			Draw.Exit()
			return

def drawButtons():
#####################################################
	global evtExport, evtExportAnim, openCamera, openEnv, openRSet, openSSet, openTmap, evtNoEvt, ExecuteCorona, DefaultExport, ExportGeom
	global sx, sy
	scissorbox=Buffer(GL_FLOAT,4)
	glGetFloatv(GL_SCISSOR_BOX,scissorbox)
	width = float(scissorbox[2])
	height = float(scissorbox[3])
	sx = width/350.0
	sy = height/225.0
	
	BGL.glColor3f(0.1,0.1,0.1)
	BGL.glRecti(0,0,int(width),int(height))
	BGL.glColor3f(0.2,0.2,0.2)
	BGL.glRectf(0,48*sy,width,180*sy)
	BGL.glColor3f(0.9,0.9,0.9)
	BGL.glRectf(0,180*sy,width,182*sy)
	BGL.glRectf(0,48*sy,width,50*sy)
	
	####### Background Image ########
	# bgfname = blenddir + "corona_radiata.png"
	# bgimg = Image.Load(bgfname)
	# xscale = width/1200
	# yscale = height/1200
	# min = xscale
	# if yscale < min: min = yscale
	# Draw.Image(bgimg,(width-1200*min)/2,(height-1200*min)/2,min,min,0,0)

	BGL.glColor3f(0.9, 0.9, 0.9) ; BGL.glRasterPos2i(int(10*sx),int(212*sy)) ; Draw.Text("corona-6: radiata v%s"%corona_export_version, "normal")
	Draw.Button("Camera", openCamera, int(col0*sx), int(185*sy), int(cwd*sx), int(18*sy), "Open camera dialog.")
	Draw.Button("Environment", openEnv, int(col1*sx), int(185*sy), int(cwd*sx), int(18*sy), "Open environment dialog.")
	Draw.Button("Renderer", openRSet, int(col2*sx), int(185*sy), int(cwd*sx), int(18*sy), "Open render settings dialog.")
	Draw.Button("Material", openMaterials, int(col3*sx), int(185*sy), int(cwd*sx), int(18*sy), "Open material dialog")
	
	ExecuteCorona = Draw.Toggle("run", evtNoEvt, int(col1*sx), int((15+12.5)*sy), int(cwd/2*sx), int(12.5*sy), ExecuteCorona.val, "Execute corona and render the saved .nra2 file.")
	DefaultExport = Draw.Toggle("def",evtNoEvt, int((col1+cwd/2)*sx), int((15+12.5)*sy), int(cwd/2*sx), int(12.5*sy), DefaultExport.val, "Use %s/default.nra2 as filename."%tmpdir) 
	ExportGeom = Draw.Toggle("geo",evtNoEvt, int(col1*sx), int(15*sy), int(cwd/2*sx), int(12.5*sy), ExportGeom.val, "Export geometry.")
	if ExportGeom.val:
		Draw.Button("clear", evtccache, int((col1+cwd/2)*sx), int(15*sy), int(cwd/2*sx), int(12.5*sy), "Clear out %s .nra2 etc ..  not pics ;)"%tmpdir)
	Draw.Button("Render", evtExport, int(col3*sx), int(15*sy), int(cwd*sx), int(25*sy), "Render single frame.")
	Draw.Button("Animate", evtExportAnim, int(col2*sx), int(15*sy), int(cwd*sx), int(25*sy), "Render animation.")
	
	# BGL.glColor3f(0.9, 0.9, 0.9) ; BGL.glRasterPos2i(int(col2*sx),int(7*sy)) ; Draw.Text("Press Q or ESC to quit.", "tiny")
	
def setEnvMap(efilename):
	EnvFile.val = efilename
	
def buttonEvt(evt):  # function that handles button events
	
	global Screen,EnvFile
	
	if evt == evtExport:
		if DefaultExport.val:
			filename = "%sdefault.nra2" % tmpdir
			render_single(filename)
		else:
			Blender.Window.FileSelector(render_single, 'Corona Export', 'default.nra2')

	if evt == evtExportAnim:
		if DefaultExport.val:
			filename = "%sdefault.nra2" % tmpdir
			render_anim_callback(filename)
		else:
			Blender.Window.FileSelector(render_anim_callback, 'Corona Render Animation', 'default.nra2')
	
	if evt == evtchangesize:
		Scene.GetCurrent().getRenderingContext().imageSizeX(SizeX.val);
		Scene.GetCurrent().getRenderingContext().imageSizeY(SizeY.val);
		
	if evt == openMaterials:
		Screen = 0
		
	if evt == openCamera:
		Screen = 1
	
	if evt == openEnv:
		Screen = 2
		
	if evt == openRSet:
		Screen = 3
	
	if evt == evtloadimg:
		Blender.Window.FileSelector(setEnvMap, "Load Environment Map", ('.hdr'))
	
	if evt == evtFocusS:
		setFocus("S")
		
	if evt == evtFocusC:
		setFocus("C")
		
	if evt == evtFocusT:
		setFocus("T")
			 
	if evt == evtgp:
		try:
			Object.Get("GroundPlane")
		except:
			scn = bpy.data.scenes.active
			gplane = Mesh.Primitives.Plane(200)
			scn.objects.new(gplane,"GroundPlane")
			mat = Material.New("GroundPLane")
			mat.spec = 0
			matl = []
			matl.append(mat)
			gplane.materials = matl
			
	if evt == evtccache:
		print "Clearing CacheDir"
		os.system("rm "+tmpdir+"*.ra2")
		os.system("rm "+tmpdir+"*.nra2")
		os.system("rm "+tmpdir+"*.uv")
		os.system("rm "+tmpdir+"*.n")
		os.system("rm "+tmpdir+"*.tga")
		os.system("rm "+tmpdir+"*.cam")
		os.system("rm "+tmpdir+"*.ppm")
		os.system("rm "+tmpdir+"previews/*.ppm")
		os.system("rm "+tmpdir+"previews/*.tga")
		print "Done"
	
	if evt == evtrpreview:
		render_preview(activemat, "preview")
	
	if evt == evtrrpreview:
		render_preview(activemat, "preview2")
		
	if evt == evtselectmat:
		FLAGS = activemat.getMode()
		if activemat.emit != 0:
			mattype.val = 3
		elif FLAGS & Material.Modes['RAYTRANSP']:
			if activemat.glossTra < 1.0:
				mattype.val = 4
			else:
				mattype.val = 2
		elif activemat.spec != 0:
			mattype.val = 1
		else:
			mattype.val = 0

	if evt == evtsetMat:
		
		activemat.name = selmats.val

		if mattype == 0: # diffuse
			activemat.setMode(50397187)
			activemat.rgbCol = diffuse.val
			activemat.ref = dmulti.val
			activemat.enableSSS = 0
			activemat.spec = 0
			activemat.emit = 0

		elif mattype == 1: # ashi
			activemat.setMode(50397187)
			activemat.rgbCol = diffuse.val
			activemat.ref = dmulti.val
			activemat.specCol = specular.val
			activemat.spec = smulti.val
			activemat.hard = int(sexpo.val/20)
			activemat.enableSSS = 0
			activemat.emit = 0
			
		elif mattype == 2: # glass
			activemat.setMode(50528259)
			activemat.emit = 0
			activemat.glossTra = 1.0
			activemat.enableSSS = absorb.val
			if archi.val:
				activemat.IOR = 1.0
				activemat.enableSSS = 0
				activemat.rgbCol = diffuse.val
				activemat.ref = dmulti.val
				activemat.specCol = specular.val
				if activemat.spec != 0:
					activemat.spec = smulti.val
				else:
					activemat.spec = 0.04
			else:
				activemat.specTransp = 1-dispers.val
				activemat.IOR = ior.val
				
			if absorb.val:
				activemat.sssRadiusRed = absr.val
				activemat.sssRadiusGreen = absg.val
				activemat.sssRadiusBlue = absb.val
				activemat.sssCol = ssscol.val
				activemat.sssScale = sssmulti.val


		elif mattype == 3:
			activemat.setMode(50397187)
			activemat.enableSSS = 0
			activemat.rgbCol = diffuse.val
			activemat.ref = 0
			activemat.emit = emulti.val/10000

		elif mattype == 4: # milkglass
			activemat.setMode(50528259)
			activemat.emit = 0
			activemat.glossTra = 0.5
			activemat.enableSSS = absorb.val
			activemat.rgbCol = diffuse.val
			activemat.ref = dmulti.val
			activemat.specCol = specular.val
			activemat.spec = smulti.val
			activemat.hard = int(sexpo.val/20)
				
			if absorb.val:
				activemat.sssRadiusRed = absr.val
				activemat.sssRadiusGreen = absg.val
				activemat.sssRadiusBlue = absb.val
				activemat.sssCol = ssscol.val
				activemat.sssScale = sssmulti.val
		
			
	Draw.Redraw()
	readwriteSettings(1)
		
if __name__=='__main__':
	
	Screen = 1
	Draw.Register(drawGUI, event, buttonEvt)

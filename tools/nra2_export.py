#!BPY

"""
Name: 'NRA2 (.nra2) ...'
Blender: 242
Group: 'Export'
Tooltip: 'Exports to nra2 format'
"""


import math
import sys
import struct, os, time, cStringIO
import Blender
from Blender import *
from cStringIO import StringIO
from Blender.Mathutils import *


def getTextureFileName(filename):
	f = filename[filename.rfind('/')+1:-4]
	f = f[f.rfind('\\')+1:]
	return f

def write(filename):
  if not filename.lower().endswith('.nra2'):
	filename += '.nra2'

  # don't messily append to existing files
  #if sys.exists(filename): return

  # leave edit mode
  editmode = Window.EditMode()
  if editmode: Window.EditMode(0)
  
  # scale light emission
  light_scale = 30.0

  # get all objects
  objects = Blender.Object.Get()

  # and open shader file
  matfile = open(filename, "wb")
  matfile.write("cloudy_sky\n")
  file_str = StringIO()

  num_mats = 0
  materials = []
  mat_num = []
  # numObjects = 0

  num_mats = 0
  camn = 1
  for o in objects:
	
	if (o.Layer & 1 != 1): continue # ignore all objects not in layer 1
		
	if o.getType() == "Camera":
		print "Camera " + o.name
		#size
		context = Scene.GetCurrent().getRenderingContext()
		xsize = context.imageSizeX()
		ysize = context.imageSizeY()
		print xsize, ysize
		#aspect ratio
	  	#  ratio = float(context.imageSizeY())/float(context.imageSizeX())
		#lens
		odata = o.getData()
		lens = odata.getLens()
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
		
		#focus distance
		tdist = 10
		for t in objects:
			if t.name == (o.name+".target"):
				t = t.loc
				c = o.loc
				tdistv = Mathutils.Vector([(t[0]-c[0]), (t[1]-c[1]), (t[2]-c[2])])
				tdist = tdistv.length
		
		if tdist == 10: print "Free Camera"
		else: print "Target Camera\nFosucs Distance = " + str(tdist)
					
		#write cam file
		camfname = Blender.sys.splitext(filename)[0] + "0" + str(camn) + '.cam'
		camn += 1
		camfile = open(camfname, 'wb')
		camfile.write(struct.pack("I", 0))
		camfile.write(struct.pack("fff", o.loc[0], o.loc[1], o.loc[2]))
		camfile.write(struct.pack("ffff",-mquat.w, mquat.x, mquat.y, mquat.z))
		camfile.write(struct.pack("f", 20))
		camfile.write(struct.pack("ii", xsize, ysize))
		camfile.write(struct.pack("fffffffff", 0,0,0,0,0,0,0,0,0))
		camfile.write(struct.pack("fffffffff", 0,0,0,0,0,0,0,0,0))
		camfile.write(struct.pack("f", tdist))
		camfile.write(struct.pack("f", 1))
		camfile.write(struct.pack("f", 1.0))
		camfile.write(struct.pack("ff", 0.24*xsize/ysize, 0.24))#0.36, 0.24))
		camfile.write(struct.pack("i", 1))
		camfile.write(struct.pack("f", (lens/110)))
		camfile.write(struct.pack("f", 0))
		camfile.write(struct.pack("i", 40))
		
		camfile.close()
		
		
	
	if o.getType() == "Mesh":
		print "Mesh -- " + o.name
	 	 # numObjects = numObjects + 1		# count objects
	  	mats = o.data.materials + o.getMaterials()
	  	# if (len(mats) > 0):  # append all materials to the export list
	  	#   mat = mats[0]
	  	for mat in mats:
			if mat.users > 0 and mat.name not in materials:
			  materials.append(mat.name)
			  # write out material, append material number
			  # TODO:
			  # eta = mat.getFresnelMirr() 
			  # if eta == 0 or eta = 1.0:
			  # TODO: if face.smooth have normal shader.
	
			  brdf_mat = num_mats
			  # file_str.write('ashi 100 100 # ' + str(num_mats) + '\n')
			  file_str.write('diffuse # ' + str(num_mats) + '\n')
			  num_mats = num_mats + 1
			  textures = mat.getTextures()
			  diffFileName = 'none'
			  specFileName = 'none'
			  bumpFileName = 'none'
			  bumpStrength = 0
	
			  for tex in textures:
				if tex == None: continue
				if tex.tex.type != Texture.Types.IMAGE: continue
				if tex.mtCol != 0:
				  diffFileName = getTextureFileName(tex.tex.image.filename)
				if tex.mtSpec != 0:
				  specFileName = getTextureFileName(tex.tex.image.filename)
				if tex.mtNor != 0:
				  bumpFileName = getTextureFileName(tex.tex.image.filename)
				  bumpStrength = tex.norfac / 5.0		   
		
			  file_str.write('normals '+ str(brdf_mat) + ' # ' + str(num_mats) + '\n')
			  num_mats = num_mats + 1
	
			  if diffFileName == 'none' :
				c = mat.rgbCol
				file_str.write('color d '+ str(round(c[0],5))+' '+str(round(c[1],5))+' '+str(round(c[2],5))+' '+ str(brdf_mat) +' # ' + str(num_mats) + '\n')
				num_mats = num_mats + 1
			  else:
				file_str.write('texture d '+ str(brdf_mat) + ' ' + diffFileName + ' # ' + str(num_mats) + '\n')
				num_mats = num_mats + 1
	
			  if mat.emit != 0 :
				# todo: use *diffuse color for colored light sources?
				file_str.write('color e '+str(round(mat.emit*light_scale,5))+' '+str(round(mat.emit*light_scale,5))+' '+str(round(mat.emit*light_scale,5))+' '+ str(brdf_mat) +' # ' + str(num_mats) + '\n')
				num_mats = num_mats + 1
				
			  ## if specFileName == 'none' :
			  ##   c = mat.rgbCol
			  ##   file_str.write('color s '+ str(round(c[0]/15,5))+' '+str(round(c[1]/15,5))+' '+str(round(c[2]/15,5))+' '+ str(brdf_mat) +' # ' + str(num_mats) + '\n')
			  ##   num_mats = num_mats + 1
			  ## else:
			  ##   file_str.write('texture s '+ str(brdf_mat) + ' ' + diffFileName +' # ' + str(num_mats) + '\n')
			  ##   num_mats = num_mats + 1
	
			  if bumpFileName != 'none' :
				file_str.write('normalmap '+ str(brdf_mat) + ' ' + bumpFileName +' # ' + str(num_mats) + '\n')
				num_mats = num_mats + 1
	
	
			  file_str.write('mult ' + str(num_mats - brdf_mat - 1) + ' ')
			  for i in range(brdf_mat+1, num_mats):
				file_str.write(str(i) + ' ')
			  file_str.write(str(brdf_mat) + ' # ' + str(num_mats) + ' '+ mat.name + '\n')
			  mat_num.append(num_mats)
			  num_mats = num_mats + 1
			
		  #else:
		  #  print "ERROR: object " + obj.name + " has mesh without material!""
		
  matfile.write(str(num_mats) + '\n')
  matfile.write(file_str.getvalue())
  file_str = StringIO()
  numObjects = 0
  # matfile.write(str(numObjects) + '\n')
  curr_mat = 0
  shaders = []
  tmpMesh = Mesh.New('sr_tmp_export')
  for o in objects:
	if (o.Layer & 1 != 1): continue # ignore all objects not in layer 1

	if o.getType() != "Mesh":
	  continue

	meshname = o.data.name
	mo = Blender.NMesh.GetRaw(meshname)
	if not mo : continue

	mats = o.data.getMaterials() + o.getMaterials();
	matrix = o.getMatrix()
	tmpMesh.getFromObject(o);
	tmpMesh.transform(o.getMatrix(), True)

	# write .ra2, .uv and .n file for each material
	for m in mats:

	  if m.users == 0: continue

	  mfilename = m.name.replace("/", "_")

	  curr_mat = 0
	  for i in range(0, len(materials)):
		if m.name == materials[i]:
		  curr_mat = mat_num[i]
	  
	  # basefilename = filename[:-5] + "/" + sys.basename(filename)[:-5] + mfilename
	  basefilename = sys.dirname(filename) + "/" + mfilename
	  if not sys.exists(basefilename + ".ra2") :
		# matfile.write(str(curr_mat) + ' ' + filename[:-5] + mfilename + '.ra2\n')
		file_str.write(str(curr_mat) + ' ' + mfilename + '.ra2\n')
		numObjects = numObjects + 1;

	  # append the data to files
	  meshfile = open(basefilename + ".ra2", "ab");
	  normfile = open(basefilename + ".n", "ab");
	  uvfile   = open(basefilename + ".uv", "ab");

	  for face in tmpMesh.faces:

		if mats[face.mat] != m: continue
		for lfn in range(2, len(face.v)):
		  for i in [0,lfn-1,lfn]:
			meshfile.write(struct.pack("fff", face.v[i].co[0], face.v[i].co[1], face.v[i].co[2]))
						
			if  tmpMesh.faceUV: uvfile.write(struct.pack("ff",face.uv[i][0], face.uv[i][1]))
			else: uvfile.write(struct.pack("ff",0,0))
							
		  # TODO: if material uses interpolated normals
		  if face.smooth:
			for i in [0,lfn-1,lfn]: 
			  normfile.write(struct.pack("fff",face.v[i].no[0], face.v[i].no[1], face.v[i].no[2]))
		  else:
			for i in [0,lfn-1,lfn]:
			  normfile.write(struct.pack("fff",face.no[0], face.no[1], face.no[2]))

	  meshfile.close()
	  normfile.close()
	  uvfile.close()


  matfile.write(str(numObjects) + '\n')
  matfile.write(file_str.getvalue())
  matfile.close()

  # restore mode
  Window.EditMode(editmode)

def main():
	Blender.Window.FileSelector(write, 'NRa2 Export', 'dreggn.nra2')

if __name__=='__main__':
	print "\nExporting Corona Files"
  	main()

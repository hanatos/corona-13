#!BPY

"""
Name: 'Target Camera'
Blender: 244
Group: 'AddMesh'
Tooltip: 'Create a Target Camera'
"""

__author__ = "Gregor Quade"
__url__ = ("http://www.r3ndr.com/mysite","blender", "elysiun")
__version__ = "1.0"
__bpydoc__ = """\

This script creates a targetcamera
"""

from Blender import Camera,Object,Scene,Constraint,Window,Draw,Mathutils

scn = Scene.GetCurrent()

convert = "TargetCam?%t|Convert Active Camera %x1|Create New Camera %x2"
result = Draw.PupMenu(convert)

if result == 1:
	camob = scn.objects.camera

else:
	cam = Camera.New('persp', 'TCam')
	camob = scn.objects.new(cam)
	scn.objects.camera = camob
	camob.setLocation(Window.GetCursorPos())

if camob:
	
	target = scn.objects.new('Empty', (camob.name + '.target'))
	matrix = camob.getMatrix()
	tvect = Mathutils.Vector([0,0,-10])
	targetv = tvect*matrix
	target.setLocation(targetv)
	
	cam = camob.data
	#cam.drawLimits = 1
	cam.drawSize = 2.0
	campos = Mathutils.Vector([camob.dloc[0],camob.dloc[1],camob.dloc[2]])
	dofdistv = campos-targetv
	dofdist = dofdistv.length
	#cam.dofDist = dofdist
	cam.clipEnd = 1500
	const = camob.constraints.append(Constraint.Type.TRACKTO)
	const[Constraint.Settings.TARGET] = target
	const[Constraint.Settings.TRACK] = Constraint.Settings.TRACKNEGZ
	const[Constraint.Settings.UP] = Constraint.Settings.UPY
	scn.update()
else:
	block = []
	block.append("")
	msg = Draw.PupBlock("No Camera in Scene !!", block)
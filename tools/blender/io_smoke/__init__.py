import bpy
import io
import os
import array
import io
import struct

from math import *

# -------------------------------------------------------------------------------

bl_info = {
    "name": "export smoke ptc",
    "description": "export the currently active object to smoke ptc if supported. have to start from unit cube [-1,1]",
    "author": "Johannes Meng",
    "version": (1, 0),
    "blender": (2, 74, 0),
    "location": "File > Export",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Import-Export"
}

def write_volume(context, obj, base_path, file_name):
    assert(obj)

    if len(obj.modifiers) == 0:
        return

    file = io.open(os.path.join(base_path, file_name), 'wb')

    for mod in obj.modifiers:
        if mod.type == 'SMOKE':
            res = mod.domain_settings.domain_resolution
            mul = 1
            if mod.domain_settings.use_high_resolution:
                mul = mod.domain_settings.amplify + 1
            hires = [mul * res[0], mul*res[1], mul*res[2]]
            sizes = [res[0], res[1], res[2], mul]
            num_cells = res[0]*res[1]*res[2]
            hires_num_cells = res[0]*res[1]*res[2]*mul*mul*mul
            arr = array.array('I', sizes)
            arr.tofile(file)
            del arr
            density = mod.domain_settings.density_grid
            # Albedo? This is RGBA, but we only use RGB.
            # color = [v for (i, v) in enumerate(mod.domain_settings.color_grid) if i % 4 != 3]
            # Temperature, convert to kelvin between ignition and max from gui:
            # flame = [1000 * (mod.domain_settings.flame_ignition + (mod.domain_settings.flame_max_temp - mod.domain_settings.flame_ignition) * f) for f in mod.domain_settings.flame_grid]
            flame = mod.domain_settings.flame_grid
            # custom blender extension:
            # vel = mod.domain_settings.velocity_grid #[:num_cells*3]

            assert(len(density) == hires_num_cells)
            assert(len(flame) == hires_num_cells)
            # assert(len(vel) == 3 * num_cells)

            # cell_size = max(mod.domain_settings.cell_size*obj.scale)
            cell_size = max(p*q for p,q in zip(mod.domain_settings.cell_size,obj.scale)) 
            arr = array.array('f', [cell_size]) # original bbox was cube from [-1,1], account for double width
            arr.tofile(file)
            del arr
            # write transform. we already scaled (voxel size), so
            # all we need is translation + rotation:
            arr = array.array('f', [obj.location[0]+obj.bound_box[0][0]*obj.scale[0], obj.location[1]+obj.bound_box[1][0]*obj.scale[1], obj.location[2]+obj.bound_box[2][0]*obj.scale[2]])
            arr.tofile(file)
            del arr
            assert(obj.rotation_mode == 'XYZ')
            arr = array.array('f', [obj.rotation_euler[0], obj.rotation_euler[1], obj.rotation_euler[2]])
            arr.tofile(file)
            del arr

            arr = array.array('f', density)
            arr.tofile(file)
            del arr
            arr = array.array('f', flame)
            arr.tofile(file)
            del arr
            # arr = array.array('f', vel)
            # arr.tofile(file)
            file.close()   

    
def write_smoke_ptc(context, base_path, file_name):
    obj = context.active_object
    if not obj:
        return {'CANCELLED'}
    for f in range(bpy.context.scene.frame_start,bpy.context.scene.frame_end):
        bpy.context.scene.frame_set(frame=f)
        fname=file_name+"%04d"%f
        write_volume(context, obj, base_path, fname)
    return {'FINISHED'}

# -------------------------------------------------------------------------------

# ExportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ExportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy.types import Operator

class ExportSmoke(Operator, ExportHelper):
    """Export the currently active object as a .ptc file"""
    bl_idname = "export_scene.smoke_ptc"  # important since its how bpy.ops.import_test.some_data is constructed
    bl_label = "Export active object as somke ptc (.ptc)" 
    # ExportHelper mixin class uses this
    filename_ext = ".ptc"

    filter_glob = StringProperty(
            default="*.ptc",
            options={'HIDDEN'},
            )

    def invoke(self, context, event):
        if not context.active_object:
            self.report({'WARNING'}, "Please activate an object.")
            return {'CANCELLED'}
        return ExportHelper.invoke(self, context, event)

    def execute(self, context):
        status = write_smoke_ptc(context, os.path.dirname(self.filepath), os.path.basename(self.filepath))
        if status == {'CANCELLED'}:
            self.report({'WARNING'}, "Please select an object.")
        return status

# -------------------------------------------------------------------------------

# Only needed if you want to add into a dynamic menu
def menu_func_export(self, context):
    self.layout.operator(ExportSmoke.bl_idname, text="Smoke (*.ptc)")


def register():
    # Register the exporter.
    bpy.utils.register_class(ExportSmoke)
    bpy.types.INFO_MT_file_export.append(menu_func_export)
    
def unregister():
    bpy.utils.unregister_class(ExportSmoke)
    bpy.types.INFO_MT_file_export.remove(menu_func_export)

if __name__ == "__main__":
    register()

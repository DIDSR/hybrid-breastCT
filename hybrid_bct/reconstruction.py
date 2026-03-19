from __future__ import annotations

import numpy as np
import tigre

from copy import deepcopy

def fxn_alter_geometry(geo, volume_with_calcs, CORoffset_mm_3D, new_vx_um, vertical_offset_mm, zstart_mm):
    geo_raytracing = deepcopy(geo)

    # Object parameters: adjusted for high resolution cropped volume
    geo_raytracing.nVoxel = np.array([volume_with_calcs.shape[0], volume_with_calcs.shape[1], volume_with_calcs.shape[2]])
    geo_raytracing.dVoxel = np.array([new_vx_um * 0.001, new_vx_um * 0.001, new_vx_um * 0.001])  # in mm
    geo_raytracing.sVoxel = geo_raytracing.nVoxel * geo_raytracing.dVoxel  # total size of the object in mm
    geo_raytracing.offOrigin = np.array([(geo_raytracing.nVoxel[0] * geo_raytracing.dVoxel[0]) / 2 + vertical_offset_mm + zstart_mm, CORoffset_mm_3D[1], CORoffset_mm_3D[2]])
    #geo_raytracing.offOrigin = np.array([(geo_raytracing.nVoxel[0] * geo_raytracing.dVoxel[0]) / 2 + vertical_offset_mm, CORoffset_mm_3D[1], CORoffset_mm_3D[2]])

    print("\nGeometry for ray tracing: ",flush=True)
    print(f"Distance from source to detector (mm): {geo_raytracing.DSD:.2f}",flush=True)
    print(f"Distance from source to object (mm): {geo_raytracing.DSO:.2f}",flush=True)
    print(f"Mode: {geo_raytracing.mode} beam",flush=True)
    print(f"Detector Matrix (px): [{geo_raytracing.nDetector[0]}, {geo_raytracing.nDetector[1]}]",flush=True)
    print(f"Detector Size (mm): [{geo_raytracing.sDetector[0]:.2f}, {geo_raytracing.sDetector[1]:.2f}]",flush=True)
    print(f"Detector Pixel Size (mm): [{geo_raytracing.dDetector[0]:.3f}, {geo_raytracing.dDetector[1]:.3f}]",flush=True)
    print(f"Voxel Size (mm): [{geo_raytracing.dVoxel[0]:.4f}, {geo_raytracing.dVoxel[1]:.4f}, {geo_raytracing.dVoxel[2]:.4f}]",flush=True)
    print(f"Object Matrix Size (vx): [{geo_raytracing.nVoxel[0]}, {geo_raytracing.nVoxel[1]}, {geo_raytracing.nVoxel[2]}]",flush=True)
    print(f"Object Image Dimensions (mm): [{geo_raytracing.sVoxel[0]:.2f}, {geo_raytracing.sVoxel[1]:.2f}, {geo_raytracing.sVoxel[2]:.2f}]",flush=True)
    print(f"Offset from Origin (mm): [{geo_raytracing.offOrigin[0]:.2f}, {geo_raytracing.offOrigin[1]:.2f}, {geo_raytracing.offOrigin[2]:.2f}]",flush=True)
    print(f"Offset of Detector (mm): [{geo_raytracing.offDetector[0]:.2f}, {geo_raytracing.offDetector[1]:.2f}]",flush=True)

    return geo_raytracing

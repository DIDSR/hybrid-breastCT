import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob
import math
import tigre
import tigre.algorithms as algs
import tigre.utilities.gpu as gpu
from scipy.interpolate import RegularGridInterpolator

def fxn_read_energy_spectrum(loc_spectra, energy_bin_size_keV):
    # Open and read the file
    with open(loc_spectra, 'r') as file:
        lines = file.readlines()

    energy_ev = []
    num_photons = []

    # Process each line
    for line in lines:
        if line.startswith('\ufeff#'):
            continue
        if line.startswith('#') or not line.strip():
            continue  # Skip comments and empty lines
        data = list(map(float, line.split()))
        energy_ev.append(data[0])
        num_photons.append(data[1])

    # Convert to numpy arrays
    energy_ev = np.array(energy_ev)
    num_photons = np.array(num_photons)

    # Create energy levels based on bin size
    energy_levels_eV = np.arange(energy_ev[0], energy_ev[-1] + 1, energy_bin_size_keV * 1000)

    # Bin the photon fluence data into the specified energy levels
    photon_fluence = np.histogram(energy_ev, bins=np.append(energy_levels_eV, np.inf), weights=num_photons)[0]

    return energy_levels_eV, photon_fluence

def fxn_read_material_file(filename, density_object, energy_levels_eV=None):
    # Open and read the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    energy_mat_eV = []
    total_mfp = []
    nom_density = None
    read_density_value_next = False

    # Process each line
    for line in lines:
        if line.strip().startswith('#'):
            if read_density_value_next:
                nom_density = float(line.split('#')[1].strip())
                read_density_value_next = False
                continue
            if 'NOMINAL DENSITY (g/cm^3)' in line:
                read_density_value_next = True
            continue  # Skip comments

        data = line.split()
        if len(data) >= 6:
            energy_mat_eV.append(float(data[0]))
            total_mfp.append(float(data[4]))
        elif len(data) == 2:
            energy_mat_eV.append(float(data[0]))
            total_mfp.append(float(data[1]))

    # Calculate mass attenuation coefficients and mu values
    energy_mat_eV = np.array(energy_mat_eV)
    total_mfp = np.array(total_mfp)
    mass_att_coeff = 1 / (total_mfp * nom_density)
    mu_cm = mass_att_coeff * density_object

    # Optionally interpolate mu values at specified energy levels
    if energy_levels_eV is not None:
        interp_mu_cm = np.interp(energy_levels_eV, energy_mat_eV, mu_cm)
        return energy_levels_eV, interp_mu_cm

    return energy_mat_eV, mu_cm

def fxn_load_seg(scanID, root_dir):
    print(f'Loading seg {scanID:04d}... ', end='',flush=True)

    # Build the file path to the first reconstructed image to read the header
    filename = os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/CTi_Corrected/CT{scanID:04d}_01.0001')
    with open(filename, 'rb') as file:
        header = np.fromfile(file, dtype='float32', count=20)

    # Extract values from the header
    Nz = int(header[14])  # Adjust index based on actual header layout if needed
    matrixsize = int(header[0])
    original_vx_um = header[9] * 1000  # Convert mm to um

    # Define the volume dimensions
    Nxyz = [Nz, matrixsize, matrixsize]
    seg_volume = np.zeros(Nxyz, dtype=np.uint8)

    # Load segmentation slices
    for islice in range(1, Nz + 1):
        slice_filename = os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/SEG/SEG{scanID:04d}_01.{islice:04d}')
        with open(slice_filename, 'rb') as file:
            IMslice = np.fromfile(file, dtype='uint8', count=matrixsize*matrixsize)
            IMslice = IMslice.reshape(matrixsize, matrixsize)
            seg_volume[islice - 1, : , :] = IMslice

    # Remap segmentation labels
    seg_volume[seg_volume == 4] = 1  # 4 is initially some voxels surrounding skin
    seg_volume[seg_volume == 3] = 2  # 2 and 3 are sparse and dense fibroglandular

    # Define labels
    labels = {
        'air': 0,
        'adipose': 1,
        'glandular': 2,
        'skin': 5
    }

    print('done.',flush=True)
    return seg_volume, Nxyz, original_vx_um, labels

def fxn_load_projections_and_geometry(scanID, root_dir, xsiz, ysiz, dexel_mm, vertical_offset_mm):
    print(f'\nLoading metadata and projection data for scan {scanID}',flush=True)

    # Define paths
    path_prj_folder = os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/SIN/')

    # Count  number of projections acquired
    nprj_OG = len(glob(os.path.join(path_prj_folder, 'dg_*')))

    # Read in first and last frames that should be used for recon
    frames_file = os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/XXX/frames{scanID:04d}.xxx')
    with open(frames_file, 'r') as file:
        line = file.readline()
        first_last_frame = np.array([int(num) for num in line.split()])

    nframes = first_last_frame[1] - first_last_frame[0] + 1
    nframes_OG = first_last_frame[1] - first_last_frame[0] + 1

    # Read in I_0
    izeros_file = os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/XXX/IZreo.txt')
    izeros = np.loadtxt(izeros_file)

    # Read angles file which has angles and projections used for recon
    angles_file = os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/XXX/angles{scanID:04d}.xxx')
    ang_number, deg = [], []
    with open(angles_file, 'r') as file:
        for line in file:
            columns = list(map(float, line.split()))
            deg.append(columns[-1])
            ang_number.append(columns[0])
    ang_OG = np.deg2rad(deg)
    ang_number = np.array(ang_number, dtype=int)
    nangles = len(ang_OG)

    # Adjust the number of frames if the number of angles is less than the original number of frames
    if nangles < nframes_OG:
        nframes = nangles

    print(f'# projections used for recon: {nframes}',flush=True)

    ## Read in recon parameters to know how much to crop raw projection
    recon_parameters = []
    with open(os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/XXX/recon_plan{scanID:04d}.xxx')) as file:
        for line in file:
            recon_parameters.append(int(line.strip()))
    recon_parameters = np.array(recon_parameters)
    xsiz_new = recon_parameters[3] #width of recon in pixels
    nSlices = recon_parameters[4]
    prjstack= np.zeros((nframes, ysiz, xsiz_new), dtype= 'float32')

    icount = 0  # Python uses 0-based indexing, adjust as necessary
    if nangles < nframes_OG:
        print(f'Will use projection {ang_number[0]} - {ang_number[-1]}.',flush=True)

        for iprj in ang_number:
            prj_filename = os.path.join(path_prj_folder, f'dg_{scanID:04d}.{iprj:04d}')
            if os.path.exists('{}dg_{:04d}.{:04d}'.format(path_prj_folder,scanID,iprj)):
                raw_prj= np.fromfile(open('{}dg_{:04d}.{:04d}'.format(path_prj_folder,scanID,iprj), 'rb'), dtype= 'float32')
                raw_prj= np.reshape(raw_prj,(ysiz, xsiz))
                raw_prj_crop = raw_prj[:,int(0+((xsiz-xsiz_new)/2)):int(xsiz-((xsiz-xsiz_new)/2))]
                #Log normalize with izero (XXX/IZreo.xxx)
                prjstack[icount,:, :] = -np.log(raw_prj_crop / izeros[icount])
            else:
                print(f'ERROR: Projection {iprj} does not exist',flush=True)
            icount += 1

            ang = ang_OG
    else:
        print(f'Will use projection {first_last_frame[0]} - {first_last_frame[1]}.',flush=True)
        for iprj in range(first_last_frame[0], first_last_frame[1] + 1):
            prj_filename = os.path.join(path_prj_folder, f'dg_{scanID:04d}.{iprj:04d}')
            if os.path.exists('{}dg_{:04d}.{:04d}'.format(path_prj_folder,scanID,iprj)):
                raw_prj= np.fromfile(open('{}dg_{:04d}.{:04d}'.format(path_prj_folder,scanID,iprj), 'rb'), dtype= 'float32')
                raw_prj= np.reshape(raw_prj,(ysiz, xsiz))
                raw_prj_crop = raw_prj[:,int(0+((xsiz-xsiz_new)/2)):int(xsiz-((xsiz-xsiz_new)/2))]
                #Log normalize with izero (XXX/IZreo.xxx)
                prjstack[icount,:, :] = -np.log(raw_prj_crop / izeros[icount])
            else:
                print(f'ERROR: Projection {iprj} does not exist',flush=True)
            icount += 1

    if ang_number[0] > first_last_frame[0]:
        ang_idx_start = ang_number[0] - first_last_frame[0]
    else:
        ang_idx_start = first_last_frame[0] - ang_number[0]
    ang_idx_end = ang_idx_start + nframes
    ang = ang_OG[ang_idx_start:ang_idx_end]

    print(f"Size of prjstack: [{prjstack.shape[0]}, {prjstack.shape[1]}, {prjstack.shape[2]}]",flush=True)
    if len(ang) != nframes:
        print('ERROR, # angles not equal to # projections',flush=True)

    # Replace NAN values with 0
    nan_mask = np.isnan(prjstack) | np.isinf(prjstack)
    prjstack[nan_mask] = 0

    # Read in calibration factors
    # Define the directory path based on scanID
    cal_files = glob(os.path.join(root_dir, f'Doheny_CT/CTD{scanID:04d}/CAL/', 'Dcalfactors*.cal'))

    if cal_files:
        cal_factors = []
        # Assuming the calibration factors are saved in a simple text format
        with open(cal_files[0], 'r') as file:
            for line in file:
                cal_factors.append(float(line.strip()))
        cal_factors = np.array(cal_factors)
        u_o = cal_factors[0]
        v_o = cal_factors[1]
        DSO = cal_factors[3]  # Source to object distance
        DSD = cal_factors[4]  # Source to detector distance
    else:
        print("Geometric calibration file 'Dcalfactor*.cal' not found",flush=True)

    # Find GPU
    listGpuNames = gpu.getGpuNames()
    if len(listGpuNames) == 0:
        print("Error: No gpu found",flush=True)
    else:
        for id in range(len(listGpuNames)):
            print("GPU in use: {}: {}".format(id, listGpuNames[id]),flush=True)

    gpuids = gpu.getGpuIds(listGpuNames[0])

    # ------ Set geometry -------------------------------------
    binning_mode=    1 #1 means 1024 x 1024 x nSlices, 2 means 2048 x 2048 x 2*nSlices
    geo=             tigre.geometry()
    geo.nVoxel =     np.array([nSlices*binning_mode, 1024*binning_mode, 1024*binning_mode])      # num voxels, vx
    geo.mode=        "cone"                          # or 'parallel'
    geo.DSD=         DSD                            # source-to-detector, mm
    geo.DSO=         DSO                            # source-to-object, mm
    geo.nDetector=   np.array([ysiz, xsiz_new])             # detector size, px
    geo.sDetector=   np.array([ysiz*dexel_mm, xsiz_new*dexel_mm]) #FOV at detector in mm
    geo.dDetector=   geo.sDetector/geo.nDetector
    geo.dVoxel=      np.array([(geo.sDetector[1]*(DSO/DSD))/(binning_mode*1024),(geo.sDetector[1]*(DSO/DSD))/(binning_mode*1024),(geo.sDetector[1]*(DSO/DSD))/(binning_mode*1024)]) #recon voxel size, isotropic. compute using FOV at isocenter (use mag factor) divided by matrix size
    geo.sVoxel=      geo.nVoxel*geo.dVoxel[0]       #size of image (mm)
    geo.accuracy=    0.5                            # samples/vx; recommended <= 0.5
    geo.COR=         0.0                            # y-displacement for COR, mm

    u_offset = -1; v_offset = 0;
    geo.offOrigin=   np.array([((geo.nVoxel[0]*geo.dVoxel[0])/2)+vertical_offset_mm, 0, 0])          # center of the phantom offset wrt COR, mm   #[offset in chest to nipple Z-direction,offset in vertical direction in recon coronal plane, offset in horizontal direction in recon coronal plane]
    # don't need to change offOrigin if your breast phantom is centered already in mc-gpu'
    geo.offDetector= np.array([v_offset+((ysiz/2)-v_o)*dexel_mm , ((xsiz/2)-u_o)*dexel_mm + u_offset])          # detector offset in uv-plane along u-direction,
    geo.rotDetector= np.array([cal_factors[2],0,0]) # [x: inplane detector tilt, y: rot around u axis, z: rot around v axis]
    print(f"Recon voxel size will be (mm): {geo.dVoxel}",flush=True)
    return prjstack, geo, ang

def fxn_crop_volume(iscan, volume, original_vx_um, Nxyz, VcropLocs, FLAGchestwall):
    rowstart = int(VcropLocs[iscan, 0])
    rowend = int(VcropLocs[iscan, 1])

    colstart = int(VcropLocs[iscan, 2])
    colend = int(VcropLocs[iscan, 3])

    zstart = int(VcropLocs[iscan, 4])
    zend = int(VcropLocs[iscan, 5])

    if FLAGchestwall == 1:
        volume_cropped = volume[ :zend, rowstart:rowend+1, colstart:colend+1] #[z, y, x]
        shift_mm_3D = [0, rowstart * (0.001 * original_vx_um), colstart * (0.001 * original_vx_um)] #[Z, Y, X] for tigre
    else: #cropped a smaller volume that does not include the chest wall
        volume_cropped = volume[zstart:zend+1, rowstart:rowend+1, colstart:colend+1] #[z, y, x]
        shift_mm_3D = [zstart * (0.001 * original_vx_um), rowstart * (0.001 * original_vx_um), colstart * (0.001 * original_vx_um)] #[Z, Y, X] for tigre
    Nxyz_cropped = volume_cropped.shape
    print(f"Cropped volume dimensions: {Nxyz_cropped[0]} x {Nxyz_cropped[1]} x {Nxyz_cropped[2]}",flush=True)
    return volume_cropped, shift_mm_3D

def fxn_upsample_volume(seg_volume, original_vx_um, new_vx_um):
    Nxyz = seg_volume.shape
    upsample_rate = original_vx_um / new_vx_um
    print(f'\nOriginal voxel size = {round(original_vx_um)} um.',flush=True)
    print(f'New voxel size      = {round(new_vx_um)} um.',flush=True)

    Nxyz_HR = [round(upsample_rate * Nxyz[0]), round(upsample_rate * Nxyz[1]), round(upsample_rate * Nxyz[2])]
    print(f'Forming upsampled cropped volume with dimensions: {Nxyz_HR[0]} x {Nxyz_HR[1]} x {Nxyz_HR[2]} ... ',flush=True)

    # Define the original coordinates
    z_orig = np.linspace(1, Nxyz[0], Nxyz[0], dtype=np.float32)
    y_orig = np.linspace(1, Nxyz[1], Nxyz[1], dtype=np.float32)
    x_orig = np.linspace(1, Nxyz[2], Nxyz[2], dtype=np.float32)

    # Define the new coordinates for the upsampled volume
    z_new = np.linspace(1, Nxyz[0], Nxyz_HR[0], dtype=np.float32)
    y_new = np.linspace(1, Nxyz[1], Nxyz_HR[1], dtype=np.float32)
    x_new = np.linspace(1, Nxyz[2], Nxyz_HR[2], dtype=np.float32)

    # Create a meshgrid for the new coordinates
    z_mesh, y_mesh, x_mesh = np.meshgrid(z_new, y_new, x_new, indexing='ij')

    # Perform the interpolation
    interp = RegularGridInterpolator((z_orig, y_orig, x_orig), seg_volume, method='nearest')
    seg_volume_HR = interp(np.array([z_mesh.flatten(), y_mesh.flatten(), x_mesh.flatten()]).T).reshape(z_mesh.shape)

    return seg_volume_HR

def fxn_upsample_volume_in_sections(seg_volume, original_vx_um, new_vx_um):
    Nxyz = seg_volume.shape
    upsample_rate = original_vx_um / new_vx_um
    print(f'\nOriginal voxel size = {round(original_vx_um)} um.',flush=True)
    print(f'New voxel size      = {round(new_vx_um)} um.',flush=True)

    Nxyz_HR = [round(upsample_rate * Nxyz[0]), round(upsample_rate * Nxyz[1]), round(upsample_rate * Nxyz[2])]
    print(f'Forming upsampled cropped volume with dimensions: {Nxyz_HR[0]} x {Nxyz_HR[1]} x {Nxyz_HR[2]} ... ',flush=True)

    slices = 50  # Adjust this value based on your available RAM
    sections = list(range(0, Nxyz[0], slices)) + [Nxyz[0]]
    num_sections = len(sections) - 1

    Nz = Nxyz[0]
    Ny = Nxyz[1]
    Nx = Nxyz[2]
    Nz_HR_original = round(upsample_rate * slices)
    Ny_HR = round(upsample_rate * Ny)
    Nx_HR = round(upsample_rate * Nx)
    seg_volume_HR = np.zeros(Nxyz_HR, dtype=np.uint8)

    print('Section ', end='',flush=True)
    for section in range(num_sections):
        print(f'{section+1}/{num_sections}, ', end='',flush=True)
        starting_slice = sections[section]
        ending_slice = sections[section+1]

        subvolume = seg_volume[starting_slice:ending_slice, :, :]

        # Define the original coordinates
        z_orig = np.linspace(1, subvolume.shape[0], subvolume.shape[0], dtype=np.float32)
        y_orig = np.linspace(1, Ny, Ny, dtype=np.float32)
        x_orig = np.linspace(1, Nx, Nx, dtype=np.float32)

        # Define the new coordinates for the upsampled volume
        Nz_HR = round(upsample_rate * subvolume.shape[0])
        z_new = np.linspace(1, subvolume.shape[0], Nz_HR, dtype=np.float32)
        y_new = np.linspace(1, Ny, Ny_HR, dtype=np.float32)
        x_new = np.linspace(1, Nx, Nx_HR, dtype=np.float32)

        # Create a meshgrid for the new coordinates
        z_mesh, y_mesh, x_mesh = np.meshgrid(z_new, y_new, x_new, indexing='ij')

        # Perform the interpolation
        interp = RegularGridInterpolator((z_orig, y_orig, x_orig), subvolume, method='nearest')
        upsampled_volume = interp(np.array([z_mesh.flatten(), y_mesh.flatten(), x_mesh.flatten()]).T).reshape(z_mesh.shape)

        if section < num_sections-1:
            seg_volume_HR[Nz_HR_original*section:Nz_HR_original*section+Nz_HR, :, :] = upsampled_volume
        else:
            # Last section: need to correct for rounding errors
            diff = upsampled_volume.shape[0] - seg_volume_HR[Nz_HR_original*section:Nz_HR_original*section+Nz_HR, :, :].shape[0]
            upsampled_volume = upsampled_volume[:(upsampled_volume.shape[0] - diff), :, :]
            seg_volume_HR[(Nz_HR_original*section):((Nz_HR_original*section)+Nz_HR), :, :] = upsampled_volume

    print('done.',flush=True)
    return seg_volume_HR

def fxn_getVOIcenters(scanID, voi_size_mm, num_SPvois_perbreast, num_SAvois_perbreast, input_folder):
    txt_filepath = os.path.join(input_folder, f"{voi_size_mm:02d}mm_voi_dim/{scanID:04d}_voicenters_mm.txt")

    if os.path.exists(txt_filepath):
        voicenters_mm_uncropped = np.loadtxt(txt_filepath)
        # Change order to match TIGRE [z, y, x]
        voicenters_mm_uncropped = voicenters_mm_uncropped[:, [2, 1, 0]]

        num_viable_voicenters = voicenters_mm_uncropped.shape[0]

        if num_viable_voicenters < (num_SPvois_perbreast + num_SAvois_perbreast):
            OG_num_SPvois_perbreast = num_SPvois_perbreast
            OG_num_SAvois_perbreast = num_SAvois_perbreast
            num_SPvois_perbreast = num_viable_voicenters // 2
            num_SAvois_perbreast = num_viable_voicenters // 2
            print(f"\nCannot generate {OG_num_SPvois_perbreast} signal present and {OG_num_SAvois_perbreast} signal absent VOIs because only {num_viable_voicenters} viable VOI centers.",flush=True)
            print(f"New # signal present VOIs: {num_SPvois_perbreast}    New # signal absent VOIS: {num_SAvois_perbreast}",flush=True)

        # Randomly shuffle viable VOI centers
        np.random.shuffle(voicenters_mm_uncropped)

        # Assign first num_SPvois_perbreast as signal-present VOI centers
        voi_centers_mm_SP = voicenters_mm_uncropped[:num_SPvois_perbreast]
        # Assign next num_SAvois_perbreast as signal-absent VOI centers
        voi_centers_mm_SA = voicenters_mm_uncropped[num_SPvois_perbreast:num_SPvois_perbreast + num_SAvois_perbreast]

        return voi_centers_mm_SP, voi_centers_mm_SA, num_SPvois_perbreast, num_SAvois_perbreast
    else:
        print(f"File {txt_filepath} does not exist in the input folder.",flush=True)

def fxn_insert_calc_cluster_new(volume_with_calcs, seg_volume_HR, true_x, true_y, true_z, logical_sphere_center, cluster_center, num_calcs, calc, half_dim, labels):
    calcs_added = 0
    while calcs_added < num_calcs:
        idx = np.random.randint(len(true_x))
        random_distance_from_cluster_center = np.array([true_z[idx], true_y[idx], true_x[idx]]) - logical_sphere_center
        calc_center = cluster_center - random_distance_from_cluster_center

        z_start = calc_center[0] - half_dim
        z_end   = calc_center[0] + (half_dim if calc.shape[0] % 2 == 0 else half_dim+1)
        y_start = calc_center[1] - half_dim
        y_end   = calc_center[1] + (half_dim if calc.shape[1] % 2 == 0 else half_dim+1)
        x_start = calc_center[2] - half_dim
        x_end   = calc_center[2] + (half_dim if calc.shape[2] % 2 == 0 else half_dim+1)

        # Check if the indices are within the valid range
        if z_start < 0 or z_end > volume_with_calcs.shape[0] or \
           y_start < 0 or y_end > volume_with_calcs.shape[1] or \
           x_start < 0 or x_end > volume_with_calcs.shape[2]:
            print(f"Skipping calc at index {idx} due to out of bounds indices",flush=True)
            continue

        if volume_with_calcs[z_start:z_end, y_start:y_end, x_start:x_end].shape != calc.shape:
            print('Shape mismatch',flush=True)
        volume_with_calcs[z_start:z_end, y_start:y_end, x_start:x_end] += calc
        calcs_added += 1

    volume_with_calcs[(seg_volume_HR == labels['glandular']) & (volume_with_calcs != 0)] = 101
    volume_with_calcs[(seg_volume_HR == labels['adipose']) & (volume_with_calcs != 0)] = 100

    return volume_with_calcs

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

def fxn_convert_mat_to_txt(input_folder, output_folder, VscanID):
    if not os.path.exists(output_folder):
        print('Output folder does not exist',flush=True)
        os.makedirs(output_folder)

    for scan_id in VscanID:
        mat_filename = f"{scan_id}_voicenters_mm.mat"
        txt_filename = f"{scan_id}_voicenters_mm.txt"

        mat_filepath = os.path.join(input_folder, mat_filename)
        txt_filepath = os.path.join(output_folder, txt_filename)

        if not os.path.exists(txt_filepath):
            if os.path.exists(mat_filepath):
                # Load the .mat file
                mat_contents = loadmat(mat_filepath)

                # Extract the data
                # Assuming the variable name inside the .mat file is 'voicenters'
                if 'voicenters_mm_uncropped' in mat_contents:
                    voicenters = mat_contents['voicenters_mm_uncropped']  # Adjust based on actual variable name
                else:
                    print(f"Variable 'voicenters_mm_uncropped' not found in {mat_filename}",flush=True)
                    continue

                # Write the data to a .txt file
                with open(txt_filepath, 'w') as txt_file:

                    for item in voicenters:
                        txt_file.write(' '.join(f"{x:.3f}" for x in item) + '\n')

                print(f"Converted {mat_filename} to {txt_filename}",flush=True)
            else:
                print(f"File {mat_filename} does not exist in the input folder.",flush=True)
        else:
            print(f"File {txt_filepath} already exists.",flush=True)

#input_folder = '/home/suhyun.lyu/Projects/BreastCT/a_CalcSimulation/Matlab/VOICenters/08mm_voi_dim/'
#output_folder = '/home/suhyun.lyu/Projects/BreastCT/a_CalcSimulation/Python/VOICenters/08mm_voi_dim/'
#fxn_convert_mat_to_txt(input_folder, output_folder, VscanID)

from scipy.interpolate import interp1d

def fxn_mtf_blur(prjstack_calcs, dexel_mm, f_MTF, mtf_MTF):
    blurred_prjstack_calcs = np.zeros(prjstack_calcs.shape, dtype=prjstack_calcs.dtype)

    # Plot MTF
    f_new = np.concatenate([-np.flipud(f_MTF), f_MTF[1:]])
    mtf_new = np.concatenate([np.flipud(mtf_MTF), mtf_MTF[1:]])
    Nx = prjstack_calcs.shape[1]
    Ny = prjstack_calcs.shape[2]

    # Optional: plot MTF
    #plt.plot(f_MTF, mtf_MTF)
    #plt.xlabel('lp/mm')
    #plt.ylabel('Detector MTF')
    #plt.show()

    # Compute FT of spatial profile
    Nyq_spatialprofile = 1 / (2 * dexel_mm)
    faxis_x = np.linspace(-Nyq_spatialprofile, Nyq_spatialprofile, Nx)  # frequency range and sampling in x
    faxis_y = np.linspace(-Nyq_spatialprofile, Nyq_spatialprofile, Ny)  # frequency range and sampling in y

    # Interpolate detector MTF to sampling rate of signal
    MTF1D_x = interp1d(f_new, mtf_new, kind='linear', bounds_error=False, fill_value=0)(faxis_x)
    MTF1D_y = interp1d(f_new, mtf_new, kind='linear', bounds_error=False, fill_value=0)(faxis_y)

    # Blur in x-y, one projection at a time
    for iprj in range(prjstack_calcs.shape[0]):
        blurred_image_x = np.zeros((Nx, Ny), dtype=prjstack_calcs.dtype)
        # Blur in X direction only
        for i in range(Ny):
            signal_x = prjstack_calcs[iprj, :, i]
            signal_f = np.fft.fftshift(np.fft.fft(signal_x))
            blurred_signal_f = MTF1D_x * signal_f
            blurred_image_x[:, i] = np.real(np.fft.ifft(np.fft.ifftshift(blurred_signal_f)))

        # Blur x-blurred image in Y direction
        for j in range(Nx):
            signal_y = blurred_image_x[j, :]
            signal_f = np.fft.fftshift(np.fft.fft(signal_y))
            blurred_signal_f = MTF1D_y * signal_f
            blurred_prjstack_calcs[iprj, j, :] = np.real(np.fft.ifft(np.fft.ifftshift(blurred_signal_f)))

    return blurred_prjstack_calcs

def fxn_extract_and_save_vois(rec, recon_alg, folder_suffix, scanID, loc_save_patches, loc_save_MIPjpgs, calc_diameter_mm, cluster_diameter_mm, num_calcs, density, num_SPvois_perbreast, voi_centers_mm_SP, voi_centers_mm_SA, recon_size_mm, voi_halfdim_vx, flagHU):
    num_SAvois_perbreast = num_SPvois_perbreast

    # Create directories for this set of parameters
    loc_parameters_MIP = f'{loc_save_patches}{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_MIP'
    loc_parameters_VOI = f'{loc_save_patches}{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_VOI'
    if not os.path.exists(loc_parameters_MIP):
        os.makedirs(loc_parameters_MIP)
    if not os.path.exists(loc_parameters_VOI):
        os.makedirs(loc_parameters_VOI)

    print('Extracting VOIs... ',flush=True)
    for ivoi in range(num_SPvois_perbreast):
        voi_local_mm = voi_centers_mm_SP[ivoi]
        voi_zloc_mm = voi_local_mm[0]
        recon_voi_center = np.round(voi_centers_mm_SP[ivoi] / recon_size_mm).astype(int)  # Find the corresponding VOI center in recon volume

        voi_recon = rec[
            recon_voi_center[0] - voi_halfdim_vx: recon_voi_center[0] + voi_halfdim_vx,
            recon_voi_center[1] - voi_halfdim_vx: recon_voi_center[1] + voi_halfdim_vx,
            recon_voi_center[2] - voi_halfdim_vx: recon_voi_center[2] + voi_halfdim_vx
        ]
        mip = np.max(voi_recon, axis=0)

        # Save as .npy
        filename = f'{loc_parameters_VOI}/{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_VOI_{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SP.npy'
        np.save(filename, voi_recon)
        filename = f'{loc_parameters_MIP}/{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_MIP_{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SP.npy'
        np.save(filename, mip)

        # just save 1 example MIPs per set of parameters, per scanID
        if ivoi == 0:
            plt.figure()
            plt.imshow(mip, cmap='gray')
            plt.axis('off')
            plt.colorbar()
            if flagHU == 1:
                plt.title(f'{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip)} HU')
            else:
                plt.title(f'{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip):.2f} cm^-1')
            filename = f'{loc_save_MIPjpgs}{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_MIP_{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SP.jpg'
            plt.savefig(filename)
            plt.close()

    # Extract signal absent VOIs
    for ivoi in range(num_SAvois_perbreast):
        voi_local_mm = voi_centers_mm_SP[ivoi]
        voi_zloc_mm = voi_local_mm[0]
        recon_voi_center = np.round(voi_centers_mm_SA[ivoi] / recon_size_mm).astype(int)  # Find the corresponding VOI center in recon volume
        voi_recon = rec[
            recon_voi_center[0] - voi_halfdim_vx: recon_voi_center[0] + voi_halfdim_vx,
            recon_voi_center[1] - voi_halfdim_vx: recon_voi_center[1] + voi_halfdim_vx,
            recon_voi_center[2] - voi_halfdim_vx: recon_voi_center[2] + voi_halfdim_vx
        ]
        mip = np.max(voi_recon, axis=0)

        # Save as .npy
        filename = f'{loc_parameters_VOI}/{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_VOI_{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SA.npy'
        np.save(filename, voi_recon)
        filename = f'{loc_parameters_MIP}/{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_MIP_{voi_zloc_mm:06.1f}_{scanID:04d}_{ivoi:03d}_SA.npy'
        np.save(filename, mip)


        if ivoi == 0:
            plt.figure()
            plt.imshow(mip, cmap='gray')
            plt.axis('off')
            plt.colorbar()
            if flagHU == 1:
                plt.title(f'{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip)} HU')
            else:
                plt.title(f'{calc_diameter_mm:.2f} mm calcs | peak: {np.max(mip):.2f} cm^-1')
            filename = f'{loc_save_MIPjpgs}{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_MIP_{scanID:04d}_{ivoi:03d}_SA.jpg'
            plt.savefig(filename)
            plt.close()

import os

def fxn_write_metadata(
    loc_save_patches, loc_spectra, effective_energy_keV, loc_material_files, density,
    energy_bin_size_keV, csI_thickness_mm, vertical_offset_mm, new_vx_um, calc_shape,
    calc_diameter_mm, num_calcs, cluster_diameter_mm, clusterCenteredFLAG,
    num_SPvois_perbreast, num_SAvois_perbreast, voi_size_mm, voi_size_vx, recon_alg, folder_suffix,
    kernel, flagHU, mu_water, home_dir, geo, ang
):
    loc_parameters_MIP = f'{loc_save_patches}{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_MIP'
    loc_parameters_VOI = f'{loc_save_patches}{calc_diameter_mm:.2f}mmCalc_{cluster_diameter_mm:.1f}mmCluster_{num_calcs:02d}Calcs_CaOx_{density["calc"]:1.2f}_{recon_alg}_{folder_suffix}_VOI'

    with open(os.path.join(loc_parameters_MIP, 'METADATA_ALLPATCHES.txt'), 'w') as file:
        file.write('METADATA FOR ALL PATCHES IN THIS DIRECTORY: \n\n')
        file.write(f"Energy spectrum: {loc_spectra}\n")
        file.write(f"Effective energy (keV): {effective_energy_keV}\n")
        file.write(f"Calc material file: {loc_material_files['calc']}\n")
        file.write(f"Calc density: {density['calc']:1.2f}\n")
        file.write(f"Adipose material file: {loc_material_files['adipose']}\n")
        file.write(f"Adipose density: {density['adipose']:1.2f}\n")
        file.write(f"Glandular material file: {loc_material_files['glandular']}\n")
        file.write(f"Glandular density: {density['glandular']:1.2f}\n")
        file.write(f"csI material file: {loc_material_files['csI']}\n")
        file.write(f"csI density: {density['csI']:1.2f}\n")
        file.write("\n")
        file.write(f"Energy bin size (keV): {energy_bin_size_keV}\n")
        file.write(f"CsI thickness (mm): {csI_thickness_mm}\n")
        file.write(f"Additional vertical offset of object (mm): {vertical_offset_mm}\n")
        file.write(f"Object voxel size (mm): {new_vx_um * .001:.3f}\n")
        file.write("\n")
        file.write(f"Calc shape: {calc_shape}\n")
        file.write(f"Calc diameter (mm): {calc_diameter_mm:.2f}\n")
        file.write(f"Num calcs: {num_calcs:02d}\n")
        file.write(f"Cluster diameter (mm): {cluster_diameter_mm:.1f}\n")
        file.write(f"Cluster centered within VOI: {clusterCenteredFLAG}\n")
        file.write("\n")
        file.write(f"# signal present VOIs per breast: {num_SPvois_perbreast}\n")
        file.write(f"# signal absent VOIs per breast: {num_SAvois_perbreast}\n")
        file.write(f"VOI size in one dimension (mm): {voi_size_mm}\n")
        file.write(f"VOI size in one dimension (vx): {voi_size_vx}\n")
        file.write("\n")
        file.write(f"Reconstruction algorithm: {recon_alg}\n")
        file.write(f"Reconstruction kernel: {kernel}\n")
        file.write(f"mu or HU: {flagHU} [0: mu  1: HU]\n")
        file.write(f"Mu water at effective energy (cm^-1): {mu_water:.3f}\n")
        file.write("\n")
        file.write(f"Detector MTF: {os.path.join(home_dir, 'a_CalcSimulation', 'Doheny_DetectorMTF_2x2_0.4mm_focalspotblur.csv')}\n")
        file.write("\n")
        file.write("Geometry for projections and reconstruction: \n")
        file.write(f"Distance from source to detector (mm): {geo.DSD:.2f}\n")
        file.write(f"Distance from source to object (mm): {geo.DSO:.2f}\n")
        file.write(f"Mode: {geo.mode} beam\n")
        file.write(f"Detector Matrix (px): [{geo.nDetector[0]}, {geo.nDetector[1]}]\n")
        file.write(f"Detector Size (mm): [{geo.sDetector[0]:.2f}, {geo.sDetector[1]:.2f}]\n")
        file.write(f"Detector Pixel Size (mm): [{geo.dDetector[0]:.3f}, {geo.dDetector[1]:.3f}]\n")
        file.write(f"Reconstructed Voxel Size (mm): [{geo.dVoxel[0]:.4f}, {geo.dVoxel[1]:.4f}, {geo.dVoxel[2]:.4f}]\n")
        file.write(f"Reconstructed Image Size (vx): [{geo.nVoxel[0]}, {geo.nVoxel[1]}, {geo.nVoxel[2]}]\n")
        file.write(f"Reconstructed Image Dimensions (mm): [{geo.sVoxel[0]:.2f}, {geo.sVoxel[1]:.2f}, {geo.sVoxel[2]:.2f}]\n")
        file.write(f"Offset from Origin (mm): [{geo.offOrigin[0]:.2f}, {geo.offOrigin[1]:.2f}, {geo.offOrigin[2]:.2f}]\n")
        file.write(f"Offset of Detector (mm): [{geo.offDetector[0]:.2f}, {geo.offDetector[1]:.2f}]\n")
        file.write(f"# projection angles: {len(ang)}\n")
        file.write('\n')
    with open(os.path.join(loc_parameters_VOI, 'METADATA_ALLPATCHES.txt'), 'w') as file:
        file.write('METADATA FOR ALL PATCHES IN THIS DIRECTORY: \n\n')
        file.write(f"Energy spectrum: {loc_spectra}\n")
        file.write(f"Effective energy (keV): {effective_energy_keV}\n")
        file.write(f"Calc material file: {loc_material_files['calc']}\n")
        file.write(f"Calc density: {density['calc']:1.2f}\n")
        file.write(f"Adipose material file: {loc_material_files['adipose']}\n")
        file.write(f"Adipose density: {density['adipose']:1.2f}\n")
        file.write(f"Glandular material file: {loc_material_files['glandular']}\n")
        file.write(f"Glandular density: {density['glandular']:1.2f}\n")
        file.write(f"csI material file: {loc_material_files['csI']}\n")
        file.write(f"csI density: {density['csI']:1.2f}\n")
        file.write("\n")
        file.write(f"Energy bin size (keV): {energy_bin_size_keV}\n")
        file.write(f"CsI thickness (mm): {csI_thickness_mm}\n")
        file.write(f"Additional vertical offset of object (mm): {vertical_offset_mm}\n")
        file.write(f"Object voxel size (mm): {new_vx_um * .001:.3f}\n")
        file.write("\n")
        file.write(f"Calc shape: {calc_shape}\n")
        file.write(f"Calc diameter (mm): {calc_diameter_mm:.2f}\n")
        file.write(f"Num calcs: {num_calcs:02d}\n")
        file.write(f"Cluster diameter (mm): {cluster_diameter_mm:.1f}\n")
        file.write(f"Cluster centered within VOI: {clusterCenteredFLAG}\n")
        file.write("\n")
        file.write(f"# signal present VOIs per breast: {num_SPvois_perbreast}\n")
        file.write(f"# signal absent VOIs per breast: {num_SAvois_perbreast}\n")
        file.write(f"VOI size in one dimension (mm): {voi_size_mm}\n")
        file.write(f"VOI size in one dimension (vx): {voi_size_vx}\n")
        file.write("\n")
        file.write(f"Reconstruction algorithm: {recon_alg}\n")
        file.write(f"Reconstruction kernel: {kernel}\n")
        file.write(f"mu or HU: {flagHU} [0: mu  1: HU]\n")
        file.write(f"Mu water at effective energy (cm^-1): {mu_water:.3f}\n")
        file.write("\n")
        file.write(f"Detector MTF: {os.path.join(home_dir, 'a_CalcSimulation', 'Doheny_DetectorMTF_2x2_0.4mm_focalspotblur.csv')}\n")
        file.write("\n")
        file.write("Geometry for projections and reconstruction: \n")
        file.write(f"Distance from source to detector (mm): {geo.DSD:.2f}\n")
        file.write(f"Distance from source to object (mm): {geo.DSO:.2f}\n")
        file.write(f"Mode: {geo.mode} beam\n")
        file.write(f"Detector Matrix (px): [{geo.nDetector[0]}, {geo.nDetector[1]}]\n")
        file.write(f"Detector Size (mm): [{geo.sDetector[0]:.2f}, {geo.sDetector[1]:.2f}]\n")
        file.write(f"Detector Pixel Size (mm): [{geo.dDetector[0]:.3f}, {geo.dDetector[1]:.3f}]\n")
        file.write(f"Reconstructed Voxel Size (mm): [{geo.dVoxel[0]:.4f}, {geo.dVoxel[1]:.4f}, {geo.dVoxel[2]:.4f}]\n")
        file.write(f"Reconstructed Image Size (vx): [{geo.nVoxel[0]}, {geo.nVoxel[1]}, {geo.nVoxel[2]}]\n")
        file.write(f"Reconstructed Image Dimensions (mm): [{geo.sVoxel[0]:.2f}, {geo.sVoxel[1]:.2f}, {geo.sVoxel[2]:.2f}]\n")
        file.write(f"Offset from Origin (mm): [{geo.offOrigin[0]:.2f}, {geo.offOrigin[1]:.2f}, {geo.offOrigin[2]:.2f}]\n")
        file.write(f"Offset of Detector (mm): [{geo.offDetector[0]:.2f}, {geo.offDetector[1]:.2f}]\n")
        file.write(f"# projection angles: {len(ang)}\n")
        file.write('\n')


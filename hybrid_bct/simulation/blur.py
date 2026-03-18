from scipy.interpolate import interp1d

def mtf_blur(prjstack_calcs, dexel_mm, f_MTF, mtf_MTF):
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

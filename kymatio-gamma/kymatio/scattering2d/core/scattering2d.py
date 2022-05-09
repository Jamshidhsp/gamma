# Authors: Georgios Exarchakis
# Scientific Ancestry: Muawiz Chaudhary, Edouard Oyallon, Laurent Sifre, Joan Bruna,

def scattering2d(x, pad, unpad, backend, J, L, phi, psi, max_order,
        out_type='array'):
 
    subsample_fourier = backend.subsample_fourier
    modulus = backend.modulus
    fft = backend.fft
    log = backend.log
    cdgmm = backend.cdgmm
    concatenate = backend.concatenate
    concat_gamma = backend.concat_gamma
    concat_gamma2 = backend.concat_gamma2

    # Define lists for output.
    out_S_0, out_S_1, out_S_2 = [], [], []
    

    N = 2**(2*J)
    
    
    U_r_c = pad(x)
    U_r = modulus(U_r_c)


    # Pareto components
    U_0_r_l = (log((U_r+1e-32) / 1e-32))

    U_0_c = fft(U_r, 'C2C')
    
    # Fourrier Transform of Pareto stuff
    U_0_c_lf = fft(U_0_r_l, 'C2C')

    # First low pass filter
    U_1_c = cdgmm(U_0_c, phi[0])
    U_1_c = subsample_fourier(U_1_c, k=2 ** J)

    S_0 = fft(U_1_c, 'C2R', inverse=True)
    S_0 = unpad(S_0)


    # First low pass filter (Pareto params)
    S_0_c_lf = cdgmm(U_0_c_lf, phi[0])
    S_0_c_lf = subsample_fourier(S_0_c_lf, k=2 ** J )

    S_0_r_l = fft(S_0_c_lf, 'C2R', inverse=True)
    S_0_r_l = unpad(S_0_r_l)




    out_S_0.append({'coef': S_0,
                    'coef_pareto': N/S_0_r_l,
                    'j': (),
                    'theta': ()})

    for n1 in range(len(psi)):
        j1 = psi[n1]['j']
        theta1 = psi[n1]['theta']

        U_1_c = cdgmm(U_0_c, psi[n1][0])
        if j1 > 0:
            U_1_c = subsample_fourier(U_1_c, k=2 ** j1)
        U_1_c = fft(U_1_c, 'C2C', inverse=True)
        U_1_c = modulus(U_1_c)
        
        # Pareto components
        U_1_c_l = (log((U_1_c+1e-32) / 1e-32))

        U_1_c = fft(U_1_c, 'C2C')
        
        # Fourrier Transform of Pareto stuff
        U_1_c_lf = fft(U_1_c_l, 'C2C')

        # Second low pass filter
        S_1_c = cdgmm(U_1_c, phi[j1])
        S_1_c = subsample_fourier(S_1_c, k=2 ** (J - j1))

        S_1_r = fft(S_1_c, 'C2R', inverse=True)
        S_1_r = unpad(S_1_r)


        # Second low pass filter (Pareto params)
        S_1_c_lf = cdgmm(U_1_c_lf, phi[j1])
        S_1_c_lf = subsample_fourier(S_1_c_lf, k=2 ** (J - j1))

        S_1_r_l = fft(S_1_c_lf, 'C2R', inverse=True)
        S_1_r_l = unpad(S_1_r_l)

        out_S_1.append({'coef': S_1_r,
                        'coef_pareto': N/S_1_r_l,
                        'j': (j1,),
                        'theta': (theta1,)})

        if max_order < 2:
            continue
        for n2 in range(len(psi)):
            j2 = psi[n2]['j']
            theta2 = psi[n2]['theta']

            if j2 <= j1:
                continue

            U_2_c = cdgmm(U_1_c, psi[n2][j1])
            U_2_c = subsample_fourier(U_2_c, k=2 ** (j2 - j1))
            U_2_c = fft(U_2_c, 'C2C', inverse=True)
            U_2_c = modulus(U_2_c)

            # Gamma stuff
            U_2_c_l = (log((U_2_c+1e-32) / 1e-32))

            U_2_c = fft(U_2_c, 'C2C')

            # Fourrier Transform of Gamma stuff
            U_2_c_lf = fft(U_2_c_l, 'C2C')
            # U_2_c_pf = fft(U_2_c_p, 'C2C')

            # Third low pass filter
            S_2_c = cdgmm(U_2_c, phi[j2])
            S_2_c = subsample_fourier(S_2_c, k=2 ** (J - j2))

            S_2_r = fft(S_2_c, 'C2R', inverse=True)
            S_2_r = unpad(S_2_r)

            # Third low pass filter (GAMMA params)
            S_2_c_lf = cdgmm(U_2_c_lf, phi[j2])
            S_2_c_lf = subsample_fourier(S_2_c_lf, k=2 ** (J - j2))

            S_2_r_l = fft(S_2_c_lf, 'C2R', inverse=True)
            S_2_r_l = unpad(S_2_r_l)

            out_S_2.append({'coef': S_2_r,
                            'coef_pareto': N/S_2_r_l,
                            'j': (j1, j2),
                            'theta': (theta1, theta2)})

    out_S = []
    out_S.extend(out_S_0)
    out_S.extend(out_S_1)
    out_S.extend(out_S_2)

    if out_type == 'array':
        
#        out_S = concat_gamma(concatenate([x['coef'] for x in out_S]), concatenate([x['coef_shape'] for x in out_S]), concatenate([x['coef_scale'] for x in out_S]))
        # out_S = concat_gamma2(concatenate([x['coef_shape'] for x in out_S]), concatenate([x['coef_scale'] for x in out_S]))
#        out_S = concat_gamma2(concatenate([x['coef_shape'] for x in out_S]), concatenate([x['coef_rate'] for x in out_S]))

       out_S = concatenate([x['coef_pareto'] for x in out_S])
        
    return out_S
#    return [U_r, U_0_r_l]


__all__ = ['scattering2d']

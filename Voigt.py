import scipy as np

def voigt2TensorStrain(strain):
    if len(strain) == 3:
        eps = np.zeros((2, 2))
        eps[0, 0] = strain[0]
        eps[1, 1] = strain[1]
        eps[0, 1] = eps[1, 0] = strain[2]/2
        return eps
    elif len(strain) == 6:
        # e_xx, e_yy, e_zz, g_xy, g_yz, g_xz
        eps = np.zeros((3, 3))
        eps[0, 0] = strain[0]
        eps[1, 1] = strain[1]
        eps[2, 2] = strain[2]
        eps[0, 1] = eps[1, 0] = strain[3]/2
        eps[1, 2] = eps[2, 1] = strain[4]/2
        eps[0, 2] = eps[2, 0] = strain[5]/2
        return eps
    else:
        return ValueError("Strain vector has to be of length 3 or 6 !")


def voigt2TensorStress(stress):
    if len(stress) == 3:
        eps = np.zeros((2, 2))
        eps[0, 0] = stress[0]
        eps[1, 1] = stress[1]
        eps[0, 1] = eps[1, 0] = stress[2]
        return eps
    elif len(stress) == 6:
        # e_xx, e_yy, e_zz, g_xy, g_yz, g_xz
        eps = np.zeros((3, 3))
        eps[0, 0] = stress[0]
        eps[1, 1] = stress[1]
        eps[2, 2] = stress[2]
        eps[0, 1] = eps[1, 0] = stress[3]
        eps[1, 2] = eps[2, 1] = stress[4]
        eps[0, 2] = eps[2, 0] = stress[5]
        return eps
    else:
        return ValueError("Stress vector has to be of length 3 or 6 !")

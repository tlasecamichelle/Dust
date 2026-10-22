import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, FloatSlider
from scipy.constants import h, c, k

#bulk densities [g cm^-3] (from Bruce Draine)
rho_gra = 2.24
rho_sil = 3.5 

#figuring out pulling the dust info:
radii_list = []
Q_abs_list = []
Q_sca_list = []

def parse_bruce_data(filepath):

    with open(filepath, 'r') as f:
        lines = f.readlines()

    radii_list = [] #empty list to store grain radius size for each data block

    wavelengths = None #placeholder, will be wavelength grid array from blocks

    Q_abs_matrix = [] #creates list holding absoprtion efficiency values Qabs

    Q_sca_matrix = [] #create list to hold scat efficiency values Qsca


    i=0 #index at 1st line of file

    while i< len(lines): #while loop, read every line through file

        line = lines[i].strip() #grabs line and cleans away useless characters


        #checking if new radius block

        if 'radius(micron)' in line: #scanning for "radius(micron)"

            parts = line.split() #extracting radius value from start of line

            current_radius = float(parts[0]) #grabbing first piece of list, converting to number(grain size)

            radii_list.append(current_radius) #adds grain radius to list


            #skip past the column headers
            i +=4


            block_waves= [] #create empty list for specific grain size

            blocks_qabs = [] #Q abs absorption coefficiency to calc mass absorp coeff kappanu, for specific grain size

            blocks_qsca = [] # Q sca scattering efficiency, for THIS grain block


            #read rows until hit end of text

            while i< len(lines): #sub-loop to read rows inside current block

                cols = lines[i].split() #split data into list of cols

                #if line doesn't have at least 4 numerical col, already hit next block
                if len(cols)<4:

                    break

                try: #sub-loop safety to convert col to floats (number), can probs ignore this was a suggested add by ai 

                    w = float(cols[0])

                    qabs = float(cols[1])

                    qsca = float(cols[2])


                    block_waves.append(w)

                    blocks_qabs.append(qabs)

                    blocks_qsca.append(qsca)

                except ValueError:

                    break

                i +=1



#save waveelength grid once for every radius block, storing and returning

            if wavelengths is None: #every grain size block uses same wavelength grid, get the array once and move on

                wavelengths = np.array(block_waves)


            Q_abs_matrix.append(blocks_qabs) #append completed list of Qabs for this grain size into matrix

            Q_sca_matrix.append(blocks_qsca) #same as above, but for Qsca


        else:

            i +=1 #if a line isnt a radius block, increments by a step and moves on 

    return wavelengths, np.array(radii_list), np.array(Q_abs_matrix), np.array(Q_sca_matrix) #once scan file, extract numpy arrays for 4 outputs


wavelengths_sil, radii_sil, Q_abs_sil, Q_sca_sil = parse_bruce_data('Sil_81')

#print(f"Silicate: Radii count: {len(radii_sil)}, Wavelengths count: {len(wavelengths_sil)}")

wavelengths_gra, radii_gra, Q_abs_gra, Q_sca_gra = parse_bruce_data('Gra_81')

#print(f"Graphite: Radii count: {len(radii_gra)}, Wavelengths count: {len(wavelengths_gra)}")

#radii of use for graphite
radii_gra_log10 = np.log10(radii_gra)
log10_mask = np.isclose(radii_gra_log10, np.round(radii_gra_log10))
lam_gra = (radii_gra[log10_mask])

#radii of use for silicate
radii_sil_log10 = np.log10(radii_sil)
log10_mask = np.isclose(radii_sil_log10, np.round(radii_sil_log10))
lam_sil = (radii_sil[log10_mask])

#'Q_abs for 1e-3: {Q_abs_gra[0]}')
#'Q_abs for 1e-2: {Q_abs_gra[20]}')
#'Q_abs for 1e-1: {Q_abs_gra[40]}')
#'Q_abs for 1e-0: {Q_abs_gra[60]}')
#'Q_abs for 1e-1: {Q_abs_gra[80]}')

#the indicies of radii_gra or radii_sil that return
Q_abs_ind = [0, 20, 40, 60, 80]

def radii_index(radii, a_um):
    idx = np.argmin(np.abs(np.log10(radii)-np.log10(a_um))) 
    '''#returns the index where difference of radii_gra/sil value and user input a_um value is smallest 
    (applies np.log10 (because the radii list is logspaced) to all of whatever radii array i specify 
    and subtract the np.log10(a_um) across all values of array)'''
    return idx

r'''
plt.figure(1)
for i in Q_abs_ind:
    plt.loglog(wavelengths_gra, Q_abs_gra[i], label = f'r={radii_gra[i]}' + r' $\mu$m')
plt.title('Q_abs of Various Graphite Dust Grain Sizes Across a Spectra')
plt.ylabel('Q_abs')
plt.xlabel(r'$\lambda$ in $\mu$m')
plt.legend()

plt.figure(2)
for i in Q_abs_ind:
    plt.loglog(wavelengths_sil, Q_abs_sil[i], label = f'r={radii_sil[i]}' + r' $\mu$m')
plt.title('Q_abs of Various Silicate Dust Grain Sizes Across a Spectra')
plt.ylabel('Q_abs')
plt.xlabel(r'$\lambda$ in $\mu$m')
plt.legend()

plt.show()
'''
#planck fxns
def B_nu(nu, T):
    L = ((2*h*nu**3)/c**2)/(np.expm1((h*nu)/(k*T)))
    return L

def B_lam(lam, T):
    L = (2*h*c**2/lam**5)/(np.expm1((h*c)/(lam*k*T)))
    return L

lam = np.logspace(-8, -2, 500)
nu = c/lam

def k_lam(Q_abs_row, a_um, rho):
    a_cm = a_um*1e-4
    k = (3*(Q_abs_row))/(4*rho*a_cm) * 0.1
    return k #rho is cm^2/g outputs to m^2/kg

def F_lam(wavelength_um, Q_abs_row, a_um, T, rho, M_solmass, d_Mpc):
    M = M_solmass*1.989e30
    d = d_Mpc*1000000 * 3.086e16
    wavelength_m = wavelength_um*1e-6
    F = (M*B_lam(wavelength_m, T)*k_lam(Q_abs_row, a_um, rho))/d**2
    return F
#M is logically given as solar masses, d in parsecs or some other large measurement

DUST = {
    'carbon':   dict(wav=wavelengths_gra, Q=Q_abs_gra, radii=radii_gra, rho=rho_gra),
    'silicate': dict(wav=wavelengths_sil, Q=Q_abs_sil, radii=radii_sil, rho=rho_sil),
}

def F_lam_component(composition, a_um, T, M_solmass, d_Mpc):
    tab = DUST[composition]
    j = radii_index(tab['radii'], a_um)
    return F_lam(tab['wav'], tab['Q'][j], tab['radii'][j], T, tab['rho'], M_solmass, d_Mpc)

def F_lam_multi(components, d_Mpc):
    return sum(F_lam_component(composition, a_um, T, M_solmass, d_Mpc) for composition, a_um, T, M_solmass in components)

#Assume NGC 2906 > 30 MPc
d_1 = 30
#1 MPc > 1,000,000 Pc > 3.086e16 m
#1 Solar Mass > 1.989e30 kg
M_1 = 1
nl = '\n'

#composition, a_um, T, M_solmass, d_Mpc
components_2017eaw_epoch1 = [('silicate', 0.1, 157, 6.2e-4),
                            ('silicate', 0.1, 2503, 4.6e-8)] 

components_2017eaw_epoch2 = [('silicate', 0.1, 157, 5.5e-4),
                            ('silicate', 0.1, 2292, 3.3e-8)] 

d_2017eaw_Mpc = 7.12

#Outputs of dustysn for 2017eaw epoch 1 into my own modified bb graph:
def semilogplot_F_lam_multi(components, d_Mpc, label):
    wav = DUST['carbon']['wav']
    F = F_lam_multi(components, d_Mpc)
    plt.semilogy(wav, F, label=label)

def loglogplot_F_lam_multi(components, d_Mpc, label):
    wav = DUST['carbon']['wav']
    F = F_lam_multi(components, d_Mpc)
    plt.loglog(wav, F, label=label)

#example graphs that vary across one variable per:
def graphite_radii():
    #lam_F of graphite across multiple radii > warm dust only (400K)
    plt.figure(1)
    plt.subplot(1, 3, 2)
    for i in Q_abs_ind:
        plt.loglog(wavelengths_gra, F_lam(wavelengths_gra, Q_abs_gra[i], radii_gra[i], 400, rho_gra, M_1, d_1), label = f'r={radii_gra[i]}' + r' $\mu$m')
    plt.title('Blackbody: 1 Solar Mass of Warm \nGraphite Dust (400K) of Varying Grain Size')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-5)
    plt.xlim(1e-1)
    plt.legend()

    #lam_F of graphite across multiple radii > hot dust only (800K)
    #plt.figure(2)
    plt.subplot(1, 3, 1)
    for i in Q_abs_ind:
        plt.loglog(wavelengths_gra, F_lam(wavelengths_gra, Q_abs_gra[i], radii_gra[i], 800, rho_gra, M_1, d_1), label = f'r={radii_gra[i]}' + r' $\mu$m')
    plt.title('Blackbody: 1 Solar Mass of Hot \nGraphite Dust (800K) of Varying Grain Size')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-4)
    plt.xlim(1e-1)
    plt.legend()

    #lam_F of graphite across multiple radii > cold dust (100K)
    #plt.figure(3)
    plt.subplot(1, 3, 3)
    for i in Q_abs_ind:
        plt.loglog(wavelengths_gra, F_lam(wavelengths_gra, Q_abs_gra[i], radii_gra[i], 100, rho_gra, M_1, d_1), label = f'r={radii_gra[i]}' + r' $\mu$m')
    plt.title('Blackbody: 1 Solar Mass of Cold \nGraphite Dust (100K) of Varying Grain Size')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-8)
    plt.xlim(1e0)
    plt.legend()
    plt.show()

def silicon_radii():
    #lam_F of silicon across multiple radii > warm dust (400K)
    plt.figure(2)
    plt.subplot(1, 3, 2)
    for i in Q_abs_ind:
        plt.loglog(wavelengths_sil, F_lam(wavelengths_sil, Q_abs_sil[i], radii_sil[i], 400, rho_sil, M_1, d_1), label = f'r={radii_sil[i]}' + r' $\mu$m')
    plt.title('Blackbody: 1 Solar Mass of Warm \nSilicon Dust (400K) of Varying Grain Size')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-6)
    plt.xlim(1e-1)
    plt.legend()

    #lam_F of silicon across multiple radii > hot dust (800K)
    plt.subplot(1, 3, 1)
    for i in Q_abs_ind:
        plt.loglog(wavelengths_sil, F_lam(wavelengths_sil, Q_abs_sil[i], radii_sil[i], 800, rho_sil, M_1, d_1), label = f'r={radii_sil[i]}' + r' $\mu$m')
    plt.title('Blackbody: 1 Solar Mass of Hot \nSilicon Dust (800K) of Varying Grain Size')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-5)
    plt.xlim(1e-1)
    plt.legend()

    #lam_F of silicon across multiple radii > cold dust (100K)
    plt.subplot(1, 3, 3)
    for i in Q_abs_ind:
        plt.loglog(wavelengths_sil, F_lam(wavelengths_sil, Q_abs_sil[i], radii_sil[i], 100, rho_sil, M_1, d_1), label = f'r={radii_sil[i]}' + r' $\mu$m')
    plt.title('Blackbody: 1 Solar Mass of Cold \nSilicon Dust (100K) of Varying Grain Size')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-9)
    plt.xlim(1e0)
    plt.legend()
    plt.show()

mass_ind = [1e-4, 1e-3, 1e-2, 1e-1, 1e0]
def graphite_dustmass():
    #lam_F graphite, but mass is changing (hot, warm, cold), assume 1e-2 micron grains > index 20
    plt.figure(3)
    plt.subplot(1, 3, 1) #hot dust
    for i in mass_ind:
        plt.loglog(wavelengths_gra, F_lam(wavelengths_gra, Q_abs_gra[20], radii_gra[20], 800, rho_gra, i, d_1), label = f'Mass={[i]}'+ r' $M_{\odot}$')
    plt.title(r'''Blackbody: 0.01 $\mu$m Hot Graphite 
    Dust (800K) of Varying $M_{\odot}$''')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-5)
    plt.xlim(1e-2)
    plt.legend()

    plt.subplot(1, 3, 2) #warm dust
    for i in mass_ind:
        plt.loglog(wavelengths_gra, F_lam(wavelengths_gra, Q_abs_gra[20], radii_gra[20], 400, rho_gra, i, d_1), label = f'Mass={[i]}'+ r' $M_{\odot}$')
    plt.title(r'''Blackbody: 0.01 $\mu$m Warm Graphite 
    Dust (400K) of Varying $M_{\odot}$''')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-6)
    plt.xlim(1e-1)
    plt.legend()

    plt.subplot(1, 3, 3) #cold dust
    for i in mass_ind:
        plt.loglog(wavelengths_gra, F_lam(wavelengths_gra, Q_abs_gra[20], radii_gra[20], 100, rho_gra, i, d_1), label = f'Mass={[i]}'+ r' $M_{\odot}$')
    plt.title(r'''Blackbody: 0.01 $\mu$m Cold Graphite 
    Dust (100K) of Varying $M_{\odot}$''')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-17, 1e-9)
    plt.xlim(1e0)
    plt.legend()
    plt.show()

def silicate_dustmass():
    #lam_F silicate, but mass is changing (hot, warm, cold), assume 1e-2 micron grains > index 20
    plt.figure(4)
    plt.subplot(1, 3, 1) #hot dust
    for i in mass_ind:
        plt.loglog(wavelengths_sil, F_lam(wavelengths_sil, Q_abs_sil[20], radii_sil[20], 800, rho_sil, i, d_1), label = f'Mass={[i]}'+ r' $M_{\odot}$')
    plt.title(r'''Blackbody: 0.01 $\mu$m Hot Silicate 
    Dust (800K) of Varying $M_{\odot}$''')
    plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-5)
    plt.xlim(1e-2)
    plt.legend()

    plt.subplot(1, 3, 2) #warm dust
    for i in mass_ind:
        plt.loglog(wavelengths_sil, F_lam(wavelengths_sil, Q_abs_sil[20], radii_sil[20], 400, rho_sil, i, d_1), label = f'Mass={[i]}'+ r' $M_{\odot}$')
    plt.title(r'''Blackbody: 0.01 $\mu$m Warm Silicate 
    Dust (400K) of Varying $M_{\odot}$''')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-16, 1e-6)
    plt.xlim(1e-1)
    plt.legend()

    plt.subplot(1, 3, 3) #cold dust
    for i in mass_ind:
        plt.loglog(wavelengths_sil, F_lam(wavelengths_sil, Q_abs_sil[20], radii_sil[20], 100, rho_sil, i, d_1), label = f'Mass={[i]}'+ r' $M_{\odot}$')
    plt.title(r'''Blackbody: 0.01 $\mu$m Cold Silicate 
    Dust (100K) of Varying $M_{\odot}$''')
    plt.xlabel(r"$\lambda$ in $\mu$m")
    plt.ylim(1e-17, 1e-9)
    plt.xlim(1e0)
    plt.legend()
    plt.show()

plt.figure(1)
plt.subplot(1,2,1)
semilogplot_F_lam_multi(components_2017eaw_epoch1, d_2017eaw_Mpc, 'Epoch 1')
semilogplot_F_lam_multi(components_2017eaw_epoch2, d_2017eaw_Mpc, 'Epoch 2')
plt.title('Mid-IR Focus')
plt.ylabel(r'$F_\lambda$ (W m$^{-2}$ m$^{-1}$ )')
plt.xlim(5,30)
plt.ylim(1e-13, 0.5e-11)
plt.legend()
plt.subplot(1,2,2)
loglogplot_F_lam_multi(components_2017eaw_epoch1, d_2017eaw_Mpc, 'Epoch 1')
loglogplot_F_lam_multi(components_2017eaw_epoch2, d_2017eaw_Mpc, 'Epoch 2')
plt.title('Optical + IR')
plt.xlim(0.400,50)
plt.ylim(1e-16,1e-9)
plt.axvline(0.7, color='k', ls='--', label='Optical>IR')
plt.suptitle('Dustysn 2017eaw fit outputs plugged into modified blackbody')
plt.xlabel(r"$\lambda$ in $\mu$m")
plt.legend()
plt.show()
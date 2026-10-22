import numpy as np

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

print(f"Silicate: Radii count: {len(radii_sil)}, Wavelengths count: {len(wavelengths_sil)}")

wavelengths_gra, radii_gra, Q_abs_gra, Q_sca_gra = parse_bruce_data('Gra_81')

print(f"Graphite: Radii count: {len(radii_gra)}, Wavelengths count: {len(wavelengths_gra)}")

DUST = {
    'graphite':   dict(wav=wavelengths_gra, Q=Q_abs_gra, radii=radii_gra, rho=rho_gra),
    'silicate': dict(wav=wavelengths_sil, Q=Q_abs_sil, radii=radii_sil, rho=rho_sil),
}

#fxn to turn the Q_abs array into k_abs array
def k_lam_cgs(Q_abs_row, a_um, rho):
    a_cm = a_um*1e-4
    k = (3*(Q_abs_row))/(4*rho*a_cm)
    return k #rho is cm^2/g and keeps outputs as such to follow dustysn

def bruce_to_dustysn(composition):
    tab = DUST[composition]
    Q_abs_ind = [0, 20, 40, 60, 80]

    wavelengths = wavelengths_gra[::-1] #reverses it be ascending cuz bruce draine files are descending
    header = ' '.join(['wavelength(mic)'] + [f'abs_k(a={tab['radii'][i]}mic)' for i in Q_abs_ind])
    cols = [wavelengths] + [k_lam_cgs(tab['Q'][i], tab['radii'][i], tab['rho'])[::-1] for i in Q_abs_ind] #outputs 6 arrays
    #print(f'Heres header: {cols}')
    #print(f'First list {totalcols[1]}')
    stack = np.column_stack(cols)
    print(f'before header {header}')
    dustysn_format_bruce_data = stack


    np.savetxt(f'k_{composition}_BruceDraine', dustysn_format_bruce_data, header=header)

bruce_to_dustysn('silicate')
bruce_to_dustysn('graphite')
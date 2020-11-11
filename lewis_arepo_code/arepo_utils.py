#!/usr/bin/python
import os
import numpy as np

k_B =1.3806581e-16
mp = 1.6726575e-24

#
# This is our arepo reader for snapshot type 2
def aread(filename):
        #
        # Define the function to read the extra 'type 2' snapshot info
        def read_tag_and_nextblock(fobj):
                skip = np.fromfile(fobj,dtype=np.int32,count=1)
                #
                # tag first!
                stringarr = np.fromfile(fobj,dtype=np.int8,count=4)
                if (np.sum(stringarr) > 0 ):
                    tag = ''.join([chr(item) for item in stringarr])
                else:
                    tag = 0
                #
                # now position of next block
                dumread = np.fromfile(fobj,dtype=np.int32,count=1)
                if (len(dumread) >0):
                    nextblock = dumread[0] 
                else:
                    nextblock = 0
                #
                # final skip to comlete the read
                skip = np.fromfile(fobj,dtype=np.int32,count=1)
                return tag, nextblock

        #
        # define a class to hold all our arepo data. The class (like a struct) is
        # what we return from the function
        class arepo_data_struct:
            pass

        print("Reading file:", filename)

        bytesleft= 256 - 6*4 - 6*8 - 2*8 - 2*4 -6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 3*4 - 4

        with open(filename,"rb") as file:
                #
                # start with the header, but read tag and record length in first
                tag, nextblock = read_tag_and_nextblock(file)
                print("Found array with tag ",tag, "and bytes to next block", nextblock)
                #
                # the first skip (the opening f77 record length)
                skip = np.fromfile(file,dtype=np.int32,count=1)
                #
                # now read the actual header 
                npart = np.fromfile(file, dtype=np.int32, count=6)
                massarr = np.fromfile(file, dtype=np.double,count=6)
                time = np.fromfile(file,dtype=np.double,count=1)[0]
                redshift = np.fromfile(file,dtype=np.double,count=1)[0]
                flag_sfr = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_feedback = np.fromfile(file,dtype=np.int32,count=1)[0]
                npartTotal = np.fromfile(file,dtype=np.int32,count=6)
                flag_cooling = np.fromfile(file,dtype=np.int32,count=1)[0]
                num_files = np.fromfile(file,dtype=np.int32,count=1)[0]
                boxsize = np.fromfile(file,dtype=np.double,count=1)[0]
                cos1 = np.fromfile(file,dtype=np.double,count=1)[0]
                cos2 = np.fromfile(file,dtype=np.double,count=1)[0]
                hubble_param = np.fromfile(file,dtype=np.double,count=1)[0]
                flag_stellarage = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_metals = np.fromfile(file,dtype=np.int32,count=1)[0]
                npartHighword = np.fromfile(file,dtype=np.int32,count=6)
                flag_entropy = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_dp = np.fromfile(file,dtype=np.int32,count=1)[0]
                flag_1pt = np.fromfile(file,dtype=np.int32,count=1)[0]
                scalefactor = np.fromfile(file,dtype=np.float32,count=1)[0]
                print("flag_dp, flag_1pt, scalefactor", flag_dp, flag_1pt, scalefactor)
                pad = np.fromfile(file,dtype=np.int32,count=(bytesleft//4)) 
                # the second skip (the closeing f77 record length)
                skip = np.fromfile(file,dtype=np.int32,count=1)

                print("npart array:", npart)
                print("mass array:", massarr)
                print("time in codeunits:", time)
                print("Total npart:", npartTotal)
                print("Header finished")

                N = int(sum(npart)) - npart[4] # type 4 is reserved for TRACER_MC
                ngas = npart[0]
                nsink = npart[5]

                #
                # Initialise our struct that holds the data                
                a = arepo_data_struct()
                a.npart = npart
                a.ngas = ngas
                a.nsink = nsink
                a.N = N
                a.time = time

                #
                # Set the units here
                a.unit_leng_cm = 1.0e+17
                a.unit_mass_g = 1.991e33
                a.unit_time_s = 2.7436898e12

                #
                # flags to indicate whether a property has been read
                igot_pot = 0
                igot_accel = 0
                igot_dt = 0
                igot_bmag = 0
                igot_chem = 0
                igot_soft = 0

                #
                # start with the header, but read tag and length to next block first
                tag, nextblock = read_tag_and_nextblock(file)
                print("Found array with tag ",tag, "and bytes to next block", nextblock)
                while (nextblock > 0):
                        if(tag=="POS "):
                                print("Reading positions")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                buf = np.fromfile(file,dtype=np.double,count=3*N)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                buf = np.reshape(buf,(-1,3))
                                a.x = buf[:, 0]
                                a.y = buf[:, 1]
                                a.z = buf[:, 2]
                        elif(tag=="VEL "):
                                print("Reading velocities")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                buf = np.fromfile(file,dtype=np.double,count=3*N)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                buf = np.reshape(buf,(-1,3))
                                a.vx = buf[:, 0]
                                a.vy = buf[:, 1]
                                a.vz = buf[:, 2]
                        elif(tag=="ID  "):
                                print("Reading particle IDs")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.partid = np.fromfile(file,dtype=np.int32,count=N)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="MASS"):
                                print("Reading particle masses")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.mass = np.fromfile(file,dtype=np.double,count=N)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="U   "):
                                print("Reading u")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.u = np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="RHO "):
                                print("Reading densities")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.rho = np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="POT "):
                                print("Reading potentials")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.potential = np.fromfile(file,dtype=np.double,count=N)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                igot_pot = 1
                        elif(tag=="DIVV"):
                                print("Reading velocity divergence")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.divv = np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="ACCE"):
                                print("Reading accelerations")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                buf = np.fromfile(file,dtype=np.double,count=N*3)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.accel = np.reshape(buf,(-1,3))
                                igot_accel = 1
                        elif(tag=="DUST"):
                                print("Reading dust temperatures")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.tdust = np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="TSTP"):
                                print("Reading timesteps")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.dt = np.fromfile(file,dtype=np.double,count=N)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                igot_dt = 1
                        elif(tag=="BFLD"):
                                print("Reading magnetic field")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                buf = np.fromfile(file,dtype=np.double,count=3*ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.bfield = np.reshape(buf,(-1,3))
                                igot_bmag = 1
                        elif(tag=="DIVB"):
                                print("Reading magnetic field divergence")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.divb = np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="SOFT"):
                                print("Reading softening")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.softening = np.fromfile(file,dtype=np.double,count=N)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                igot_soft = 1
                        elif(tag=="CHEM"):
                                print("Reading chemistry")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                num_species = (nextblock - 8) // ngas // 8
                                print('There appear to be ', num_species, ' chemical species in the network')
                                buf = np.fromfile(file,dtype=np.double,count=num_species*ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.chem = np.reshape(buf,(-1,num_species))
                                igot_chem = 1
                        elif(tag=="ROTV"):
                                print("Reading the magnitude of the curl")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.curlvel = np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="VDIS"):
                                print("Reading velocity dispersion")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.veldisp = np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=="PEAK"):
                                print("Reading the peaks of the potential")
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.peak = np.fromfile(file,dtype=np.int32,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        elif(tag=='FACA'):
                                print('Reading face angles')
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                                a.angle=np.fromfile(file,dtype=np.double,count=ngas)
                                skip = np.fromfile(file,dtype=np.int32,count=1)
                        else:
                                print("Skipping through property",tag," with record length", nextblock-8)
                                dummy = np.fromfile(file,dtype=np.int32,count=(nextblock//4))

                        #
                        # read the new tag and length to next block
                        tag, nextblock = read_tag_and_nextblock(file)
                        print("Found array with tag ",tag, "and bytes to next block", nextblock)


        #
        # if we have sinks present, let them go in separate variables
        if(a.nsink > 0):
            #
            # first, copy the sink properites to their own arrays
            print("Sinks read. Making sink arrays.")
            n_not_sink=int(sum(npart)-npart[5])
            a.idsink = np.linspace(0,nsink-1,nsink,dtype='int32') + n_not_sink #ngas
            a.sinkx = a.x[a.idsink]
            a.sinky = a.y[a.idsink]
            a.sinkz = a.z[a.idsink]
            a.sinkvx = a.vx[a.idsink]
            a.sinkvy = a.vy[a.idsink]
            a.sinkvz = a.vz[a.idsink]
            a.sinkmass = a.mass[a.idsink]
            a.sinkid = a.partid[a.idsink]

        #
        # now remove the sinks from the other arrays (so pos, vel, mass, id etc have
        # length ngas, not N)
        i_not_gas = [1, 2, 3, 5] # type 4 is reservered for TRACER_MC
        if(sum(npart[i_not_gas]) > 0):
            igas = np.linspace(0,ngas-1,ngas,dtype='int32')
            a.x = a.x[igas]
            a.y = a.y[igas]
            a.z = a.z[igas]
            a.vx = a.vx[igas]
            a.vy = a.vy[igas]
            a.vz = a.vz[igas]
            a.partid = a.partid[igas]
            a.mass = a.mass[igas]
            if (igot_pot == 1): 
                a.potential = a.potential[igas]
            if (igot_accel == 1):
                a.accel = a.accel[igas, :]
            if (igot_dt == 1):
                a.dt = a.dt[igas]
            if (igot_soft == 1):
                a.softening= a.softening[igas]


        #
        # if we have chemistry, then we need to provide the temperature
        if (igot_chem > 0):
            print('Creating an array with T [K] from specific energies')
            ABHE = 0.1
            uenergy = 2.64481e+42
            ulength = 1e17
            udensity = 1.991e-18 
            yn = a.rho*udensity / ((1.0 + 4.0 * ABHE) * mp)
            energy = a.u * a.rho * uenergy / ulength**3
            yntot = (1.0 + ABHE - a.chem[:, 0] + a.chem[:, 1]) * yn
            a.temp = 2.0 * energy / (3.0 * yntot * k_B)

        #
        # get radius from densest point, com, or first sink
        if (a.nsink > 0):
           a.rad = np.sqrt((a.x - a.sinkx[0])**2 + (a.y - a.sinky[0])**2 + (a.z - a.sinkz[0])**2)


        #
        #
        if(igot_bmag > 0):
           a.l = (a.mass / a.rho)**(1./3.)
           a.bmag = np.sqrt(a.bfield[:,0]**2 + a.bfield[:,1]**2 + a.bfield[:,2]**2)   
 
        print("Finished reading file:", filename,"\n")
        
        return a  # we simply return the struct a.--- that holds all the data! 


#
# The function for making quick images - based on the 2d histogram
def arepoimage(x, y, weight):
        import matplotlib.pyplot as plt
        import numpy as np
        #
        # output some potentially useful info
        print('xrange', min(x), max(x))
        print('yrange', min(y), max(y))
        #
        # make 2 histograms, one with the weights, one without
        hist_weighted, xb, yb = np.histogram2d(y, x, weights=weight, bins=(500, 500))
        hist_numbers, xb, yb = np.histogram2d(y, x, bins=(500, 500))
        #
        # divide the weighted by the normal histogram to get an image of mean values
        hist_final = hist_weighted / hist_numbers
        hist_final = np.ma.masked_where(hist_numbers < 1, hist_final)
        #
        # do we need to log the image?
        ip = np.where(hist_numbers > 0)
        print('len', len(ip[0]))
        print(ip)
        max_image = np.max(hist_final[ip])
        min_image = np.min(hist_final[ip])
        if ( (max_image/min_image) > 50 ):
            hist_final = np.nan_to_num(np.log10(hist_final))
            print('The image range is > 50 so it will be logged')
        #
        # make the plot
        plt.clf()
        image = plt.imshow(hist_final, aspect='auto', cmap='plasma', origin='lower', extent=[yb[0],yb[-1],xb[0],xb[-1]])
        plt.colorbar()
        plt.show()



def plot_imf(mass):
        import matplotlib.pyplot as plt
        import numpy as np

        ip = np.where(mass > 0)
        log_mass = np.log10(mass[ip])

        plot_line_masses = np.zeros(3)
        plot_line_masses[0] = -1.0
        plot_line_masses[1] = 0.0
        plot_line_masses[2] = 1.0
        # make a straight line to overplot
        num_exp_eqmass  = 10.0**(-1.0*plot_line_masses +  1)
        num_exp_salpeter = 10.0**(-1.35*plot_line_masses + 1)

        # make a plot!
        n, bins, patches = plt.hist(log_mass, 20, log=True, alpha=0.5)

        plt.plot(plot_line_masses, num_exp_eqmass)
        plt.plot(plot_line_masses, num_exp_salpeter)
        plt.show()

import inputclass
import os 
from convolution import * 

def runManyCowan(input: inputclass.Input):
    
    '''
    Runs the cowan code many times. generates an input deck.
    '''
    
    if input.runCowan:
        print('Commencing Cowan-Configuration Average Run.')
        
        for i in range(0,len(input.shell_list)):
            print('Ionizing shell ',input.shell_list[i])
            run_calc(i,input)

    if input.runRates:
        print('Convolving cross sections into Maxwellian rates.')
        print(' ADF-style format in "ratesFile.dat".')
        print(' Another format in "plotRates.dat".')

        produceRates(input)
    
    return 0 

def generateInput(ionizedIndex,
                  input:inputclass.Input):

    base_string = '   {:2}   {:2}  {:2}+{:1} init        '.format(input.atomic_num,
                                                                  input.ionstage,
                                                                  input.symbol,
                                                                  input.ionstage-1)
    ion_string  = '   {:2}   {:2}  {:2}+{:1} ion         '.format(input.atomic_num,
                                                                  input.ionstage+1,
                                                                  input.symbol,
                                                                  input.ionstage)
    
    #making the configs.
    for (i,shell) in enumerate(input.shell_list):
        base_string += shell + '{:<2} '.format(input.init_occup[i])
        if (i ==ionizedIndex):
            ion_string  += shell + '{:<2} '.format(input.init_occup[i]-1)
        else: 
            ion_string  += shell + '{:<2} '.format(input.init_occup[i])
    
    blank = ' '*(5 * (8-len(input.shell_list)))
    base_string += blank
    ion_string += blank
    base_string += input.shell_list[ionizedIndex]+'\n'
    ion_string += '\n'
    f = open('crcxin','w')
    # I do not know what these do.
    f.write('200-51 1 2  01.   1.   5.0E-08   1.0E-11-2  0130 0 1.00 0.6510.00 0.5       0.70\n')
    f.write(base_string)
    f.write(ion_string)

    f.write('   -1 \n')
    f.close()
    
    
    return 0 

def run(index,input:inputclass.Input):
    os.system("./crcxl.x")
    os.system("./caiona.x")
    
    os.system("cp crcxout crcxout{}".format(input.shell_list[index]))
    os.system("cp caxout caxout{}".format(input.shell_list[index]))

    return 0

def extract(index):
    
    #extracts ionizaiton potential in the least efficient way possible 
    
    os.system('grep "ionization potential" caxout > iop')
    os.system('grep "incident energy=" caxout > en')
    os.system('grep "total cross section=" caxout > csa')
    os.system('tail -n 17 caxout > csaExtra')

    f = open('iop','r')
    x=  f.readline()
    iop_ev = float(x.split()[2].replace('D','E'))
    f.close()
    
    print(iop_ev)
    
    f = open('csa','r')
    x = f.readlines()
    csa = []
    for string in x:
        csa.append(float(string.split()[3].replace('D','E')))
    f.close()
    
    f = open('en','r')
    x = f.readlines()
    en = []
    for string in x:
        en.append(float(string.split()[2].replace('D','E')))
    f.close()
    print(csa,en)
    
    #f = open('csaExtra','r')
    #x = f.readlines()
    #for string in x:
    #    en. append(float(string.split()[1].replace('D','E'))*iop_ev)
    #    csa.append(float(string.split()[6].replace('D','E')))
#
    #f.close()
    
    f = open('shell'+str(index),'w')
    f.write("# {}\n".format(iop_ev))
    for ii in range(0,len(en)):
        #ignore numerical failures in the fitting/radial mesh
        if csa[ii] > 0:
            f.write("{} {}\n".format(en[ii],csa[ii]))
    f.close()
    
    return 0

def generateIonization(input: inputclass.Input):
    
    for i in range(0,len(input.shell_list)):
        
        run_calc(i,input)
    
    return 0

def run_calc(shellIndex,input: inputclass.Input):
    generateInput(shellIndex,input)
    run(shellIndex,input)
    extract(shellIndex)
    return 0

def produceRates(input: inputclass.Input):
    import numpy as np 
    
    temperature_grid = input.temperatureGridKelvin
    
    rates_vector = [] 
    
    shells_done = []
    for i in range(0,len(input.shell_list)):
        filename= 'shell'+str(i)
        try:
            data = np.loadtxt('shell'+str(i))
            shells_done.append(input.shell_list[i])

        except:
            print('error looking for shell ',i)
            
        f = open(filename)
        x = f.readline()
        f.close()
        ipot = float(x.split()[1])
        #print(ipot)

        #processCSA takes in the final electron energy - i.e after impact.
        #Therefore we subtract the ionization potential.
        #Slightly redundant as I add it back within this routine, which was originally intended
        #to process FAC data which ouputs final electorn energy.
        #
        rates = processCSA(data[:,0]-ipot,data[:,1],temperature_grid,ipot)
        rates_vector.append(rates)

        formatString = "{:10.2e}"
    
    totalRate = np.zeros_like(temperature_grid)
    for rates in rates_vector:
        totalRate = totalRate + rates
    
    tempHeader = "  "
    
    f = open('ratesFile.dat','w')
    
    for ii in range(0,len(temperature_grid)):
        tempHeader += " " + formatString.format(temperature_grid[ii]) 
    f.write(tempHeader+"\n")

    for ii in range(0,len(shells_done)):
        string = input.shell_list[ii]
        for jj in range(0,len(temperature_grid)):
            
            string += " "+formatString.format(rates_vector[ii][jj])
        f.write(string+"\n")
        
    string = 'TO'
    for jj in range(0,len(temperature_grid)):        
        string += " "+formatString.format(totalRate[jj])
    f.write(string+"\n")
            #string += " " + formatString.format(rates_vector[]) 
            #print(shells_done[ii])
    f.close()
    
    f = open('plotRates.dat','w')
    
    header = "#"+9*' '
    for shell in shells_done:
        header += ' '*8+shell 
    header+= '\n' 
    f.write(header)
    
    for jj in range(0,len(temperature_grid)):
        string = formatString.format(temperature_grid[jj]) 
        for ii in range(0,len(shells_done)):
            string+= formatString.format(rates_vector[ii][jj])
        string+='\n'
        f.write(string)
    f.close()
    
    return 0 
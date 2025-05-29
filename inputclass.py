import numpy as np 
import sys 
#Class that stores the input deck and any defaults.

class Input:
    def __init__(self,
                    num_shells=4,
                    shell_list = ['4p','4d','5s','5p'],
                    init_occup = [  6,  10,   2,   4],
                    atomic_num = 52,
                    ionstage = 1,
                    symbol = 'Te',
                    temperatureGridKelvin = [100,1000,2000,5000,10000],
                    runCowan = True,
                    runRates = True
                 ):
        
            self.num_shells=num_shells
            self.shell_list = shell_list
            self.init_occup = init_occup
            self.atomic_num = atomic_num
            self.ionstage = ionstage
            self.symbol = symbol
            
            if len(symbol) > 2:
                print('ATOMIC SYMBOL MUST BE 2 CHARACTERS IN LENGTH.')
                sys.exit()
                
            if num_shells > 8:
                print('Cowan can only have 8 shells, amend structure.')
                sys.exit()
            
            self.temperatureGridKelvin = temperatureGridKelvin
            self.runCowan = runCowan 
            self.runRates = runRates 
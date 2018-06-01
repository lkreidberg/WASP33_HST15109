import sys
sys.path.insert(0, './models')
from constant import constant
from polynomial1 import polynomial1

class Functions:
    def __init__(self, data, funcs):
        self.astro = []
        self.astro_porder = []
        self.sys = []
        self.sys_porder = []

        for f in funcs:
            if f == "constant":
                self.sys.append(constant)
                self.sys_porder.append(data.par_order['c']*data.nvisit)
            elif f == "polynomial1":
                self.sys.append(polynomial1)
                self.sys_porder.append(data.par_order['v']*data.nvisit)
            else:
                #FIXME return error here
                return 0


    #def modelramp(self, t, params):
                    

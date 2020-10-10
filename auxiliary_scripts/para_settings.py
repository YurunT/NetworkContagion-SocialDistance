class para_setting():    
    def __init__(self,):
        self.m = 0.45
        self.T = 0.6
        self.tm1 = 0.3
        self.tm2 = 0.7
        self.msg = 'test'
        self.modelname = 'mask'
        self.itemname = 'es'
        self.change = 0
    
    def print_paras(self,):
        print("m:", self.m)
        print("T:", self.T)
        print("tm1:", self.tm1)
        print("tm2:", self.tm2)
        print("msg:", self.msg)
        print("modelname:", self.modelname)
        print("itemname:", self.itemname)
        print("change:", self.change)
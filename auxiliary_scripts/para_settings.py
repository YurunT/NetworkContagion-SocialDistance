class para_setting():    
    def __init__(self,):
        self.m = 0.45
        self.T = 0.6
        self.tm1 = [0.3, 1]
        self.tm2 = [0.7, 1]
        self.msg = 'test'
        self.modelname = 'mask'
        self.itemname = 'es'
        self.change = 0
        self.n = 50
        self.e = 10
        self.cp = 5
    
    def print_paras(self,):
        print("m:", self.m)
        print("T:", self.T)
        print("tm1:", self.tm1)
        print("tm2:", self.tm2)
        print("msg:", self.msg)
        print("modelname:", self.modelname)
        print("itemname:", self.itemname)
        print("change:", self.change)
        print("n:", self.n)
        print("e:", self.e)
        print("cp:", self.cp)
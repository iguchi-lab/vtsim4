YCC_L:  int = 1
YCC_M:  int = 2
YCC_H:  int = 3
RARA0:  int = 4
RARA1:  int = 5
RARA2:  int = 6
SUMIKA: int = 7

class fan_spec:
    def __init__(self):
        self.qmax = {}
        self.pmax = {}
        self.q1   = {}
        self.p1   = {}

    def set_spec(self, num, qmax, pmax, q1, p1):
        self.qmax[num] = qmax
        self.pmax[num] = pmax
        self.q1[num]   = q1
        self.p1[num]   = p1

    def fan_dict(self, df, dic):
        qmax_d, pmax_d, q1_d, p1_d = {}, {}, {}, {}
        for k in dic:  
            qmax_d[k] = self.qmax[dic[k]]
            pmax_d[k] = self.pmax[dic[k]]
            q1_d[k]   = self.q1[dic[k]]
            p1_d[k]   = self.p1[dic[k]]

        qmax = df.replace(qmax_d).to_list()
        pmax = df.replace(pmax_d).to_list()
        q1   = df.replace(q1_d).to_list()
        p1   = df.replace(p1_d).to_list()
        
        return {'qmax': qmax, 'pmax': pmax, 'q1': q1, 'p1': p1}

fs  = fan_spec()

fs.set_spec(YCC_L,  100 / 3600, 110.0, 100 / 3600,  95.0)
fs.set_spec(YCC_M,  200 / 3600, 110.0, 200 / 3600,  80.0)
fs.set_spec(YCC_H,  390 / 3600, 110.0, 200 / 3600,  80.0)

fs.set_spec(RARA0,  100 / 3600,  10.0,   0 / 3600, 10.0)
fs.set_spec(RARA1,  150 / 3600,  80.0,   0 / 3600, 80.0)
fs.set_spec(RARA2,  210 / 3600, 120.0,   0 / 3600, 120.0)

fs.set_spec(SUMIKA, 200 / 3600, 100.0, 200 / 3600, 100.0)
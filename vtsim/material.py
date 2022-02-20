class Material:
    def __init__(self):
        self.material = {}

    def set_material(self, name, lambda_, capa_):
        self.material[name] = {'lambda': lambda_, 'capa': capa_}

class Wall:
    def __init__(self, mat, wall):
        
        r = [0.0] * len(wall)
        u = [0.0] * len(wall)
        c = [0.0] * len(wall)
        
        for i, l in enumerate(wall):
            for w in l[0]:
                if w[0] == '通気層': 
                    r[i] += 0.110
                    c[i] += w[1] * 1.006 * 1.2 * 1000            #m  *  kJ / (kg・K)  *  kg / m3  *  J / kJ
                elif w[0] == '中空層':
                    r[i] += 0.090
                    c[i] += w[1] * 1.006 * 1.2 * 1000            #m  *  kJ / (kg・K)  *  ka / m3  *  J / kJ
                else:
                    r[i] += w[1] / mat.material[w[0]]['lambda']
                    c[i] += w[1] * mat.material[w[0]]['capa'] * 1000 #J/(m2・K)

            u[i] = 1 / r[i] * l[1]
            c[i] = c[i] * l[1]

        self.u_value = sum(u)
        self.capa_w  = sum(c)
    
    def spec(self):
        return {'U_w': self.u_value, 'capa_w': self.capa_w, 'eta_w': 0.8}

#lambda = W/mK, capa = kJ/(m3・K)
material = {
    '硬質ウレタンフォーム保温板1種2号':     {'lambda': 0.024, 'capa':   43.953},
    '押出法ポリスチレンフォーム3種':        {'lambda': 0.028, 'capa':   25.116},
    '木片セメント板':                       {'lambda': 0.170, 'capa': 1678.59 },
    '合板':                                 {'lambda': 0.160, 'capa':  715.806},
    '吹込用セルロースファイバー断熱材1':    {'lambda': 0.040, 'capa':   37.674},
    'せっこうボード':                       {'lambda': 0.220, 'capa':  904.176},
    '天然木材1類（桧、杉、えぞ松等）':      {'lambda': 0.120, 'capa':  519.064},
    'セメント・モルタル':                   {'lambda': 1.500, 'capa': 1599.05 },
    'コンクリート':                         {'lambda': 1.600, 'capa': 1896.26 }
}

wall_composition= {
    '木造_床_無断熱': [
        [[['合板',                              0.120 ]], 1.00]
    ],
    '木造_間仕切壁_2重中空': [
        [[['せっこうボード',                    0.125 ],
          ['中空層',                            0.100 ],
          ['せっこうボード',                    0.125 ]], 1.00]
    ],
    '木造_天井_無断熱': [
        [[['せっこうボード',                    0.0125]], 1.00]
    ],
    'FPJ_RC_基礎壁_内断熱50+30': [
        [[['押出法ポリスチレンフォーム3種',     0.030 ],
          ['コンクリート',                      0.150 ],
          ['押出法ポリスチレンフォーム3種',     0.050 ],
          ['セメント・モルタル',                0.150 ]], 1.00]
    ],
    'FPJ_RC_床_内断熱30': [
        [[['押出法ポリスチレンフォーム3種',     0.030 ],
          ['コンクリート',                      0.150 ]], 1.00]
    ],
    'FPJ_木造（在来）_外壁_充填+外断熱': [
        [[['せっこうボード',                    0.0125],
          ['吹込用セルロースファイバー断熱材1', 0.120 ],
          ['通気層',                            0.180 ],
          ['木片セメント板',                    0.015 ],
          ['硬質ウレタンフォーム保温板1種2号',  0.030 ]], 0.83],
        [[['せっこうボード',                    0.0125],
          ['天然木材1類（桧、杉、えぞ松等）',   0.120 ],
          ['通気層',                            0.180 ],
          ['木片セメント板',                    0.015 ],
          ['硬質ウレタンフォーム保温板1種2号',  0.030 ]], 0.17]
    ],
    'FPJ_木造（在来）_外壁_充填': [
        [[['せっこうボード',                    0.0125],
          ['吹込用セルロースファイバー断熱材1', 0.120 ],
          ['通気層',                            0.180 ],
          ['木片セメント板',                    0.015 ]], 0.83],
        [[['せっこうボード',                    0.0125],
          ['天然木材1類（桧、杉、えぞ松等）',   0.120 ],
          ['通気層',                            0.180 ],
          ['木片セメント板',                    0.015 ]], 0.17]
    ],        
    'FPJ_木造_天井_吹込': [
        [[['せっこうボード',                    0.0125],
          ['吹込用セルロースファイバー断熱材1', 0.200 ]], 1.00]
    ],
    'FPJ_木造_間仕切壁': [
        [[['押出法ポリスチレンフォーム3種',     0.050 ]], 0.83],
        [[['天然木材1類（桧、杉、えぞ松等）',   0.120 ],
        ['押出法ポリスチレンフォーム3種',       0.050 ]], 0.17]
    ],
    'FPJ_木造_屋根_充填': [
        [[['吹込用セルロースファイバー断熱材1', 0.200 ],
          ['合板',                              0.012 ],
          ['通気層',                            0.030 ],
          ['合板',                              0.012 ]], 0.86],
        [[['天然木材1類（桧、杉、えぞ松等）',   0.200 ],
          ['合板',                              0.012 ],
          ['通気層',                            0.030 ],
          ['合板',                              0.012 ]], 0.14]
    ],
    'FPJ_RC_床_内断熱': [
        [[['押出法ポリスチレンフォーム3種',     0.015 ]], 1.00]
    ]
}

mat = Material()
wall = {}

for m in material:              mat.set_material(m, material[m]['lambda'], material[m]['capa'])
for wl in wall_composition:     wall[wl] = Wall(mat, wall_composition[wl])

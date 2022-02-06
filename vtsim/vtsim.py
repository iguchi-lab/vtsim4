###############################################################################
# import
###############################################################################

import numpy as np
import pandas as pd
import time
import json

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import vtsimc as vt

###############################################################################
# define const
###############################################################################

SOLVE_LU:   int = vt.SOLVE_LU
SOLVE_SOR:  int = vt.SOLVE_SOR

OPT_DF:     int  = 0              #DataFrameを出力
OPT_CSV:    int  = 1              #上記に加え、csvファイルを出力
OPT_GRAPH:  int  = 2              #上記に加えグラフを描画

SN_NONE:    int = vt.SN_NONE
SN_CALC:    int = vt.SN_CALC
SN_FIX:     int = vt.SN_FIX
SN_DLY:     int = vt.SN_DLY
    
VN_SIMPLE:  int = vt.VN_SIMPLE
VN_GAP:     int = vt.VN_GAP
VN_FIX:     int = vt.VN_FIX
VN_AIRCON:  int = vt.VN_AIRCON
VN_FAN:     int = vt.VN_FAN
    
TN_SIMPLE:  int = vt.TN_SIMPLE
TN_AIRCON:  int = vt.TN_AIRCON
TN_SOLAR:   int = vt.TN_SOLAR
TN_GROUND:  int = vt.TN_GROUND
TN_HEATER:  int = vt.TN_HEATER
    
AC_AUTO:    int = vt.AC_AUTO
AC_HEATING: int = vt.AC_HEATING
AC_COOLING: int = vt.AC_COOLING
AC_STOP:    int = vt.AC_STOP

###############################################################################
# define lambda
###############################################################################
read_csv = lambda fn:   pd.read_csv(fn, index_col = 0, parse_dates = True).fillna(method = 'bfill')\
                                                                          .fillna(method = 'ffill')     #csvファイルの読み込み

#index = lambda df:     df.index.strftime('%Y/%m/%d %H:%M:%S').to_list()
#data  = lambda df:     df.to_list()

df   = lambda length:   pd.DataFrame(index = pd.date_range(datetime(2021, 1, 1, 0, 0, 0), 
                                     datetime(2021, 1, 1, 0, 0, 0) + timedelta(seconds = length), 
                                     freq='1s'))                                                        #長さlength、1s毎の時刻
                                                          
d_node  = lambda name:  name + '_c'                                                                     #遅延ノードの名前作成

calc = vt.VTSim()

def encode(object):
    if isinstance(object, pd.core.indexes.datetimes.DatetimeIndex): return(object.strftime('%Y/%m/%d %H:%M:%S').to_list())
    if isinstance(object, pd.core.series.Series):                   return(object.to_list())

def write_json(input, fn):
    with open(fn, 'w') as f:
        json.dump(input, f, default = encode, ensure_ascii = False, indent = 4)

def to_json(input):
    return(json.dumps(input, default = encode, ensure_ascii = False, indent = 4))

df_p   = pd.DataFrame()
df_c   = pd.DataFrame()
df_t   = pd.DataFrame()
df_qv  = pd.DataFrame()
df_qt1 = pd.DataFrame()
df_qt2 = pd.DataFrame()

###############################################################################
# define function
###############################################################################
def to_list_f(v):
    if   type(v) == list:                   return(v)
    elif type(v) == np.ndarray:             return(v)
    elif type(v) == pd.core.series.Series:  return(np.array(v))
    else:                                   return[float(v)] * calc.sts.length

def to_list_i(v):
    if   type(v) == list:                   return(v)
    elif type(v) == np.ndarray:             return(v)
    elif type(v) == pd.core.series.Series:  return(np.array(v))
    else:                                   return[int(v)] * calc.sts.length

def run_calc(inp):                                                     #はじめに呼び出される関数    
    input = json.loads(inp)
    
    print('Set calc status.')
    if 'index' in input:    set_calc_status(input)  
    else:                   raise Exception('ERROR: index does not exist!')

    print('Add Capacity')
    if 'sn' in input:       input = add_capa(input)
    else:                   raise Exception('ERROR: sn does not exist!')

    with open('calc.json', 'w') as f:
        json.dump(input, f, ensure_ascii = False, indent = 4)


    print('Set SimNode.')
    if 'sn' in input:       set_sim_node(input['sn'])
    else:                   raise Exception('ERROR: sn does not exist!')

    print('Set VentNet.')
    if 'vn' in input:       set_vent_net(input['vn'])

    print('Set ThrmNet.')
    if 'tn' in input:       set_thrm_net(input['tn'])

    print('ready')
    print('sts     ', calc.sts)
    print('sn      ', calc.sn)
    print('node    ', calc.node)
    print('vn      ', calc.vn)
    print('tn      ', calc.tn)
    print('v_idc   ', calc.v_idc)
    print('c_idc   ', calc.c_idc)
    print('t_idc   ', calc.t_idc)
    print('i_vn_ac ', calc.i_vn_ac)
    print('i_tn_ac ', calc.i_tn_ac)

    print('Start vtsim calc.')
    s_time = time.time()
    calc.calc()
    e_time = time.time() - s_time    
    print('Finish vtsim calc.')
    print("calc time = {0}".format(e_time * 1000) + "[ms]")

    opt = input['opt'] if 'opt' in input else OPT_GRAPH

    return output_calc(input['index'], opt, calc.result(), input['sn'].keys(), input['vn'].keys(), input['tn'].keys())

def set_calc_status(input):
    sts  = vt.CalcStatus()

    ix = pd.to_datetime(input['index'])
    sts.length = len(ix)
    sts.t_step = (ix[1] - ix[0]).seconds + (ix[1] - ix[0]).microseconds / 1000000   #t_stepの読み込み    

    if 'calc_sts' in input: calc_sts = input['calc_sts']
    else:                   calc_sts = {}
    
    if 'solve'     in calc_sts:   sts.solve     = calc_sts['solve']
    if 'step_p'    in calc_sts:   sts.step_p    = calc_sts['step_p']
    if 'vent_err'  in calc_sts:   sts.vent_err  = calc_sts['vent_err']
    if 'step_t'    in calc_sts:   sts.step_t    = calc_sts['step_t']
    if 'thrm_err'  in calc_sts:   sts.thrm_err  = calc_sts['thrm_err']
    if 'conv_err'  in calc_sts:   sts.conv_err  = calc_sts['conv_err']
    if 'sor_ratio' in calc_sts:   sts.sor_ratio = calc_sts['sor_ratio']
    if 'sor_err'   in calc_sts:   sts.sor_err   = calc_sts['sor_err']

    calc.setup(sts)

def add_capa(input):
    if 'tn' not in input:   input['tn'] = {}
    
    for n in [n for n in input['sn'] if 'capa' in input['sn'][n]]:                              #熱容量の設定のあるノード
        input['sn'][d_node(n)] = {}
        input['sn'][d_node(n)]['t_flag'] = vt.SN_DLY                                            #計算フラグ
        input['sn'][d_node(n)]['s_i']    = n                                                    #親ノードの設定
        if 't' in input['sn'][n]:   input['sn'][d_node(n)]['t'] = input['sn'][n]['t']           #初期温度の継承

        input['tn'][n + ' -> ' + d_node(n)] = {}
        input['tn'][n + ' -> ' + d_node(n)]['type'] = vt.TN_SIMPLE                              #熱容量の設定
        input['tn'][n + ' -> ' + d_node(n)]['cdtc'] = input['sn'][n]['capa'] / calc.sts.t_step  #コンダクタンス（熱容量）     

    return input

def set_sim_node(sn):
    for i, n in enumerate(sn):
        calc.set_node(n, i)
        v_flag = sn[n]['v_flag'] if 'v_flag' in sn[n] else vt.SN_NONE
        c_flag = sn[n]['c_flag'] if 'c_flag' in sn[n] else vt.SN_NONE
        t_flag = sn[n]['t_flag'] if 't_flag' in sn[n] else vt.SN_NONE
        calc.sn_add(i, [v_flag, c_flag, t_flag])

        if 'p'     in sn[n]:    calc.sn[i].p     = to_list_f(sn[n]['p'])                    #圧力、行列で設定可能
        if 'c'     in sn[n]:    calc.sn[i].c     = to_list_f(sn[n]['c'])                    #濃度、行列で設定可能
        if 't'     in sn[n]:    calc.sn[i].t     = to_list_f(sn[n]['t'])                    #温度、行列で設定可能
        if 'h_sr'  in sn[n]:    calc.sn[i].h_sr  = to_list_f(sn[n]['h_sr'])                 #日射量、行列で設定可能
        if 'h_inp' in sn[n]:    calc.sn[i].h_inp = to_list_f(sn[n]['h_inp'])                #発熱、行列で設定可能
        if 'v'     in sn[n]:    calc.sn[i].v     = to_list_f(sn[n]['v'])                    #気積、行列で設定可能
        if 'm'     in sn[n]:    calc.sn[i].m     = to_list_f(sn[n]['m'])                    #発生量、行列で設定可能
        if 'beta'  in sn[n]:    calc.sn[i].beta  = to_list_f(sn[n]['beta'])                 #濃度減少率、行列で設定可能
        if 's_i'   in sn[n]:    calc.sn[i].s_i   = calc.node[sn[n]['s_i']]

def set_vent_net(vn):
    for i, nt in enumerate(vn):
        #vn_type = vn[nt]['type'] if 'type' in vn[nt] else vt.VN_FIX
        
        if 'type' not in vn[nt]:
            if   ('alpha' in vn[nt]) and ('area' in vn[nt]):  vn_type = vt.VN_SIMPLE
            elif ('a'     in vn[nt]) and ('n'    in vn[nt]):  vn_type = vt.VN_GAP
            elif ('qmax'  in vn[nt]) and ('pmax' in vn[nt]):  vn_type = vt.VN_FAN
            elif  'vol'   in vn[nt]:                          vn_type = vt.VN_FIX
            else:                                             raise Exception('ERROR: vn_type dose not exist! ' + nt)    
        else:
            vn_type = vn[nt]['type']

        h1 = to_list_f(vn[nt]['h1']) if 'h1' in vn[nt] else to_list_f(0.0)                          #高さ1、行列設定不可
        h2 = to_list_f(vn[nt]['h2']) if 'h2' in vn[nt] else to_list_f(0.0)                          #高さ2、行列設定不可
        
        s = nt.replace(' ', '')
        if s.find('->') == -1:  raise Exception('ERROR: -> does not exist!')
        if s.find(':')  == -1:  n1, n2 = s[:s.find('->')], s[s.find('->') + 2:]
        else:                   n1, n2 = s[:s.find('->')], s[s.find('->') + 2: s.find(':')]

        calc.vn_add(i, calc.node[n1], calc.node[n2], vn_type, h1, h2)
        
        if (vn_type == vt.VN_FIX) or (vn_type == vt.VN_AIRCON):       
            calc.vn[i].qv = to_list_f(vn[nt]['vol'])                                                #風量固定値、行列で設定可能
        if vn_type == vt.VN_SIMPLE:    
            calc.vn[i].alpha = to_list_f(vn[nt]['alpha'])
            calc.vn[i].area  = to_list_f(vn[nt]['area'])                                            #単純開口、行列で設定可能
        if vn_type == vt.VN_GAP:           
            calc.vn[i].a     = to_list_f(vn[nt]['a'])
            calc.vn[i].n     = to_list_f(vn[nt]['n'])                                               #隙間、行列で設定可能
        if vn_type == vt.VN_FAN:           
            calc.vn[i].q_max = to_list_f(vn[nt]['qmax']) 
            calc.vn[i].p_max = to_list_f(vn[nt]['pmax']) 
            calc.vn[i].q1    = to_list_f(vn[nt]['q1'])
            calc.vn[i].p1    = to_list_f(vn[nt]['p1'])                                              #ファン、行列で設定可能
        if vn_type == vt.VN_AIRCON:
            calc.i_vn_ac = i
        calc.vn[i].eta = to_list_f(vn[nt]['eta']) if 'eta' in vn[nt] else to_list_f(0.0)            #粉じん除去率、行列で設定可能

def set_thrm_net(tn):
    for i, nt in enumerate(tn):
        #tn_type = tn[nt]['type'] if 'type' in tn[nt] else vt.TN_SIMPLE

        if 'type' not in tn[nt]:
            if    'ms'       in tn[nt]:                            tn_type = vt.TN_SOLAR
            elif ('phi_0'    in tn[nt]) and ('cof_r'   in tn[nt]): tn_type = vt.TN_GROUND
            elif ('ac_mode'  in tn[nt]) and ('pre_tmp' in tn[nt]): tn_type = vt.TN_AIRCON
            elif  'cdtc'     in tn[nt]:                            tn_type = vt.TN_SIMPLE
            else:                                                  raise Exception('ERROR: tn_type dose not exist! ' + nt) 
        else:
            tn_type = tn[nt]['type']

        s = nt.replace(' ', '')
        calc.tn_add(i, calc.node[s[:s.find('->')]], calc.node[s[s.find('->') + 2:]], tn_type)

        if tn_type == vt.TN_SIMPLE:     
            calc.tn[i].cdtc = to_list_f(tn[nt]['cdtc'])                                #コンダクタンス、行列設定可能
        if tn_type == vt.TN_AIRCON:     
            calc.tn[i].ac_mode = to_list_i(tn[nt]['ac_mode']) 
            calc.tn[i].pre_tmp = to_list_f(tn[nt]['pre_tmp'])                          #エアコン運転モード
            calc.i_tn_ac = i
        if tn_type == vt.TN_SOLAR:       
            calc.tn[i].ms      = to_list_f(tn[nt]['ms'])                               #日射熱取得率、行列設定可能
        if tn_type == vt.TN_GROUND:     
            calc.tn[i].area    = to_list_f(tn[nt]['area'])           
            calc.tn[i].rg      = to_list_f(tn[nt]['rg']) 
            calc.tn[i].phi_0   = tn[nt]['phi_0']
            calc.tn[i].cof_r   = tn[nt]['cof_r']
            calc.tn[i].cof_phi = tn[nt]['cof_phi']                                     #地盤熱応答、行列設定不可（面積と断熱性能はOK）         

def output_calc(ix, opt, res, sn_c, vn_c, tn_c):
    dat_list  = [{'df': pd.DataFrame(), 'columns': sn_c, 'fn': 'vent_p.csv',   'title': '圧力',  'unit': '[Pa]'},
                 {'df': pd.DataFrame(), 'columns': sn_c, 'fn': 'vent_c.csv',   'title': '濃度',  'unit': '[個/L]'},
                 {'df': pd.DataFrame(), 'columns': sn_c, 'fn': 'them_t.csv',   'title': '温度',  'unit': '[℃]'},
                 {'df': pd.DataFrame(), 'columns': vn_c, 'fn': 'vent_qv.csv',  'title': '風量',  'unit': '[m3/s]'},
                 {'df': pd.DataFrame(), 'columns': vn_c, 'fn': 'thrm_qt1.csv', 'title': '熱量1', 'unit': '[W]'},
                 {'df': pd.DataFrame(), 'columns': tn_c, 'fn': 'thrm_qt2.csv', 'title': '熱量2', 'unit': '[W]'}]
    
    ix = pd.to_datetime(ix, format='%Y/%m/%d %H:%M:%S')

    for i, d in enumerate(dat_list):
        if len(res[i]) != 0: d['df'] = pd.DataFrame(np.array(res[i]).T,  index = ix, columns = d['columns'])

    if opt >= OPT_CSV:
        print('Output csv files.')
        for d in dat_list:      d['df'].to_csv(d['fn'], encoding = 'utf_8_sig')

    if opt >= OPT_GRAPH:
        print('Draw Graphs.')
        fig = plt.figure(facecolor = 'w', figsize = (18, len(dat_list) * 4))
        fig.subplots_adjust(wspace = -0.1, hspace=0.9)

        for i, d in enumerate(dat_list):
            a = fig.add_subplot(len(dat_list), 1, i + 1)
            for cl in d['df'].columns:
                a.plot(d['df'][cl], linewidth = 1.0, label = cl)
            a.legend(ncol = 5, bbox_to_anchor = (0, 1.05, 1, 0), 
                     loc = 'lower right', borderaxespad = 0, facecolor = 'w', edgecolor = 'k')
            a.set_title(d['title'], loc='left')
            a.set_ylabel(d['unit'])

    global df_p, df_c, df_t, df_qv, df_qt1, df_qt2
    df_p   = dat_list[0]['df']
    df_c   = dat_list[1]['df']
    df_t   = dat_list[2]['df'] 
    df_qv  = dat_list[3]['df']
    df_qt1 = dat_list[4]['df']
    df_qt2 = dat_list[5]['df']

    return dat_list[0]['df'], dat_list[1]['df'], dat_list[2]['df'], dat_list[3]['df'], dat_list[4]['df'], dat_list[5]['df']
###############################################################################
# import
###############################################################################
import numpy as np
import pandas as pd
import time
import json
import logging

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import vtsimc as vt

logger = logging.getLogger(__name__)
streamhandler = logging.StreamHandler()
logger.setLevel(logging.INFO)
format = '[%(asctime)s][%(levelname)s][%(message)s]'
datefmt='%Y/%m/%d %I:%M:%S'
formatter = logging.Formatter(format, datefmt)
streamhandler.setFormatter(formatter)


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
# define lambda etc.
###############################################################################
read_csv = lambda fn:   pd.read_csv(fn, index_col = 0, parse_dates = True).fillna(method = 'bfill')\
                                                                          .fillna(method = 'ffill')     #csvファイルの読み込み

df   = lambda length:   pd.DataFrame(index = pd.date_range(datetime(2021, 1, 1, 0, 0, 0), 
                                     datetime(2021, 1, 1, 0, 0, 0) + timedelta(seconds = length), 
                                     freq='1s'))                                                        #長さlength、1s毎の時刻

def encode(object):                                                                                     #JSON形式に変換する際のエンコード
    if isinstance(object, pd.core.indexes.datetimes.DatetimeIndex): 
        return(object.strftime('%Y/%m/%d %H:%M:%S').to_list())
    if isinstance(object, pd.core.series.Series):                   
        return(object.to_list())

def read_json(fn):                                                                                      #JSON形式を読み込みdict型で返す
    with open(fn) as f:
        input = json.load(f)
    return input

def write_json(input, fn):                                                                              #dict型をJSONファイルに書き出し
    input['version'] = '4.0.0'
    with open(fn, 'w') as f:
        json.dump(input, f, default = encode, ensure_ascii = False, indent = 4)

def to_json(input):                                                                                     #dict型をJSON形式に変換
    input['version'] = '4.0.0'
    return(json.dumps(input, default = encode, ensure_ascii = False, indent = 4))                       

calc = vt.VTSim()

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

def run_calc(input):                                                                    #はじめに呼び出される関数    
    if   type(input) == dict: input = to_json(input)                                    #辞書型であれば、JSON形式に変換
    elif type(input) != str: raise Exception('ERROR: inputは、辞書型かJSON形式である必要がります。')      #文字列（JSON形式)で無ければエラー
    
    input = json.loads(input)                                                           #JSON形式を辞書型に変換
    if 'vn' not in input:   input['vn'] = {}
    if 'tn' not in input:   input['tn'] = {}

    logger.info('Set calc status.')
    if 'index' in input:    set_calc_status(input)                                      #計算条件を設定
    else:                   raise Exception('ERROR: index が存在しません。')             #indexが無ければエラー

    logger.info('Add Capacity.')   
    if 'sn' in input:       input = add_capa(input)                                     #熱容量を設定
    else:                   raise Exception('ERROR: ノード(sn)が存在しません。')         #sn（ノード）が無ければエラー

    logger.info('Set Aircon.')
    if 'aircon' in input:   input = set_aircon(input)                                   #エアコンをセット

    with open('calc.json', 'w') as f:                                                   #計算入力を　calc.jsonに格納
        json.dump(input, f, ensure_ascii = False, indent = 4)

    logger.info('Set SimNode.')
    if 'sn' in input:       set_sim_node(input['sn'])                                   #sn（ノード）の設定
    else:                   raise Exception('ERROR: ノード(sn)が存在しません。')            #sn（ノード）が無ければエラー

    logger.info('Set VentNet.')
    if 'vn' in input:       set_vent_net(input['vn'])                                   #vn（換気回路網）の設定

    logger.info('Set ThrmNet.')
    if 'tn' in input:       set_thrm_net(input['tn'])                                   #tn（熱回路網）の設定

    logger.info('ready')
    logger.info('sts     ', calc.sts)
    logger.info('sn      ', calc.sn)
    logger.info('node    ', calc.node)
    logger.info('vn      ', calc.vn)
    logger.info('tn      ', calc.tn)
    logger.info('v_idc   ', calc.v_idc)
    logger.info('c_idc   ', calc.c_idc)
    logger.info('t_idc   ', calc.t_idc)
    logger.info('i_vn_ac ', calc.i_vn_ac)
    logger.info('i_tn_ac ', calc.i_tn_ac)

    logger.info('Start vtsim calc.')
    s_time = time.time()
    calc.calc()                                                                         #計算
    e_time = time.time() - s_time    
    logger.info('Finish vtsim calc.')
    logger.info("calc time = {0}".format(e_time * 1000) + "[ms]")

    opt = input['opt'] if 'opt' in input else OPT_GRAPH

    dat_list  = [{'columns': input['sn'].keys(), 'fn': 'vent_p.csv',   'title': '圧力',  'unit': '[Pa]'},
                 {'columns': input['sn'].keys(), 'fn': 'vent_c.csv',   'title': '濃度',  'unit': '[個/L]'},
                 {'columns': input['sn'].keys(), 'fn': 'them_t.csv',   'title': '温度',  'unit': '[℃]'},
                 {'columns': input['vn'].keys(), 'fn': 'vent_qv.csv',  'title': '風量',  'unit': '[m3/s]'},
                 {'columns': input['vn'].keys(), 'fn': 'thrm_qt1.csv', 'title': '熱量1', 'unit': '[W]'},
                 {'columns': input['tn'].keys(), 'fn': 'thrm_qt2.csv', 'title': '熱量2', 'unit': '[W]'}]

    ix = pd.to_datetime(input['index'], format='%Y/%m/%d %H:%M:%S')
    res = calc.result()

    for i, d in enumerate(dat_list):
        if len(res[i]) != 0: d['df'] = pd.DataFrame(np.array(res[i]).T,  index = ix, columns = d['columns'])

    global df_p, df_c, df_t, df_qv, df_qt1, df_qt2
    df_p   = dat_list[0]['df']
    df_c   = dat_list[1]['df']
    df_t   = dat_list[2]['df'] 
    df_qv  = dat_list[3]['df']
    df_qt1 = dat_list[4]['df']
    df_qt2 = dat_list[5]['df']

    output_calc(dat_list, opt)

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
    for n in [n for n in input['sn'] if 'capa' in input['sn'][n]]:                              #熱容量の設定のあるノード
        nc = n + '_c'
        
        input['sn'][nc] = {}
        input['sn'][nc]['t_flag'] = vt.SN_DLY                                            #計算フラグ
        input['sn'][nc]['s_i']    = n                                                    #親ノードの設定
        if 't' in input['sn'][n]:   input['sn'][nc]['t'] = input['sn'][n]['t']           #初期温度の継承

        input['tn'][n + ' -> ' + nc] = {}
        input['tn'][n + ' -> ' + nc]['type'] = vt.TN_SIMPLE                              #熱容量の設定
        input['tn'][n + ' -> ' + nc]['cdtc'] = input['sn'][n]['capa'] / calc.sts.t_step  #コンダクタンス（熱容量）     

    return input

def set_aircon(input):
    aircon = input['aircon']

    for i, ac in enumerate([ac for ac in input['aircon']]):
    
        if i == 1:  raise Exception('ERROR: エアコンは2台以上設置できません')

        ac_in, ac_out = ac + '_in', ac + '_out'  
        n1, n2, n3 = aircon[ac]['in'], aircon[ac]['out'], aircon[ac]['set']

        input['sn'][ac_in]            = {}
        input['sn'][ac_out]           = {'t_flag': vt.SN_CALC}

        vol = aircon[ac]['vol'] if 'vol' in aircon[ac] else 1000 / 3600

        input['vn'][n1     + ' -> ' + ac_in]  = {'vol': vol}
        input['vn'][ac_in  + ' -> ' + ac_out] = {'ac_vol': vol}
        input['vn'][ac_out + ' -> ' + n2]     = {'vol': vol}

        input['tn'][n3  + ' -> ' + ac_out] = {'ac_mode': aircon[ac]['ac_mode'], 'pre_tmp': aircon[ac]['pre_tmp']}
    
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

def get_n1n2(nt):  
    s = nt.replace(' ', '')
    if s.find('->') == -1:  raise Exception('ERROR: vnもしくはtnのキーに -> が存在しません')
    if s.find(':')  == -1:  n1, n2 = s[:s.find('->')], s[s.find('->') + 2:]
    else:                   n1, n2 = s[:s.find('->')], s[s.find('->') + 2: s.find(':')]

    if n1 not in calc.node: raise Exception('ERROR: ノード(sn)の中に ' + n1 + ' がありません。')
    if n2 not in calc.node: raise Exception('ERROR: ノード(sn)の中に ' + n2 + ' がありません。')

    return n1, n2

def set_vent_net(vn):
    for i, nt in enumerate(vn):
        
        if 'type' not in vn[nt]:
            if   ('alpha'  in vn[nt]) and ('area' in vn[nt]):  vn_type = vt.VN_SIMPLE
            elif ('a'      in vn[nt]) and ('n'    in vn[nt]):  vn_type = vt.VN_GAP
            elif ('qmax'   in vn[nt]) and ('pmax' in vn[nt]):  vn_type = vt.VN_FAN
            elif  'vol'    in vn[nt]:                          vn_type = vt.VN_FIX
            elif  'ac_vol' in vn[nt]:                          vn_type = vt.VN_AIRCON
            else:                                             raise Exception('ERROR: ' + nt + 'のvn_typeを認識できません')    
        else:
            vn_type = vn[nt]['type']

        h1 = to_list_f(vn[nt]['h1']) if 'h1' in vn[nt] else to_list_f(0.0)                          #高さ1、行列設定不可
        h2 = to_list_f(vn[nt]['h2']) if 'h2' in vn[nt] else to_list_f(0.0)                          #高さ2、行列設定不可
        
        n1, n2 = get_n1n2(nt)
        calc.vn_add(i, calc.node[n1], calc.node[n2], vn_type, h1, h2)
        
        if vn_type == vt.VN_FIX:       
            calc.vn[i].qv = to_list_f(vn[nt]['vol'])                                                #風量固定値、行列で設定可能
        if vn_type == vt.VN_AIRCON:
            calc.vn[i].qv = to_list_f(vn[nt]['ac_vol'])                                             #風量固定値、行列で設定可能
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

        if 'type' not in tn[nt]:
            if    'ms'       in tn[nt]:                            tn_type = vt.TN_SOLAR
            elif ('phi_0'    in tn[nt]) and ('cof_r'   in tn[nt]): tn_type = vt.TN_GROUND
            elif ('ac_mode'  in tn[nt]) and ('pre_tmp' in tn[nt]): tn_type = vt.TN_AIRCON
            elif  'cdtc'     in tn[nt]:                            tn_type = vt.TN_SIMPLE
            else:                                                  raise Exception('ERROR: ' + nt + 'のtn_typeを認識できません') 
        else:
            tn_type = tn[nt]['type']

        n1, n2 = get_n1n2(nt)
        calc.tn_add(i, calc.node[n1], calc.node[n2], tn_type)

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

def output_calc(dat_list, opt):

    if opt == OPT_CSV or opt == OPT_GRAPH:
        logger.info('Output csv files.')
        for d in dat_list:      d['df'].to_csv(d['fn'], encoding = 'utf_8_sig')

    if opt == OPT_GRAPH:
        logger.info('Draw Graphs.')
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
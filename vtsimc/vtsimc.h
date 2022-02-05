#include "vtsim_set.h"              //計算の設定　ヘッダーファイルの読み込み
#include "node.h"                   //ノード　ヘッダーファイルの読み込み
#include "vent_net.h"               //換気回路網　ヘッダーファイルの読み込み
#include "thrm_net.h"               //熱回路網　ヘッダーファイルの読み込み
#include "mymath.h"
#include <iostream>
#include <fstream>

//#define DEBUG_ON

#ifdef  DEBUG_ON
#define LOG_PRINT(...)     ofs << __FILE__ << " (" << __LINE__ << ") " << __func__ << ":" << __VA_ARGS__
#define LOG_CONTENTS(...)  ofs << __VA_ARGS__
const char *logfileName = "log.txt";
ofstream ofs(logfileName);
#else
#define LOG_PRINT(...)
#define LOG_CONTENTS(...)
#endif

using namespace::std;

class CalcStatus{
public:
    long length      = 0;
    double t_step    = 3600;
    int solve        = SOLVE_LU;
    double step_p    = STEP_P, 
           vent_err  = VENT_ERR, 
           step_t    = STEP_T, 
           thrm_err  = THRM_ERR, 
           conv_err  = CONV_ERR, 
           sor_ratio = SOR_RATIO, 
           sor_err   = SOR_ERR;
};

class VTSim{
public:
    CalcStatus sts;
    vector<Node> sn;                                        //ノード
    map<string, int> node;
    vector<Vent_Net> vn;                                    //換気回路網
    vector<Thrm_Net> tn;                                    //熱回路網
    vector<int> v_idc, c_idc, t_idc;
    int i_vn_ac = -1, i_tn_ac = -1;

    void setup(CalcStatus sts_){
        sts = sts_;
        sn.clear();
        node.clear();
        vn.clear();
        tn.clear();
        v_idc.clear();
        c_idc.clear();
        t_idc.clear();
        i_vn_ac = -1;
        i_tn_ac = -1;
    }

    void sn_add(int i, tuple<int, int, int> flag){
        sn.push_back(Node(sts.length, i, flag));
        if(get<0>(flag) == SN_CALC)     v_idc.push_back(i);
        if(get<1>(flag) == SN_CALC)     c_idc.push_back(i);
        if(get<2>(flag) == SN_CALC)     t_idc.push_back(i);
    }

    void set_node(string name, int i){
        node[name] = i;
    }

    void vn_add(int i, int i1, int i2, int vn_type, vector<double> h1, vector<double> h2){
        vn.push_back(Vent_Net(sts.length, i, i1, i2, vn_type, h1, h2));
    }

    void tn_add(int i, int i1, int i2, int tn_type){
        tn.push_back(Thrm_Net(sts.length, i, i1, i2, tn_type));
    }

    void change_sn_t_flag(int i, int flag_){
        get<2>(sn[i].flag) = flag_;
        t_idc.clear();
        for(unsigned int i = 0; i < sn.size(); i++)   
            if(get<2>(sn[i].flag) == SN_CALC)     t_idc.push_back(sn[i].i);
    }

    void set_aircon_status(long ts, int flag_){
        if(tn[i_tn_ac].aircon_on[ts] == AC_OFF && flag_ == AC_ON){
            tn[i_tn_ac].aircon_on[ts] = AC_ON;
            change_sn_t_flag(tn[i_tn_ac].i1, SN_FIX);
            vn[i_vn_ac].vn_type = VN_AIRCON;
        }
        else if(tn[i_tn_ac].aircon_on[ts] == AC_ON && flag_ == AC_OFF){
            tn[i_tn_ac].aircon_on[ts] = AC_OFF;
            change_sn_t_flag(tn[i_tn_ac].i1, SN_CALC);  
            vn[i_vn_ac].vn_type = VN_FIX;
        }
    }

    vector<double> qv_sum(vector<double> p, long ts, int flag){
        vector<double> qvsum(sn.size(), 0.0);                                                                             //風量収支の初期化                                 
        for(unsigned int i = 0; i < vn.size(); i++){
            double rgh1 = get_rho(sn[vn[i].i1].t[ts]) * G * vn[i].h1[ts];
            double rgh2 = get_rho(sn[vn[i].i2].t[ts]) * G * vn[i].h2[ts];
            double qv   = vn[i].get_qv((p[vn[i].i1] - rgh1) - (p[vn[i].i2] - rgh2), ts);
            qvsum[vn[i].i1] -= qv;                                                                                        //風量収支の計算
            qvsum[vn[i].i2] += qv;                                                                                        //風量収支の計算      
        }

        for(unsigned int i = 0; i << sn.size(); i++)
            if(flag == 0)   LOG_PRINT("i = " << i << " : " << qvsum[i] << endl);
        return qvsum;
    }

    double calc_vent(long ts){
        vector<double>          p0(sn.size()), qvsum_0(sn.size()), qvsum_d(sn.size());
        vector<vector<double>>  a(v_idc.size(), vector<double>(v_idc.size()));
        vector<double>          b(v_idc.size()), dp(v_idc.size());
        double                  rmse;
        int                     itr = 0;
        
        for(unsigned int i = 0; i < sn.size(); i++)    p0[sn[i].i] = sn[i].p[ts];                                                //圧力の初期化       
        do{
            rmse = 0.0;
            qvsum_0 = qv_sum(p0, ts, 0);
            
            for(unsigned int j = 0; j < v_idc.size(); j++){
                p0[v_idc[j]] += sts.step_p;                                                                             //ダミー圧力の作成       
                qvsum_d = qv_sum(p0, ts, 1);                                                                               //ダ三－風量収支の計算
                for(unsigned int i = 0; i < v_idc.size(); i++)  a[i][j] = (qvsum_d[v_idc[i]] - qvsum_0[v_idc[i]]) / sts.step_p;  //aの計算
                b[j] = -qvsum_0[v_idc[j]];
                p0[v_idc[j]] -= sts.step_p; 
            }

            if(sts.solve == SOLVE_SOR)  dp = SOR(a, b, v_idc.size(), sts.sor_ratio, sts.sor_err);                       //SOR法による計算     
            else                        dp = LU(a, b, v_idc.size()); 

            for(unsigned int i = 0; i < v_idc.size(); i++)   p0[v_idc[i]] += dp[i];                                              //圧力の更新
            qvsum_0 = qv_sum(p0, ts, 2);    
            for(unsigned int i = 0; i < v_idc.size(); i++)   rmse += pow(qvsum_0[v_idc[i]], 2.0) / v_idc.size();
            LOG_PRINT(itr << ": ts = " << ts << ": rmse = " << sqrt(rmse) << endl);
            itr++;
        }while(sts.vent_err < sqrt(rmse));
        for(unsigned int i = 0; i < v_idc.size(); i++)   sn[v_idc[i]].p[ts] = p0[v_idc[i]];                                      //圧力の計算                                                               
        return rmse;
    }

    vector<double> qt_sum(vector<double> t, long ts){
        vector<double> qtsum(sn.size(), 0.0);           
       
        for(unsigned int i = 0; i < vn.size(); i++){                                                                             //移流に伴う熱移動
            if(vn[i].qv[ts] > 0)            qtsum[vn[i].i2] += vn[i].get_qt(t[vn[i].i1] - t[vn[i].i2], ts);
            else                            qtsum[vn[i].i1] += vn[i].get_qt(t[vn[i].i1] - t[vn[i].i2], ts);
        }

        for(unsigned int i = 0; i < tn.size(); i++){                                                                             //貫流、日射、発熱による熱移動
            switch(tn[i].tn_type){
                case TN_SIMPLE:
                case TN_GROUND:
                    qtsum[tn[i].i1] -= tn[i].get_qt(t[tn[i].i1] - t[tn[i].i2], ts);     
                    qtsum[tn[i].i2] += tn[i].get_qt(t[tn[i].i1] - t[tn[i].i2], ts);
                    break;
                case TN_SOLAR:
                    qtsum[tn[i].i1] += tn[i].get_qt_s(sn[tn[i].i2].h_sr[ts], ts);
                    break;
                case TN_HEATER:
                    qtsum[tn[i].i1] += tn[i].get_qt_h(sn[tn[i].i2].h_inp[ts], ts);
                    break;
            }    
        }

        if(i_tn_ac != -1 && tn[i_tn_ac].aircon_on[ts] == AC_ON){                                                        //エアコンによる温度調整
            qtsum[tn[i_tn_ac].i2] -= qtsum[tn[i_tn_ac].i1];
            qtsum[tn[i_tn_ac].i1]  = 0.0;
        }
        return qtsum;
    }

    double calc_thrm(long ts){
        vector<double>          t0(sn.size()), qtsum_0(sn.size()), qtsum_d(sn.size());
        vector<vector<double>>  a(t_idc.size() + 1, vector<double>(t_idc.size() + 1));                                  //エアコン分も余計に確保
        vector<double>          b(t_idc.size() + 1), dt(t_idc.size() + 1);                                              //エアコン分も余計に確保
        double                  rmse;        
        int                     itr = 0;

        for(unsigned int i = 0; i < sn.size(); i++)    t0[sn[i].i] = sn[i].t[ts];                                                //温度の初期化
        
        do{
            if(i_tn_ac != -1){
                set_aircon_status(ts, AC_ON);                                                                           //全部ONになってない？？
                if(itr == 0)    set_aircon_status(ts, AC_OFF);
                else{
                    switch(tn[i_tn_ac].ac_mode[ts]){
                        case AC_AUTO:       set_aircon_status(ts, AC_ON);
                                            break;
                        case AC_COOLING:    if(t0[tn[i_tn_ac].i1] >= tn[i_tn_ac].pre_tmp[ts])   set_aircon_status(ts, AC_ON);
                                            break;
                        case AC_HEATING:    if(t0[tn[i_tn_ac].i1] <= tn[i_tn_ac].pre_tmp[ts])   set_aircon_status(ts, AC_ON);
                                            break;
                    }
                }
            }

            rmse = 0.0;
            qtsum_0 = qt_sum(t0, ts);                                                                                   //熱量収支の計算

            for(unsigned int j = 0; j < t_idc.size(); j++){
                t0[t_idc[j]] += sts.step_t;                                                                             //ダミー温度の作成
                qtsum_d = qt_sum(t0, ts);                                                                               //ダミー熱量の計算
                for(unsigned int i = 0; i < t_idc.size(); i++)  a[i][j] = (qtsum_d[t_idc[i]] - qtsum_0[t_idc[i]]) / sts.step_t;  //aの計算
                b[j] = -qtsum_0[t_idc[j]];                                                                              //bの計算
                t0[t_idc[j]] -= sts.step_t;                                                                             //ダミー温度を戻す
            }

            if(sts.solve == SOLVE_SOR)  dt = SOR(a, b, t_idc.size(), sts.sor_ratio, sts.sor_err);                                   //SOR法による計算
            else                        dt = LU(a, b, t_idc.size());

            for(unsigned int i = 0; i < t_idc.size(); i++)   t0[t_idc[i]] += dt[i];                                              //温度の更新
            qtsum_0 = qt_sum(t0, ts);
            for(unsigned int i = 0; i < t_idc.size(); i++)   rmse += pow(qtsum_0[t_idc[i]], 2.0) / t_idc.size(); 

            if(i_tn_ac != -1){
                switch(tn[i_tn_ac].ac_mode[ts]){
                    case AC_AUTO:       if(abs(t0[tn[i_tn_ac].i1] - tn[i_tn_ac].pre_tmp[ts]) > sts.thrm_err){
                                            t0[tn[i_tn_ac].i1] = tn[i_tn_ac].pre_tmp[ts];
                                            rmse = 999;
                                        }
                                        break;
                    case AC_HEATING:    if(t0[tn[i_tn_ac].i1] < tn[i_tn_ac].pre_tmp[ts]){
                                            t0[tn[i_tn_ac].i1] = tn[i_tn_ac].pre_tmp[ts];
                                            rmse = 999;
                                        }                   
                                        break;
                    case AC_COOLING:    if(t0[tn[i_tn_ac].i1] > tn[i_tn_ac].pre_tmp[ts]){
                                            t0[tn[i_tn_ac].i1] = tn[i_tn_ac].pre_tmp[ts];
                                            rmse = 999;
                                        }
                                        break;
                }       
            }
            LOG_PRINT(itr << ": ts = " << ts << ": rmse = " << sqrt(rmse) << endl);
            itr++;
        }while(sts.thrm_err < sqrt(rmse));

        for(unsigned int i = 0; i < t_idc.size(); i++)   sn[t_idc[i]].t[ts] = t0[t_idc[i]];                                          //温度の計算
        return rmse;
    }

    void calc_qv(long ts1, long ts2){
        //LOG_PRINT("ts1 = " << ts1 << " <<<<---- " << "ts2 = " <<  ts2 << endl);
        for(unsigned int i = 0; i < vn.size(); i++){    
            double rgh1 = get_rho(sn[vn[i].i1].t[ts2]) * G * vn[i].h1[ts2];
            double rgh2 = get_rho(sn[vn[i].i2].t[ts2]) * G * vn[i].h2[ts2];
            vn[i].qv[ts1] = vn[i].get_qv((sn[vn[i].i1].p[ts2] - rgh1) - (sn[vn[i].i2].p[ts2] - rgh2), ts2);                 //風量の計算
        }
    }

    void calc_qt(long ts1, long ts2){
        //LOG_PRINT("ts1 = " << ts1 << " <<<<---- " << "ts2 = " << ts2 << endl);
        for(unsigned int i = 0; i < tn.size(); i++){    
            switch(tn[i].tn_type){
                case TN_SIMPLE: 
                case TN_GROUND:         tn[i].qt[ts1] =   tn[i].get_qt(sn[tn[i].i1].t[ts2] - sn[tn[i].i2].t[ts2], ts2);
                                        break;
                case TN_SOLAR:          tn[i].qt[ts1] = - tn[i].get_qt_s(sn[tn[i].i2].h_sr[ts2], ts2);
                                        break;
                case TN_HEATER:         tn[i].qt[ts1] = - tn[i].get_qt_h(sn[tn[i].i2].h_inp[ts2], ts2);
                                        break;
            }
        }
        for(unsigned int i = 0; i < vn.size(); i++)
            vn[i].qt[ts1] = vn[i].get_qt(sn[vn[i].i1].t[ts2] - sn[vn[i].i2].t[ts2], ts2);                                   //熱量の計算
    }

    int calc(){
        vector<double>  pre_p(sn.size(), 0.0), pre_t(sn.size(), 0.0);
        double          delta_p, delta_t; 

        LOG_PRINT("******************************************************************************Start calc!" << endl);     

        LOG_PRINT("length:     " << sts.length << endl);
        LOG_PRINT("sts-t_step: " << sts.t_step << endl);
        LOG_PRINT("sts-solve:  " << sts.solve << endl);
        LOG_PRINT("step_p:     " << sts.step_p << endl);
        LOG_PRINT("vent_err:   " << sts.vent_err << endl);
        LOG_PRINT("step_t:     " << sts.step_t << endl);
        LOG_PRINT("thrm_err:   " << sts.thrm_err << endl);
        LOG_PRINT("sor_ratio:  " << sts.sor_ratio << endl);
        LOG_PRINT("sor_err:    " << sts.sor_err << endl);
        LOG_CONTENTS(endl);

        for(unsigned int i = 0; i < sn.size(); i++){
            LOG_PRINT("sn[" << i << "] ="); 
            LOG_CONTENTS(get<0>(sn[i].flag) << "," << get<1>(sn[i].flag)  << "," << get<2>(sn[i].flag) << ") ");
            LOG_CONTENTS("i=" << sn[i].i << " ,s_i=" << sn[i].s_i << ", p[0]=" << sn[i].p[0] << ", c[0]" << sn[i].c[0] << ", t[0]" << sn[i].t[0]);
            if(sn[i].m.size()     != 0) LOG_CONTENTS(", m[0]"       << sn[i].m[0]);
            if(sn[i].h_sr.size()  != 0) LOG_CONTENTS(", h_sr[0]="   << sn[i].h_sr[0]);
            if(sn[i].h_inp.size() != 0) LOG_CONTENTS(", h_inp[0]="  << sn[i].h_inp[0]);
            if(sn[i].v.size()     != 0) LOG_CONTENTS(", v[0]="      << sn[i].v[0]);
            if(sn[i].beta.size()  != 0) LOG_CONTENTS(", beta[0]="   << sn[i].beta[0]);
            LOG_CONTENTS(endl);
        }
        for(unsigned int i = 0; i < vn.size(); i++){
            LOG_PRINT("vn[" << i << "] = " << vn[i].vn_type << " (" << vn[i].i1 << "," << vn[i].i2  << ") " << vn[i].h1[0] << " - " << vn[i].h2[0]);
            LOG_CONTENTS(",qv[0]=" << vn[i].qv[0] << ",qt[0]=" << vn[i].qt[0]); 
            if(vn[i].alpha.size() != 0) LOG_CONTENTS(", alpha[0]=" << vn[i].alpha[0]);
            if(vn[i].area.size()  != 0) LOG_CONTENTS(", area[0]="  << vn[i].area[0]);
            if(vn[i].a.size()     != 0) LOG_CONTENTS(", a[0]="     << vn[i].a[0]);
            if(vn[i].n.size()     != 0) LOG_CONTENTS(", n[0]="     << vn[i].n[0]);
            if(vn[i].eta.size()   != 0) LOG_CONTENTS(", eta[0]="   << vn[i].eta[0]);
            if(vn[i].q_max.size() != 0) LOG_CONTENTS(", p_max[0]=" << vn[i].q_max[0]);
            if(vn[i].p_max.size() != 0) LOG_CONTENTS(", q_max[0]=" << vn[i].p_max[0]);
            if(vn[i].q1.size()    != 0) LOG_CONTENTS(", q1[0]="    << vn[i].q1[0]);
            if(vn[i].p1.size()    != 0) LOG_CONTENTS(", p1[0]="    << vn[i].p1[0]);
            LOG_CONTENTS(endl);
        }

        for(unsigned int i = 0; i < tn.size(); i++){
            LOG_PRINT("tn[" << i << "] = " << tn[i].tn_type << " (" << tn[i].i1 << "," << tn[i].i2 << ")");
            LOG_CONTENTS(",qt[0]=" << tn[i].qt[0]);
            if(tn[i].cdtc.size()      != 0) LOG_CONTENTS(", cdtc[0]="      << tn[i].cdtc[0]);
            if(tn[i].ms.size()        != 0) LOG_CONTENTS(", ms[0]="        << tn[i].ms[0]);
            if(tn[i].area.size()      != 0) LOG_CONTENTS(", area[0]="      << tn[i].area[0]);
            if(tn[i].rg.size()        != 0) LOG_CONTENTS(", rg[0]="        << tn[i].rg[0]);
            if(tn[i].cof_r.size()     != 0) LOG_CONTENTS(", cof_r[0]="     << tn[i].cof_r[0]);
            if(tn[i].cof_phi.size()   != 0) LOG_CONTENTS(", cof_phi[0]="   << tn[i].cof_phi[0]);
            if(tn[i].t_dash_gs.size() != 0) LOG_CONTENTS(", t_dash_gs[0]=" << tn[i].t_dash_gs[0]);
            if(tn[i].aircon_on.size() != 0) LOG_CONTENTS(", aircon_on[0]=" << tn[i].aircon_on[0]);
            if(tn[i].ac_mode.size()   != 0) LOG_CONTENTS(", alpha[0]="     << tn[i].ac_mode[0]);
            if(tn[i].pre_tmp.size()   != 0) LOG_CONTENTS(", alpha[0]="     << tn[i].pre_tmp[0]);
            LOG_CONTENTS(endl);
        }

        LOG_PRINT("i_vn_ac = " << i_vn_ac << endl);
        LOG_PRINT("i_tn_ac = " << i_tn_ac << endl);
        LOG_PRINT("t_idc.size() = " << t_idc.size() << endl);
        LOG_PRINT("c_idc.size() = " << c_idc.size() << endl);
        LOG_PRINT("v_idc.size() = " << v_idc.size() << endl); 
        LOG_CONTENTS(endl);

        for(long ts = 1; ts < sts.length; ts++){
            if(v_idc.size() > 0){
                for(unsigned int i = 0; i < sn.size(); i++)
                    if(get<0>(sn[i].flag) == SN_CALC)   sn[i].p[ts] = sn[i].p[ts - 1];
                calc_qv(ts, ts - 1);
            }

            if(t_idc.size() > 0){    
                for(unsigned int i = 0; i < sn.size(); i++){
                    if(get<2>(sn[i].flag) == SN_CALC)   sn[i].t[ts] = sn[i].t[ts - 1];
                    if(get<2>(sn[i].flag) == SN_DLY)    sn[i].t[ts] = sn[sn[i].s_i].t[ts - 1];
                }
                calc_qt(ts, ts - 1);
                for(unsigned int i = 0; i < tn.size(); i++)
                    if(tn[i].tn_type == TN_GROUND)      tn[i].refresh(sn[tn[i].i1].t[ts], sn[tn[i].i2].t[ts], ts);   //地盤のリフレッシュ
            }

            do{
                delta_p = 0.0;
                delta_t = 0.0;

                for(unsigned int i = 0; i < sn.size(); i++){
                    if(get<0>(sn[i].flag) == SN_CALC)   pre_p[i] = sn[i].p[ts];
                    if(get<2>(sn[i].flag) == SN_CALC)   pre_t[i] = sn[i].t[ts];
                }
                calc_vent(ts);
                calc_thrm(ts);
                
                for(unsigned int i = 0; i < sn.size(); i++){
                    if(get<0>(sn[i].flag) == SN_CALC)   delta_p += pow(sn[i].p[ts] - pre_p[i], 2.0) / v_idc.size();
                    if(get<2>(sn[i].flag) == SN_CALC)   delta_t += pow(sn[i].t[ts] - pre_t[i], 2.0) / t_idc.size();
                }
                LOG_PRINT("ts = " << ts << " delta_p = " << delta_p << ", delta_t = " << delta_t << " >>>>> " << sqrt((delta_p + delta_t)/ 2) << endl);

            }while(sts.conv_err < sqrt((delta_p + delta_t)/ 2));
            
            if(v_idc.size() > 0){
                calc_qv(ts, ts);
                LOG_CONTENTS("p ");
                for(unsigned int i = 0; i < sn.size(); i++)    LOG_CONTENTS(sn[i].i << ": " << sn[i].p[ts] << "Pa, ");
                LOG_CONTENTS(endl << "qv ");
                for(unsigned int i = 0; i < vn.size(); i++)    LOG_CONTENTS(vn[i].i << ": " << vn[i].qv[ts] * 3600 << "m3/h, ");
                LOG_CONTENTS(endl);    
            }

            if(t_idc.size()> 0){
                calc_qt(ts, ts);
                LOG_CONTENTS("t ");
                for(unsigned int i = 0; i < sn.size(); i++)    LOG_CONTENTS(sn[i].i << ": " << sn[i].t[ts] << "deg, ");
                LOG_CONTENTS(endl << "qt1 ");
                for(unsigned int i = 0; i < vn.size(); i++)    LOG_CONTENTS(vn[i].i << ": " << vn[i].qt[ts] << "W, ");
                LOG_CONTENTS(endl << "qt2 ");
                for(unsigned int i = 0; i < tn.size(); i++)    LOG_CONTENTS(tn[i].i << ": " << tn[i].qt[ts] << "W, ");
                LOG_CONTENTS(endl);
            }

            if(c_idc.size() > 0){
                for(unsigned int i = 0; i < c_idc.size(); i++){
                    sn[c_idc[i]].c[ts] =  sn[c_idc[i]].c[ts - 1] * exp(-sn[c_idc[i]].beta[ts] * sts.t_step);
                    sn[c_idc[i]].c[ts] += sn[c_idc[i]].m[ts] / (sn[c_idc[i]].beta[ts] * sn[c_idc[i]].v[ts]) * (1 - exp(-sn[c_idc[i]].beta[ts] * sts.t_step));
                }
                for(unsigned int i = 0; i < vn.size(); i++){
                    if(vn[i].qv[ts] > 0 && get<1>(sn[vn[i].i2].flag) == SN_CALC)    
                        sn[vn[i].i2].c[ts] += vn[i].qv[ts] * (sn[vn[i].i1].c[ts] * (1 - vn[i].eta[ts]) - sn[vn[i].i2].c[ts]) * sts.t_step / sn[vn[i].i2].v[ts];
                    if(vn[i].qv[ts] < 0 && get<1>(sn[vn[i].i1].flag) == SN_CALC)    
                        sn[vn[i].i1].c[ts] += vn[i].qv[ts] * (sn[vn[i].i2].c[ts] * (1 - vn[i].eta[ts]) - sn[vn[i].i1].c[ts]) * sts.t_step / sn[vn[i].i1].v[ts];
                }
                for(unsigned int i = 0; i < sn.size(); i++)    LOG_CONTENTS("c" << sn[i].i << ": " << sn[i].c[ts] << "num/L, ");
                LOG_CONTENTS(endl);
            }
            LOG_CONTENTS(endl);
        }
        LOG_PRINT("******************************************************************************Finish calc!" << endl << endl);
        return 0;
    }

    tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, 
          vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> result(){

        vector<vector<double>>   r_p, r_c, r_t, r_qv, r_qt1, r_qt2;
    
        for(Node n: sn){
            r_p.push_back(n.p);
            r_c.push_back(n.c);
            r_t.push_back(n.t);
        }
        for(Vent_Net nt: vn){
            r_qv.push_back(nt.qv);
            r_qt1.push_back(nt.qt);
        }
        for(Thrm_Net nt: tn)   
            r_qt2.push_back(nt.qt);

        return forward_as_tuple(r_p, r_c, r_t, r_qv, r_qt1, r_qt2);
    }

};

using namespace std;

class Thrm_Net{
public:
    int i, i1, i2, tn_type;                                                                 //入口出口のノード番号、熱回路網の種類
    
    vector<double>      cdtc;                                                               //コンダクタンス
    vector<double>      ms;                                                                 //日射取得率
    vector<double>      area, rg;                                                           //面積、表面熱抵抗
    double              phi_0;                                                              //応答係数                                                        
    vector<double>      cof_r, cof_phi, t_dash_gs;
    vector<double>      qt;                                                                 //熱流
    vector<int>         aircon_on, ac_mode;                                                 //エアコンのON/OFF、エアコン運転モード
    vector<double>      pre_tmp;                                                            //エアコン設定温度

    Thrm_Net(long length, int i_, int i1_ , int i2_, int tn_type_){
        i       = i_;
        i1      = i1_;
        i2      = i2_;
        tn_type = tn_type_;
        qt.assign(length, 0.0);
        t_dash_gs.assign(10, 0.0);
        aircon_on.assign(length, 0);                                                        //熱量を0に初期化        
    }              

    double get_qt(double dt, long ts){
        double sum_t_dash_gs = 0.0;                                                         //項別成分の合計
        switch(tn_type){
            case TN_SIMPLE:    
                return cdtc[ts] * dt;                                                       //単純熱回路網の場合
            case TN_GROUND:                                                                 //地盤の場合  
                for(int i = 0; i < 10; i++)     sum_t_dash_gs += t_dash_gs[i];              //合計の計算
                return area[ts] / (rg[ts] + phi_0) * (dt - sum_t_dash_gs);                  //熱流の計算
            default:
                return 0.0;
        }
    }

    double get_qt_s(double h_sr, long ts)  {return ms[ts] * h_sr;}                          //日射取得
    double get_qt_h(double h_inp, long ts) {return h_inp;}                                  //発熱
    
    void refresh(double t_uf_t_1, double t_g_t_1, long ts){
        double sum_t_dash_gs = 0.0;                                                                 //項別成分の合計
        for(int i = 0; i < 10; i++)         sum_t_dash_gs += t_dash_gs[i];                          //合計の計算
        double t_gs = (phi_0 * t_uf_t_1 + rg[ts] * (sum_t_dash_gs + t_g_t_1)) / (rg[ts] + phi_0);   //表面温度の計算
        for(int i = 0; i < 10; i++)
            t_dash_gs[i] = cof_phi[i] * (t_uf_t_1 - t_gs) / rg[ts] + cof_r[i] * t_dash_gs[i];       //吸熱応答の項別成分の計算
    }
};
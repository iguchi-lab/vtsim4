using namespace std;

class Node{
public:
    int i, s_i;
    tuple<int, int, int> flag;                              //換気、濃度、熱計算フラグ                
    vector<double> p, c, t;                                 //圧力、濃度、温度
    vector<double> m;                                       //発生量
    vector<double> h_sr, h_inp;                             //日射量、発熱量
    vector<double> v;                                       //気積     
    vector<double> beta;                                    //沈着率

    Node(long length, int i_, tuple<int, int, int> flag_){
        i = i_;
        flag = flag_;
        p.assign(length, 0.0);
        c.assign(length, 0.0);
        t.assign(length, 20.0);
        if(get<1>(flag) == SN_CALC)   m.assign(length, 0.0);
        if(get<2>(flag) == SN_CALC)   h_inp.assign(length, 0.0);
    }

};
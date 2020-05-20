#ifndef TI_OBJECT_H
#define TI_OBJECT_H

struct Data {
    Data() = default;
    Data(int a, double b) :
    _i(a), _d(b) {};
    int _i;
    double _d;
};

struct CTI_Object {
    CTI_Object();
    void Try();
    void Try2();
    void Try3();
    Data ** Try4();
};



#endif //TI_OBJECT_H
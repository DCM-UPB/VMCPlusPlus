#include "vmc/TwoBodyPseudoPotential.hpp"



double TwoBodyPseudoPotential::u(const double * r1, const double * r2){
    return ur(_metric->dist(r1, r2));
}


void TwoBodyPseudoPotential::computeAllDerivatives(const double * r1, const double * r2){
    const double ud1 = urD1(_metric->dist(r1, r2));
    const double ud2 = urD2(_metric->dist(r1, r2));

    _metric->distD1(r1, r2, _foo);
    _metric->distD2(r1, r2, _foo2);

    if (_flag_vd1) urVD1(_metric->dist(r1, r2), _vfoo);
    if (_flag_d1vd1 || _flag_d2vd1) urD1VD1(_metric->dist(r1, r2), _vfoo1);
    if (_flag_d2vd1) urD2VD1(_metric->dist(r1, r2), _vfoo2);

    for (int i=0; i<_ndim2; ++i){
        _d1[i] = _foo[i] * ud1;
    }

    for (int i=0; i<_ndim2; ++i){
        _d2[i] = _foo2[i] * ud1 + _foo[i] * _foo[i] * ud2;
    }

    if (_flag_vd1){
        for (int i=0; i<_nvp; ++i){
            _vd1[i] = _vfoo[i];
        }
    }

    if (_flag_d1vd1){
        for (int i=0; i<_ndim2; ++i){
            for (int j=0; j<_nvp; ++j){
                _d1vd1[i][j] = _foo[i] * _vfoo1[j];
            }
        }
    }

    if (_flag_d2vd1){
        for (int i=0; i<_ndim2; ++i){
            for (int j=0; j<_nvp; ++j){
                _d2vd1[i][j] = _foo2[i] * _vfoo1[j] + _foo[i] * _foo[i] * _vfoo2[j];
            }
        }
    }
}


TwoBodyPseudoPotential::TwoBodyPseudoPotential(Metric * metric, const int &nvp, bool flag_vd1, bool flag_d1vd1, bool flag_d2vd1){
    _metric = metric;
    _ndim2 = 2*_metric->getNSpaceDim();
    _nvp = nvp;

    _flag_vd1 = flag_vd1;
    _flag_d1vd1 = flag_d1vd1;
    _flag_d2vd1 = flag_d2vd1;

    _d1 = new double[_ndim2];
    _d2 = new double[_ndim2];

    _vd1 = 0;
    if (flag_vd1){
        _vd1 = new double[_nvp];
    }

    _d1vd1 = 0;
    if (flag_d1vd1){
        _d1vd1 = new double*[_ndim2];
        for (int id1=0; id1<_ndim2; ++id1) _d1vd1[id1] = new double[_nvp];
    }

    _d2vd1 = 0;
    if (flag_d2vd1){
        _d2vd1 = new double*[_ndim2];
        for (int id2=0; id2<_ndim2; ++id2) _d2vd1[id2] = new double[_nvp];
    }

    _foo = new double[_ndim2];
    _foo2 = new double[_ndim2];
    _vfoo = 0;
    if (flag_vd1 || flag_d1vd1 || flag_d2vd1) _vfoo = new double[_nvp];
    _vfoo1 = 0;
    if (flag_d1vd1 || flag_d2vd1) _vfoo1 = new double[_nvp];
    _vfoo2 = 0;
    if (flag_d2vd1) _vfoo2 = new double[_nvp];
}


TwoBodyPseudoPotential::~TwoBodyPseudoPotential(){
    _metric = 0;

    delete[] _d1; _d1=0;
    delete[] _d2; _d2=0;
    if (_vd1 != 0){
        delete[] _vd1; _vd1 = 0;
    }
    if (_d1vd1 != 0){
        for (int i=0; i<_ndim2; ++i) delete[] _d1vd1[i];
        delete[] _d1vd1; _d1vd1 = 0;
    }
    if (_d2vd1 != 0){
        for (int i=0; i<_ndim2; ++i) delete[] _d2vd1[i];
        delete[] _d2vd1; _d2vd1 = 0;
    }

    delete[] _foo; _foo = 0;
    delete[] _foo2; _foo2 = 0;
    if (_vfoo){
        delete[] _vfoo; _vfoo = 0;
    }
    if (_vfoo1){
        delete[] _vfoo1; _vfoo1 = 0;
    }
    if (_vfoo2){
        delete[] _vfoo2; _vfoo2 = 0;
    }
}

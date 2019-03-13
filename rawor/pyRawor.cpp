#include <vector>
#include <cmath>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bp = boost::python;
namespace np = boost::python::numpy;

class rawor{
    double nbar_ran, Delta_r, r_max, r_min, V_box;
    int N_parts, N_rans, N_shells;
    std::vector<double> rs, w, x;
    
    void initVectors();
    
    void swapIfGreater(double &a, double &b);
    
    double sphereOverlapVolume(double d, double R, double r);
    
    double crossSectionVolume(double r1, double r2, double r3);
    
    int getPermutations(double r1, double r2, double r3);
    
    double sphericalShellVolume(double r);
    
    double nbarData(unsigned long *DD, double r, double r1);
    
    double gaussQuadCrossSection(double r1, double r2, double r3);
    
    double gaussQuadCrossSectionDDR(unsigned long *DD, double r1, double r2, double r3);
    
    public:
        rawor(int numParticles, int numRandoms, int numShells, double VolBox, double rMax, double rMin = 0);
        
        void setNumParts(int numParticles);
        
        void setNumRans(int numRandoms);
        
        void setNumShells(int numShells);
        
        void setRMax(double rMax);
        
        void setRMin(double rMin);
        
        void setVBox(double VBox);
        
        int getNumParts();
        
        int getNumRans();
        
        int getNumShells();
        
        double getRMax();
        
        double getRMin();
        
        double getVBox();
        
        np::ndarray getRRR();
        
        np::ndarray getDRR();
        
        np::ndarray getDDR(np::ndarray const &dd);
};

void rawor::initVectors() {
    for (int i = 0; i < rawor::N_shells; ++i) {
        rawor::rs.push_back(rawor::r_min + (i + 0.5)*rawor::Delta_r);
    }
    
    rawor::w = {0.8888888888888888, 0.5555555555555556, 0.5555555555555556};
    
    rawor::x = {0.0000000000000000, -0.7745966692414834, 0.7745966692414834};
}

void rawor::swapIfGreater(double &a, double &b) {
    if (a > b) {
        double temp = a;
        a = b;
        b = temp;
    }
}

double rawor::sphereOverlapVolume(double d, double R, double r) {
    double V = 0;
    swapIfGreater(r, R);
    if (d < R + r) {
        if (d > R - r) {
            V = (M_PI*(R + r - d)*(R + r - d)*(d*d + 2.0*d*r - 3.0*r*r + 2.0*d*R + 6.0*r*R - 3.0*R*R))/(12.0*d);
        } else {
            V = (4.0*M_PI/3.0)*r*r*r;
        }
    }
    return V;
}

double rawor::crossSectionVolume(double r1, double r2, double r3) {
    double V_oo = sphereOverlapVolume(r1, r3 + 0.5*rawor::Delta_r, r2 + 0.5*rawor::Delta_r);
    double V_oi = sphereOverlapVolume(r1, r3 + 0.5*rawor::Delta_r, r2 - 0.5*rawor::Delta_r);
    double V_io = sphereOverlapVolume(r1, r3 - 0.5*rawor::Delta_r, r2 + 0.5*rawor::Delta_r);
    double V_ii = sphereOverlapVolume(r1, r3 - 0.5*rawor::Delta_r, r2 - 0.5*rawor::Delta_r);
    
    return V_oo - V_oi - V_io + V_ii;
}

int rawor::getPermutations(double r1, double r2, double r3) {
    int perm = 1;
    if (r1 != r2 && r1 != r3 && r2 != r3) {
        perm = 6;
    } else if ((r1 == r2 && r1 != r3) || (r1 == r3 && r1 != r2) || (r2 == r3 && r2 != r1)) {
        perm = 3;
    }
    return perm;
}

double rawor::sphericalShellVolume(double r) {
    double r_o = r + 0.5*rawor::Delta_r;
    double r_i = r - 0.5*rawor::Delta_r;
    return 4.0*M_PI*(r_o*r_o*r_o - r_i*r_i*r_i)/3.0;
}

double rawor::nbarData(unsigned long *DD, double r, double r1) {
    int bin = r/rawor::Delta_r;
    double nbar = DD[bin]/(rawor::N_parts*sphericalShellVolume(r1));
    int num_bins = rawor::N_shells;
    if (r <= (bin + 0.5)*rawor::Delta_r) {
        if (bin != 0) {
            double n1 = DD[bin]/(rawor::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin - 1]/(rawor::N_parts*sphericalShellVolume(r1 - rawor::Delta_r));
            double b = n1 - ((n1 - n2)/rawor::Delta_r)*r1;
            nbar = ((n1 - n2)/rawor::Delta_r)*r + b;
        } else {
            double n1 = DD[bin]/(rawor::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin + 1]/(rawor::N_parts*sphericalShellVolume(r1 + rawor::Delta_r));
            double b = n1 - ((n2 - n1)/rawor::Delta_r)*r1;
            nbar = ((n2 - n1)/rawor::Delta_r)*r + b;
        }
    } else {
        if (bin != num_bins - 1) {
            double n1 = DD[bin]/(rawor::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin + 1]/(rawor::N_parts*sphericalShellVolume(r1 + rawor::Delta_r));
            double b = n1 - ((n2 - n1)/rawor::Delta_r)*r1;
            nbar = ((n2 - n1)/rawor::Delta_r)*r + b;
        } else {
            double n1 = DD[bin]/(rawor::N_parts*sphericalShellVolume(r1));
            double n2 = DD[bin - 1]/(rawor::N_parts*sphericalShellVolume(r1 - rawor::Delta_r));
            double b = n1 - ((n1 - n2)/rawor::Delta_r)*r1;
            nbar = ((n1 - n2)/rawor::Delta_r)*r + b;
        }
    }
    return nbar;
}

double rawor::gaussQuadCrossSection(double r1, double r2, double r3) {
    double result = 0.0;
    for (int i = 0; i < rawor::w.size(); ++i) {
        double r_1 = r1 + 0.5*rawor::Delta_r*rawor::x[i];
        result += 0.5*rawor::Delta_r*rawor::w[i]*crossSectionVolume(r_1, r2, r3)*r_1*r_1;
    }
    return result;
}

double rawor::gaussQuadCrossSectionDDR(unsigned long *DD, double r1, double r2, double r3) {
    double result = 0.0;
    for (int i = 0; i < rawor::w.size(); ++i) {
        double r_1 = r1 + 0.5*rawor::Delta_r*rawor::x[i];
        double nbar = nbarData(DD, r_1, r1);
        result += 0.5*rawor::Delta_r*rawor::w[i]*crossSectionVolume(r_1, r2, r3)*r_1*r_1*nbar;
    }
    return result;
}

rawor::rawor(int numParticles, int numRandoms, int numShells, double VolBox, double rMax, double rMin) {
    rawor::N_parts = numParticles;
    rawor::N_rans = numRandoms;
    rawor::N_shells = numShells;
    rawor::r_max = rMax;
    rawor::r_min = rMin;
    rawor::V_box = VolBox;
    rawor::Delta_r = (rMax - rMin)/numShells;
    rawor::nbar_ran = numRandoms/VolBox;
    rawor::initVectors();
}

void rawor::setNumParts(int numParticles) {
    rawor::N_parts = numParticles;
}

void rawor::setNumRans(int numRandoms) {
    rawor::N_rans = numRandoms;
    rawor::nbar_ran = numRandoms/rawor::V_box;
}

void rawor::setNumShells(int numShells) {
    rawor::N_shells = numShells;
    rawor::Delta_r = (rawor::r_max - rawor::r_min)/rawor::N_shells;
}

void rawor::setRMax(double rMax) {
    rawor::r_max = rMax;
    rawor::Delta_r = (rawor::r_max - rawor::r_min)/rawor::N_shells;
}

void rawor::setRMin(double rMin) {
    rawor::r_min = rMin;
    rawor::Delta_r = (rawor::r_max - rawor::r_min)/rawor::N_shells;
}

void rawor::setVBox(double VBox) {
    rawor::V_box = VBox;
    rawor::nbar_ran = rawor::N_rans/rawor::V_box;
}

int rawor::getNumParts() {
    return rawor::N_parts;
}

int rawor::getNumRans() {
    return rawor::N_rans;
}

int rawor::getNumShells() {
    return rawor::N_shells;
}

double rawor::getVBox() {
    return rawor::V_box;
}

double rawor::getRMin() {
    return rawor::r_min;
}

double rawor::getRMax() {
    return rawor::r_max;
}

np::ndarray rawor::getRRR() {
    std::vector<int> N;
    for (int i = 0; i < rawor::N_shells; ++i) {
        for (int j = i; j < rawor::N_shells; ++j) {
            for (int k = j; k < rawor::N_shells; ++k) {
                if (rawor::rs[k] <= rawor::rs[i] + rawor::rs[j]) {
                    int index = k + rawor::N_shells*(j + rawor::N_shells*i);
                    double V = rawor::gaussQuadCrossSection(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    int n_perm = rawor::getPermutations(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    N.push_back(int(4.0*M_PI*n_perm*rawor::nbar_ran*rawor::nbar_ran*V*rawor::N_rans));
                }
            }
        }
    }
    np::dtype dt = np::dtype::get_builtin<int>();
    np::ndarray n = np::zeros(bp::make_tuple(N.size()), dt);
    std::copy(N.begin(), N.end(), reinterpret_cast<int*>(n.get_data()));
    return n;
}

np::ndarray rawor::getDRR() {
    std::vector<int> N;
    for (int i = 0; i < rawor::N_shells; ++i) {
        for (int j = i; j < rawor::N_shells; ++j) {
            for (int k = j; k < rawor::N_shells; ++k) {
                if (rawor::rs[k] <= rawor::rs[i] + rawor::rs[j]) {
                    int index = k + rawor::N_shells*(j + rawor::N_shells*i);
                    double V = rawor::gaussQuadCrossSection(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    int n_perm = rawor::getPermutations(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    N.push_back(int(4.0*M_PI*n_perm*rawor::nbar_ran*rawor::nbar_ran*V*rawor::N_parts));
                }
            }
        }
    }
    np::dtype dt = np::dtype::get_builtin<int>();
    np::ndarray n = np::zeros(bp::make_tuple(N.size()), dt);
    std::copy(N.begin(), N.end(), reinterpret_cast<int*>(n.get_data()));
    return n;
}

np::ndarray rawor::getDDR(np::ndarray const &dd) {
    unsigned long *DD = reinterpret_cast<unsigned long *>(dd.get_data());
    std::vector<int> N;
    for (int i = 0; i < rawor::N_shells; ++i) {
        double r1 = rawor::rs[i];
        for (int j = i; j < rawor::N_shells; ++j) {
            double r2 = rawor::rs[j];
            for (int k = j; k < rawor::N_shells; ++k) {
                double r3 = rawor::rs[k];
                if (rawor::rs[k] <= rawor::rs[i] + rawor::rs[j]) {
                   int index = k + rawor::N_shells*(j + rawor::N_shells*i);
                   double V = rawor::gaussQuadCrossSectionDDR(DD, r1, r2, r3);
                   double N_temp = 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                   if (r1 != r2 && r1 != r3 && r2 != r3) {
                       V = rawor::gaussQuadCrossSectionDDR(DD, r2, r3, r1);
                       N_temp += 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                       V = rawor::gaussQuadCrossSectionDDR(DD, r3, r1, r2);
                       N_temp += 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                       V = rawor::gaussQuadCrossSectionDDR(DD, r1, r3, r2);
                       N_temp += 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                       V = rawor::gaussQuadCrossSectionDDR(DD, r2, r1, r3);
                       N_temp += 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                       V = rawor::gaussQuadCrossSectionDDR(DD, r3, r2, r1);
                       N_temp += 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                   } else if ((r1 == r2 && r1 != r3) || (r1 == r3 && r1 != r2) || (r2 == r3 && r2 != r1)) {
                       V = rawor::gaussQuadCrossSectionDDR(DD, r2, r3, r1);
                       N_temp += 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                       V = rawor::gaussQuadCrossSectionDDR(DD, r3, r1, r2);
                       N_temp += 4.0*M_PI*rawor::nbar_ran*V*rawor::N_parts;
                   }
                   N.push_back(int(floor(N_temp + 0.5)));
                }
            }
        }
    }
    np::dtype dt = np::dtype::get_builtin<int>();
    np::ndarray n = np::zeros(bp::make_tuple(N.size()), dt);
    std::copy(N.begin(), N.end(), reinterpret_cast<int*>(n.get_data()));
    delete[] DD;
    return n;
}

BOOST_PYTHON_MODULE(rawor) {
    np::initialize();
    using namespace boost::python;
    
    class_<rawor>("rawor", init<int, int, int, double, double, double>())
        .def("set_num_parts", &rawor::setNumParts)
        .def("set_num_rans", &rawor::setNumRans)
        .def("set_num_shells", &rawor::setNumShells)
        .def("set_r_max", &rawor::setRMax)
        .def("set_r_min", &rawor::setRMin)
        .def("set_V_box", &rawor::setVBox)
        .def("get_num_parts", &rawor::getNumParts)
        .def("get_num_rans", &rawor::getNumRans)
        .def("get_num_shells", &rawor::getNumShells)
        .def("get_V_box", &rawor::getVBox)
        .def("get_r_min", &rawor::getRMin)
        .def("get_r_max", &rawor::getRMax)
        .def("get_RRR", &rawor::getRRR)
        .def("get_DRR", &rawor::getDRR)
        .def("get_DDR", &rawor::getDDR)
    ;
}

#ifndef _RAWOR_HPP_
#define _RAWOR_HPP_

#include <vector>
#include <cmath>

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
    
    double nbarData(std::vector<int> &DD, double r, double r1);
    
    double gaussQuadCrossSection(double r1, double r2, double r3);
    
    double gaussQuadCrossSectionDDR(std::vector<int> &DD, double r1, double r2, double r3);
    
    public:
        rawor(int numParticles, int numRandoms, int numShells, double VolBox, double rMax, double rMin = 0);
        
        void updateNumParts(int numParticles);
        
        void updateNumRans(int numRandoms);
        
        std::vector<int> getRRR();
        
        std::vector<int> getDRR();
        
        std::vector<int> getDDR(std::vector<int> &DD);
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

double rawor::nbarData(std::vector<int> &DD, double r, double r1) {
    int bin = r/rawor::Delta_r;
    double nbar = DD[bin]/(rawor::N_parts*sphericalShellVolume(r1));
    int num_bins = DD.size();
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

double rawor::gaussQuadCrossSectionDDR(std::vector<int> &DD, double r1, double r2, double r3) {
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

void rawor::updateNumParts(int numParticles) {
    rawor::N_parts = numParticles;
}

void rawor::updateNumRans(int numRandoms) {
    rawor::N_rans = numRandoms;
    rawor::nbar_ran = numRandoms/rawor::V_box;
}

std::vector<int> rawor::getRRR() {
    std::vector<int> N(rawor::N_shells*rawor::N_shells*rawor::N_shells);
    for (int i = 0; i < rawor::N_shells; ++i) {
        for (int j = i; j < rawor::N_shells; ++j) {
            for (int k = j; k < rawor::N_shells; ++k) {
                if (rawor::rs[k] <= rawor::rs[i] + rawor::rs[j]) {
                    int index = k + rawor::N_shells*(j + rawor::N_shells*i);
                    double V = rawor::gaussQuadCrossSection(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    int n_perm = rawor::getPermutations(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    N[index] = int(4.0*M_PI*n_perm*rawor::nbar_ran*rawor::nbar_ran*V*rawor::N_rans);
                }
            }
        }
    }
    return N;
}

std::vector<int> rawor::getDRR() {
    std::vector<int> N(rawor::N_shells*rawor::N_shells*rawor::N_shells);
    for (int i = 0; i < rawor::N_shells; ++i) {
        for (int j = i; j < rawor::N_shells; ++j) {
            for (int k = j; k < rawor::N_shells; ++k) {
                if (rawor::rs[k] <= rawor::rs[i] + rawor::rs[j]) {
                    int index = k + rawor::N_shells*(j + rawor::N_shells*i);
                    double V = rawor::gaussQuadCrossSection(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    int n_perm = rawor::getPermutations(rawor::rs[i], rawor::rs[j], rawor::rs[k]);
                    N[index] = int(4.0*M_PI*n_perm*rawor::nbar_ran*rawor::nbar_ran*V*rawor::N_parts);
                }
            }
        }
    }
    return N;
}

std::vector<int> rawor::getDDR(std::vector<int> &DD) {
    std::vector<int> N(rawor::N_shells*rawor::N_shells*rawor::N_shells);
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
                   N[index] = int(floor(N_temp + 0.5));
                }
            }
        }
    }
    return N;
}

#endif

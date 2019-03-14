/* rawor.hpp (RAndoms WithOut Randoms)
 * 
 * MIT License
 * 
 * Copyright (c) 2019 David W. Pearson (DWP)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * This header only library will compute predict the RRR, DRR and DDR counts
 * for use in the calculation of the three point correlation function (3PCF) of
 * points in periodic boxes.
 */

#ifndef _RAWOR_HPP_
#define _RAWOR_HPP_

#include <vector>
#include <cmath>

/* The rawor class object is initialized with the number of data particles (e.g.
 * dark matter particles, halos, or galaxies), the number of random particles
 * the number of spherical shells in which pairs of particles are binned, the
 * volume of the box, the maximum search radius for pairs, and optionally the 
 * minimum search radius (default value of 0).
 * 
 * The public functions allow the user to change the number of data or random
 * particles, as well as obtain arrays (std::vector<int>) of the RRR, DRR, and
 * DDR counts. These arrays will contain many unused elements which correspond
 * to permutations of the side lengths of triangles. This was done for 
 * compatibility with other software that was written by DWP.
 * 
 * Private data members:
 *      nbar_ran    = The average number density of random points. This is computed
 *                    on initialization from the number of random points and the 
 *                    volume of the box.
 *      Delta_r     = The width of the spherical shells used for binning particles.
 *                    This is computed on initialization from the maximum and 
 *                    minimum search radii along with the number of shells
 *      r_max       = Maximum search radius
 *      r_min       = Minimum search radius (optional, default value 0)
 *      V_box       = Volume of the box containing the data particles
 *      N_parts     = Number of data particles
 *      N_rans      = Number of random particles
 *      N_shells    = Number of spherical shells
 *      rs          = The central radii of the spherical shells
 *      w           = Weights needed for the Gaussian quadrature integration
 *      x           = Abscissae needed for Gaussian quadrature integration
 */
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
        
        std::vector<int> getRRR();
        
        std::vector<int> getDRR();
        
        std::vector<int> getDDR(std::vector<int> &DD);
};

// (Private) Simple function that initializes rs, w, and x at the time the rawor object is 
// initialized.
void rawor::initVectors() {
    for (int i = 0; i < rawor::N_shells; ++i) {
        rawor::rs.push_back(rawor::r_min + (i + 0.5)*rawor::Delta_r);
    }
    
    rawor::w = {0.8888888888888888, 0.5555555555555556, 0.5555555555555556};
    
    rawor::x = {0.0000000000000000, -0.7745966692414834, 0.7745966692414834};
}

// (Private) Simple function that swaps two doubles if the first is larger than the second. 
// This is used to ensure that the r which is feed into rawor::sphereOverlapVolume is in 
// fact smaller than R so that volume evaluates correctly.
void rawor::swapIfGreater(double &a, double &b) {
    if (a > b) {
        double temp = a;
        a = b;
        b = temp;
    }
}

// (Private) Function for computing the volume of overlap of two spheres. The conditional 
// statements ensure that the volume is calculated correctly for cases where the two spheres
// are not touching, and when the smaller sphere is completely within the larger sphere.
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

// (Private) Function for computing the overlap volume of two spherical shells from the 
// overlap volumes of the various spheres that make up the shells.
double rawor::crossSectionVolume(double r1, double r2, double r3) {
    double V_oo = sphereOverlapVolume(r1, r3 + 0.5*rawor::Delta_r, r2 + 0.5*rawor::Delta_r);
    double V_oi = sphereOverlapVolume(r1, r3 + 0.5*rawor::Delta_r, r2 - 0.5*rawor::Delta_r);
    double V_io = sphereOverlapVolume(r1, r3 - 0.5*rawor::Delta_r, r2 + 0.5*rawor::Delta_r);
    double V_ii = sphereOverlapVolume(r1, r3 - 0.5*rawor::Delta_r, r2 - 0.5*rawor::Delta_r);
    
    return V_oo - V_oi - V_io + V_ii;
}

// (Private) Function that returns the appropriate number of permutations for a bin based 
// on the side lengths.
int rawor::getPermutations(double r1, double r2, double r3) {
    int perm = 1;
    if (r1 != r2 && r1 != r3 && r2 != r3) {
        perm = 6;
    } else if ((r1 == r2 && r1 != r3) || (r1 == r3 && r1 != r2) || (r2 == r3 && r2 != r1)) {
        perm = 3;
    }
    return perm;
}


// (Private) Function that returns the volume of a spherical shell given the shells central
// radius and the shell width stored in class data member Delta_r
double rawor::sphericalShellVolume(double r) {
    double r_o = r + 0.5*rawor::Delta_r;
    double r_i = r - 0.5*rawor::Delta_r;
    return 4.0*M_PI*(r_o*r_o*r_o - r_i*r_i*r_i)/3.0;
}

// (Private) Function that takes the data-data pair counts and used them to linearly 
// interpolate the number density to a specific location r. The r1 value should be 
// the central radius of the shell that r falls in. The pair counts should include
// double counting but no self pairs.
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

// (Private) Function that performs the Gaussian quadrature numerical integration for the RRR
// and DRR predicted counts.
double rawor::gaussQuadCrossSection(double r1, double r2, double r3) {
    double result = 0.0;
    for (int i = 0; i < rawor::w.size(); ++i) {
        double r_1 = r1 + 0.5*rawor::Delta_r*rawor::x[i];
        result += 0.5*rawor::Delta_r*rawor::w[i]*crossSectionVolume(r_1, r2, r3)*r_1*r_1;
    }
    return result;
}

// (Private) Function that performs the Gaussian quadrature numerical integration for the DDR 
// predicted counts. This is a separate function given the need to allow for a variable number
// density in this case.
double rawor::gaussQuadCrossSectionDDR(std::vector<int> &DD, double r1, double r2, double r3) {
    double result = 0.0;
    for (int i = 0; i < rawor::w.size(); ++i) {
        double r_1 = r1 + 0.5*rawor::Delta_r*rawor::x[i];
        double nbar = nbarData(DD, r_1, r1);
        result += 0.5*rawor::Delta_r*rawor::w[i]*crossSectionVolume(r_1, r2, r3)*r_1*r_1*nbar;
    }
    return result;
}

// (Public) Initializer function that sets all of the data members to the provided values. The 
// value of rMin can be omitted if one is using the default value of 0.
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

// (Public) Function that resets the value of the data member N_parts. This should be useful
// when analyzing a large number of simulations where the exact number of particles may 
// vary slightly between the realizations. Since the values of the numbers of shells, volume
// of the box and search radii would be unlikely to change, the use can simply use this
// to update the number of particles instead of declaring a new class object for each
// simulation.
void rawor::setNumParts(int numParticles) {
    rawor::N_parts = numParticles;
}

// (Public) Function that resets the value of the data member N_rans. This should be useful
// for the same reasons outlined for rawor::updateNumParts. If the number of particles 
// changes between realizations and one wants to keep the ratio of the number of particles
// to the number of randoms constant, it would be necessary to also change this value.
void rawor::setNumRans(int numRandoms) {
    rawor::N_rans = numRandoms;
    rawor::nbar_ran = numRandoms/rawor::V_box;
}

// (Public) Function that resets the value of the data member N_shells. This will force
// a recalculation of the Delta_r data member.
void rawor::setNumShells(int numShells) {
    rawor::N_shells = numShells;
    rawor::Delta_r = (rawor::r_max - rawor::r_min)/rawor::N_shells;
}

// (Public) Function that resets the value of the data member r_max. This will force
// a recalculation of the Delta_r data member.
void rawor::setRMax(double rMax) {
    rawor::r_max = rMax;
    rawor::Delta_r = (rawor::r_max - rawor::r_min)/rawor::N_shells;
}

// (Public) Function that resets the value of the data member r_min. This will force
// a recalculation of the Delta_r data member.
void rawor::setRMin(double rMin) {
    rawor::r_min = rMin;
    rawor::Delta_r = (rawor::r_max - rawor::r_min)/rawor::N_shells;
}

// (Public) Function that resets the value of the data member V_box. This will force
// a recalculation of the nbar_ran data member.
void rawor::setVBox(double VBox) {
    rawor::V_box = VBox;
    rawor::nbar_ran = rawor::N_rans/rawor::V_box;
}

// (Public) Returns the current value of N_parts
int rawor::getNumParts() {
    return rawor::N_parts;
}

// (Public) Returns the current value of N_rans
int rawor::getNumRans() {
    return rawor::N_rans;
}

// (Public) Returns the current value of N_shells
int rawor::getNumShells() {
    return rawor::N_shells;
}

// (Public) Returns the current value of V_box
double rawor::getVBox() {
    return rawor::V_box;
}

// (Public) Returns the current value of r_min
double rawor::getRMin() {
    return rawor::r_min;
}

// (Public) Returns the current value of r_max
double rawor::getRMax() {
    return rawor::r_max;
}

// (Public) Function that returns the RRR predicted counts. The array size will be N_shells^3
// even though many elements (those which are permutations or the filled in r1, r2, r3 bins)
// will be zero. There are two reasons this decision was made. First, this allows one to 
// very simply calculate the associated bin index given the values of r1, r2, and r3 (sorted
// so that r1 <= r2 <= r3 to access the non-zero bin). Second, because of the simplicity in
// calculating the index for such a padded array, this was necessary for compatibility with
// other software written by DWP for counting the DDD triangles.
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

// (Public) Function that returns the DRR predicted counts. The array size will be N_shells^3.
// See the comment/documentation of rawor::getRRR() for the reasons behind this decision.
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

// (Public) Function that returns the DDR predicted counts. The array size will be N_shells^3
// for the same reasons given previously in the comment/documentation for rawor::getRRR(). The 
// main differences between this function and the other two is that the permutations have to
// calculated explicitly since the number density will change based on which side you set as
// the distance between the data-data pair, and that you actually have to call this function
// with additional information, namely the data-data pair counts. These pair counts should
// include double counting (you can still use tricks that avoid explicit double counting for
// speed if needed, just make sure to count up by 2 instead of 1). The reason for needing
// the double counting again comes down to compatibility issues with other software written
// by DWP (their pair counting algorithm was implemented on the GPU, where memory access
// pecularities made it more efficient to double count than not).
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

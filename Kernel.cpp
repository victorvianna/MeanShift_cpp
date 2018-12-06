//
// Created by victorvianna on 06/12/18.
//

#include "Kernel.h"

double epanechnikov_kernel(double distance, double kernel_bandwidth) {
    double temp = (distance <= kernel_bandwidth) ? 1 - (distance * distance) / (kernel_bandwidth * kernel_bandwidth)
                                                 : 0;
    return temp;
}

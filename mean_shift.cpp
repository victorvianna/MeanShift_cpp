#include <stdio.h>
#include <math.h>
#include "mean_shift.h"
#include<string>

using namespace std;

#define CLUSTER_EPSILON 0.5

double MeanShift::get_distance_squared(const vector<double> &point_a, const vector<double> &point_b){
    int n = point_a.size();
    if(n != point_b.size() || (metric_weights != vector<double>{} && n != metric_weights.size())){
        throw std::string("Invalid vector sizes for metric");
    }
    vector<double> weights = (metric_weights == vector<double>{}) ? vector<double>(n, 1) : metric_weights;
    double total = 0;
    for(int i=0; i<n; i++){
        const double temp = (point_a[i] - point_b[i]);
        total += weights[i]*temp*temp;
    }
    return total;
}

double MeanShift::get_distance(const vector<double> &point_a, const vector<double> &point_b){
    return sqrt(get_distance_squared(point_a, point_b));
}

double MeanShift::gaussian_kernel(double distance, double kernel_bandwidth){
    double temp =  exp(-1.0/2.0 * (distance*distance) / (kernel_bandwidth*kernel_bandwidth));
    return temp;
}

MeanShift::MeanShift(double (*_kernel_func)(double,double),
        std::vector<double> _metric_weights){
    metric_weights = _metric_weights;
    kernel_func = _kernel_func;
}


void MeanShift::shift_point(const Point &point,
                            const std::vector<Point> &points,
                            double kernel_bandwidth,
                            Point &shifted_point) {
    shifted_point.resize( point.size() ) ;
    for(int dim = 0; dim<shifted_point.size(); dim++){
        shifted_point[dim] = 0;
    }
    double total_weight = 0;
    for(int i=0; i<points.size(); i++){
        const Point& temp_point = points[i];
        double distance = get_distance(point, temp_point);
        double weight = kernel_func(distance, kernel_bandwidth);
        for(int j=0; j<shifted_point.size(); j++){
            shifted_point[j] += temp_point[j] * weight;
        }
        total_weight += weight;
    }

    const double total_weight_inv = 1.0/total_weight;
    for(int i=0; i<shifted_point.size(); i++){
        shifted_point[i] *= total_weight_inv;
    }
}

std::vector<MeanShift::Point> MeanShift::meanshift(const std::vector<Point> &points,
                                             double kernel_bandwidth,
                                             double EPSILON){
    const double EPSILON_SQR = EPSILON*EPSILON;
    vector<bool> stop_moving(points.size(), false);
    vector<Point> shifted_points = points;
    double max_shift_distance;
    Point point_new;
    do {
        max_shift_distance = 0;
        for(int i=0; i<points.size(); i++){
            if (!stop_moving[i]) {
                shift_point(shifted_points[i], points, kernel_bandwidth, point_new);
                double shift_distance_sqr = get_distance_squared(point_new, shifted_points[i]);
                if(shift_distance_sqr > max_shift_distance){
                    max_shift_distance = shift_distance_sqr;
                }
                if(shift_distance_sqr <= EPSILON_SQR) {
                    stop_moving[i] = true;
                }
                shifted_points[i] = point_new;
            }
        }
        //printf("max_shift_distance: %f\n", sqrt(max_shift_distance));
    } while (max_shift_distance > EPSILON_SQR);
    return shifted_points;
}

vector<Cluster> MeanShift::cluster(const std::vector<Point> &points,
    const std::vector<Point> &shifted_points)
{
    vector<Cluster> clusters;

    for (int i = 0; i < shifted_points.size(); i++) {

        int c = 0;
        for (; c < clusters.size(); c++) {
            if (get_distance(shifted_points[i], clusters[c].mode) <= CLUSTER_EPSILON) {
                break;
            }
        }

        if (c == clusters.size()) {
            Cluster clus;
            clus.mode = shifted_points[i];
            clusters.push_back(clus);
        }

        clusters[c].original_points.push_back(points[i]);
        clusters[c].shifted_points.push_back(shifted_points[i]);
    }

    return clusters;
}

vector<Cluster> MeanShift::cluster(const std::vector<Point> &points, double kernel_bandwidth){
    vector<Point> shifted_points = meanshift(points, kernel_bandwidth);
    return cluster(points, shifted_points);
}

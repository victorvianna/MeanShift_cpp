#pragma once 

#include <vector>

struct Cluster {
    std::vector<double> mode;
    std::vector<std::vector<double> > original_points;
    std::vector<std::vector<double> > shifted_points;
};

class MeanShift {
public:
    typedef std::vector<double> Point;
    MeanShift(double (*_kernel_func)(double,double) = gaussian_kernel,
            std::vector<double> _metric_weights = {});
    std::vector<Point> meanshift(const std::vector<Point> & points,
                                                double kernel_bandwidth,
                                                double EPSILON = 0.00001);
    std::vector<Cluster> cluster(const std::vector<Point> &, double);

private:
    static double gaussian_kernel(double distance, double kernel_bandwidth);
    double (*kernel_func)(double,double);
    std::vector<double> metric_weights;
    void shift_point(const Point&, const std::vector<Point> &, double, Point&);
    std::vector<Cluster> cluster(const std::vector<Point> &, const std::vector<Point> &);
    double get_distance_squared(const std::vector<double> &point_a, const std::vector<double> &point_b);
    double get_distance(const std::vector<double> &point_a, const std::vector<double> &point_b);
};

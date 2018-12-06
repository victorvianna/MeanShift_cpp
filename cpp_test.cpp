#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include "MeanShift.h"
#include "Kernel.h"

using namespace std;

vector<vector<double> > load_points(const char *filename) {
    vector<vector<double> > points;
    ifstream stream_file(filename);
    while (!stream_file.eof()) {
        string line;
        getline(stream_file, line);
        istringstream stream_line(line);
        vector<double> point;
        while(!stream_line.eof()){
            double coord; char separator;
            stream_line>>coord;
            stream_line>>separator;
            point.push_back(coord);
        };
        points.push_back(point);
    }
    stream_file.close();
    return points;
}

void print_points(vector<vector<double> > points){
    for(int i=0; i<points.size(); i++){
        for(int dim = 0; dim<points[i].size(); dim++) {
            printf("%f ", points[i][dim]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    MeanShift *msp = new MeanShift(epanechnikov_kernel);
    double kernel_bandwidth = 3;

    vector<vector<double> > points = load_points("test.csv");
    vector<Cluster> clusters = msp->cluster(points, kernel_bandwidth);

    FILE *fp = fopen("result.csv", "w");
    if(!fp){
        perror("Couldn't write result.csv");
        exit(0);
    }

    printf("\n====================\n");
    printf("Found %lu clusters\n", clusters.size());
    printf("====================\n\n");
    for(int cluster = 0; cluster < clusters.size(); cluster++) {
      printf("Cluster %i:\n", cluster);
      for(int point = 0; point < clusters[cluster].original_points.size(); point++){
        for(int dim = 0; dim < clusters[cluster].original_points[point].size(); dim++) {
          printf("%f ", clusters[cluster].original_points[point][dim]);
          fprintf(fp, dim?",%f":"%f", clusters[cluster].original_points[point][dim]);
        }
        printf(" -> ");
        for(int dim = 0; dim < clusters[cluster].shifted_points[point].size(); dim++) {
          printf("%f ", clusters[cluster].shifted_points[point][dim]);
        }
        printf("\n");
        fprintf(fp, "\n");
      }
      printf("\n");
    }
    fclose(fp);

    return 0;
}

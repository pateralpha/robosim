#include <stdio.h>
#include <vector>
#include <array>
#include <cmath>
#include "la.h"
#include "circle.h"

int main(void){

    std::vector<fvector<4>> xd = {
        {0, 0, 0, 0},
        {2, 1.5, pi/4, pi/2},
        {2, 2.5, pi/4, pi/2},
        {1, 4, 0, 3*pi/4}
    };

    std::vector<xd_t> ret;
    printf("%d\r\n", make_trj(xd, ret));

    for(auto v : ret)v.print();

    FILE *fp = fopen("c.dat", "w");
    if(!fp){
        printf("file cannot open.\r\n");
        exit(EXIT_FAILURE);
    }

    printf("save to file...%s\r\n", "c.dat");
    fprintf(fp, "# control points\r\n");

    for(auto v : ret){
        if(v.type == xd_t::Line){
            fprintf(fp, "0 %6.3f %6.3f 0 0 0 ", v.xd[0], v.xd[1]);
        }
        else {
            fprintf(fp, "1 %6.3f %6.3f %6.3f %6.3f %6.3f ", v.arc[0], v.arc[1], v.arc[2], v.arc[3], v.arc[4]);
        }
        fprintf(fp, "\r\n");
    }

    fclose(fp);
    printf("file closed.\r\n");

    return 0;
}
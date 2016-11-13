/**
 *  main.cpp
 *  the main procedure to use c++ get the data and use matlab to plot the image
 *  Shock/Entropy Wave Interaction
 *
 *  Created by David Wang on 16/11/5.
 *
 */

#include "solver.hpp"

int main(int argc, const char * argv[]) {
    solver solver;
    
    // 对比有无MUSCL Limiter
    solver.setLimiter(0);
    solver.solve();
    solver.output("nonMuscl.txt");
    solver.setLimiter(1);
    solver.solve();
    solver.output("hasMuscl.txt");
    
    // 对比Kappa值
    solver.setKappa(-1.0);
    solver.solve();
    solver.output("kappa-1.txt");
    solver.setKappa(0.0);
    solver.solve();
    solver.output("kappa0.txt");
    solver.setKappa(1.0/3.0);
    solver.solve();
    solver.output("kappa13.txt");
    solver.setKappa(1.0);
    solver.solve();
    solver.output("kappa1.txt");
    
    // 对比Limiter类型
    solver.setLimiter(1);
    solver.solve();
    solver.output("vanLeer.txt");
    solver.setLimiter(2);
    solver.solve();
    solver.output("vanAlbada.txt");
    solver.setLimiter(3);
    solver.solve();
    solver.output("minmod.txt");
    solver.setLimiter(4);
    solver.solve();
    solver.output("superbee.txt");
    
    return 0;
}

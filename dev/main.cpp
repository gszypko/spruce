#include "equationset.hpp"
#include "idealmhd.hpp"
#include "plasmadomain.hpp"
#include <iostream>
#include <memory>

int main(int argc, char *argv[]){
    std::cout << "foo\n";
    PlasmaDomain pd("/Users/gszypko/Desktop/mhdtoy/dev/foobar/","/Users/gszypko/Desktop/mhdtoy/mhs_corona.config","/Users/gszypko/Desktop/mhdtoy/mhs_ar.state",false,true);
    std::cout << "bar\n";
    std::unique_ptr<EquationSet> m_eq(new IdealMHD(pd));
    std::cout << "how\n";
    std::cout << m_eq->grid("rho") << std::endl;
    std::cout << "now\n";
    return 0;
}
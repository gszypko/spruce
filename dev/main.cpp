#include "equationset.hpp"
#include "idealmhd.hpp"
#include "plasmadomain.hpp"
#include <iostream>
#include <memory>

int main(int argc, char *argv[]){
    PlasmaDomain pd("./foobar_path","./mhs_corona.config","./mhs_ar.state",false,true);
    std::unique_ptr<EquationSet> m_eq(new IdealMHD(pd));
    return 0;
}
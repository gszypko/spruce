#include "boundaryoutflow.hpp"
#include "plasmadomain.hpp"
#include "solarutils.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <sstream>
#include <iostream>
#include <cmath>

BoundaryOutflow::BoundaryOutflow(PlasmaDomain &pd) : Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    accel_template = Grid::Zero(1,1);
}

void BoundaryOutflow::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "max_accel") max_accel = std::stod(this_rhs);
        else if(this_lhs == "falloff_length") falloff_length = std::stod(this_rhs);
        else if(this_lhs == "boundary") boundary = this_rhs;
        else if(this_lhs == "falloff_shape") falloff_shape = this_rhs;
        else if(this_lhs == "feather_length") feather_length = std::stod(this_rhs);
        else if(this_lhs == "field_aligned_mode") field_aligned_mode = (this_rhs == "true");        
        else if(this_lhs == "dynamic_mode") dynamic_mode = (this_rhs == "true");
        else if(this_lhs == "dynamic_time") dynamic_time = std::stod(this_rhs);
        else if(this_lhs == "dynamic_target_speed") dynamic_target_speed = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void BoundaryOutflow::setupModule(){
    assert((falloff_shape == "exp" || falloff_shape == "gaussian" || falloff_shape == "flat") && "BoundaryOutflow shape must be exp or gaussian or flat");
    accel_template = constructBoundaryAccel(falloff_length,feather_length);
    mean_outflow = computeMeanOutflow();
}

void BoundaryOutflow::postIterateModule(double dt){
    curr_accel = max_accel;
    mean_outflow = computeMeanOutflow();
    if(dynamic_mode){
        curr_accel = (dynamic_target_speed - mean_outflow)/dynamic_time;
        curr_accel = std::min(curr_accel,max_accel);
        curr_accel = std::max(curr_accel,0.0);
    }
    Grid accel = curr_accel*accel_template;
    if (boundary=="x_bound_1" || boundary=="x_bound_2"){
        if(field_aligned_mode){
            m_pd.m_eqs->grid(IdealMHD::mom_x) += (dt*accel*m_pd.m_eqs->grid(IdealMHD::b_hat_x)*m_pd.m_eqs->grid(IdealMHD::rho));
            m_pd.m_eqs->grid(IdealMHD::mom_y) += (dt*accel*m_pd.m_eqs->grid(IdealMHD::b_hat_y)*m_pd.m_eqs->grid(IdealMHD::rho));
        }
        else{
            m_pd.m_eqs->grid(IdealMHD::mom_x) += (dt*accel*m_pd.m_eqs->grid(IdealMHD::rho));
        }
    }
    else if (boundary=="y_bound_1" || boundary=="y_bound_2"){
        if(field_aligned_mode){
            m_pd.m_eqs->grid(IdealMHD::mom_x) += (dt*accel*m_pd.m_eqs->grid(IdealMHD::b_hat_x)*m_pd.m_eqs->grid(IdealMHD::rho));
            m_pd.m_eqs->grid(IdealMHD::mom_y) += (dt*accel*m_pd.m_eqs->grid(IdealMHD::b_hat_y)*m_pd.m_eqs->grid(IdealMHD::rho));
        }
        else{
            m_pd.m_eqs->grid(IdealMHD::mom_y) += (dt*accel*m_pd.m_eqs->grid(IdealMHD::rho));
        }
    }
    m_pd.m_eqs->propagateChanges();
}

std::string BoundaryOutflow::commandLineMessage() const
{
    std::string result = boundary + " boundary outflow enforced (max " + std::to_string(mean_outflow) + " cm/s outflow)";
    result = result + " (accel. " + std::to_string(curr_accel) +" cm/s^2)";
    return result;
}

Grid BoundaryOutflow::constructBoundaryAccel(double length, double feather) const
{
    const Grid& x = m_pd.grid("pos_x");
    const Grid& y = m_pd.grid("pos_y");

    Grid ghost_mask_withboundary = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    int xl = m_pd.m_xl, xu = m_pd.m_xu, yl = m_pd.m_yl, yu = m_pd.m_yu;
    if (boundary=="x_bound_1") xl = m_pd.m_xl_dt;
    else if (boundary=="x_bound_2") xu = m_pd.m_xu_dt;
    else if (boundary=="y_bound_2") yu = m_pd.m_yu_dt;
    else if (boundary=="y_bound_1") yl = m_pd.m_yl_dt;
    else assert(false && "BoundaryOutflow boundary config must be {x,y}_bound_{1,2}");
    for(int i = xl; i <= xu; i++){
        for(int j = yl; j <= yu; j++){
        ghost_mask_withboundary(i,j) = 1.0;
        }
    }

    Grid result = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    if(falloff_shape == "exp"){
        if (boundary=="x_bound_1") result = (-2.3*(x - x.min())/length).exp();
        else if (boundary=="x_bound_2") result = (2.3*(x - x.max())/length).exp();
        else if (boundary=="y_bound_2") result = (2.3*(y - y.max())/length).exp();
        else if (boundary=="y_bound_1") result = (-2.3*(y - y.min())/length).exp();
    } else if(falloff_shape == "gaussian"){
        if (boundary=="x_bound_1") result = (-2.3*((x - x.min())/length).square()).exp();
        else if (boundary=="x_bound_2") result = (-2.3*((-x + x.max())/length).square()).exp();
        else if (boundary=="y_bound_2") result = (-2.3*((-y + y.max())/length).square()).exp();
        else if (boundary=="y_bound_1") result = (-2.3*((y - y.min())/length).square()).exp();
    } else if(falloff_shape == "flat"){
        double x_extremum = 0.0, y_extremum = 0.0;
        if (boundary=="x_bound_1") x_extremum = x.min();
        else if (boundary=="x_bound_2") x_extremum = x.max();
        else if (boundary=="y_bound_2") y_extremum = y.max();
        else if (boundary=="y_bound_1") y_extremum = y.min();
        if (boundary=="x_bound_1" || boundary=="x_bound_2") {
            for(int i = xl; i <= xu; i++){
                for(int j = yl; j <= yu; j++){
                    if(std::abs(x(i,j) - x_extremum) <= length) result(i,j) = 1.0;
                }
            }
        }
        if (boundary=="y_bound_1" || boundary=="y_bound_2") {
            for(int i = xl; i <= xu; i++){
                for(int j = yl; j <= yu; j++){
                    if(std::abs(y(i,j) - y_extremum) <= length) result(i,j) = 1.0;
                }
            }
        }
    } else {
        assert(false && "Unexpected option for BoundaryOutflow falloff_shape");
    }

    if(feather > 0.0){
        if (boundary=="x_bound_1" || boundary=="x_bound_2") result *= ((-2.3*((y - (y.max()-2.0*feather)).max(0.0)/feather).square()).exp()
                                                                        *(-2.3*((y - (y.min()+2.0*feather)).min(0.0)/feather).square()).exp()
                                                                          -0.01).max(0.0);
        else if (boundary=="y_bound_1" || boundary=="y_bound_2") result *= ((-2.3*((x - (x.max()-2.0*feather)).max(0.0)/feather).square()).exp()
                                                                            *(-2.3*((x - (x.min()+2.0*feather)).min(0.0)/feather).square()).exp()
                                                                             -0.01).max(0.0);
    }
    return ghost_mask_withboundary*result;
}

double BoundaryOutflow::computeMeanOutflow() const
{
    const Grid& x = m_pd.grid("pos_x");
    const Grid& y = m_pd.grid("pos_y");

    int xl = m_pd.m_xl, xu = m_pd.m_xu, yl = m_pd.m_yl, yu = m_pd.m_yu;
    if (boundary=="x_bound_1") xl = m_pd.m_xl_dt;
    else if (boundary=="x_bound_2") xu = m_pd.m_xu_dt;
    else if (boundary=="y_bound_2") yu = m_pd.m_yu_dt;
    else if (boundary=="y_bound_1") yl = m_pd.m_yl_dt;

    if(boundary=="x_bound_1" || boundary=="x_bound_2"){
        //apply feathering offsets
        for(int j = yl; j <= yu; j++){
            if(y(0,j) - y(0,yl) >= feather_length){
                yl = j;
                break;
            }
        }
        for(int j = yu; j >= yl; j--){
            if(y(0,yu) - y(0,j) >= feather_length){
                yu = j;
                break;
            }
        }
        //apply falloff offset
        if(boundary=="x_bound_1"){
            for(int i = xl; i <= xu; i++){
                if(x(i,0) - x(xl,0) >= falloff_length){
                    xu = i;
                    break;
                }
            }
        }
        else{ //(boundary=="x_bound_2")
            for(int i = xu; i >= xl; i--){
                if(x(xu,0) - x(i,0) >= falloff_length){
                    xl = i;
                    break;
                }
            }
        }
    }
    else{
        assert(boundary=="y_bound_1" || boundary=="y_bound_2" && "Unexpected boundary specified");
        //apply feathering offset
        for(int i = xl; i <= xu; i++){
            if(x(i,0) - x(xl,0) >= feather_length){
                xl = i;
                break;
            }
        }
        for(int i = xu; i >= xl; i--){
            if(x(xu,0) - x(i,0) >= feather_length){
                xu = i;
                break;
            }
        }
        //apply falloff offset
        if(boundary=="y_bound_1"){
            for(int j = yl; j <= yu; j++){
                if(y(0,j) - y(0,yl) >= falloff_length){
                    yu = j;
                    break;
                }
            }
        }
        else{ //(boundary=="y_bound_2")
            for(int j = yu; j >= yl; j--){
                if(y(0,yu) - y(0,j) >= falloff_length){
                    yl = j;
                    break;
                }
            }
        }
    }
    
    const Grid& b_hat_x = m_pd.m_eqs->grid(IdealMHD::b_hat_x);
    const Grid& b_hat_y = m_pd.m_eqs->grid(IdealMHD::b_hat_y);
    const Grid& v_x = m_pd.m_eqs->grid(IdealMHD::v_x);
    const Grid& v_y = m_pd.m_eqs->grid(IdealMHD::v_y);

    double max = -1.0*dynamic_target_speed;
    for(int i = xl; i <= xu; i++){
        for(int j = yl; j <= yu; j++){
            double curr = b_hat_x(i,j)*v_x(i,j) + b_hat_y(i,j)*v_y(i,j);
            //correct for orientation of field to ensure we are targeting outflow
            if (boundary=="x_bound_1" && b_hat_x(i,j) > 0.0) curr *= -1.0;
            else if (boundary=="x_bound_2" && b_hat_x(i,j) < 0.0) curr *= -1.0;
            else if (boundary=="y_bound_2" && b_hat_y(i,j) < 0.0) curr *= -1.0;
            else if (boundary=="y_bound_1" && b_hat_y(i,j) > 0.0) curr *= -1.0;
            max = std::max(max,curr);
        }
    }
    return max;
}


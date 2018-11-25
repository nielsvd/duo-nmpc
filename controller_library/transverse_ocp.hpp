/*
 *    This file is part of duo-nmpc.
 *
 *    Duo-nmpc
 *    Copyright (C) 2018 Niels van Duijkeren, KU Leuven.
 *
 *    Duo-nmpc is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    Duo-nmpc is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with duo-nmpc; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef TRANSVERSE_OCP_H_
#define TRANSVERSE_OCP_H_

#include "ocp.hpp"

extern "C" {
#include "casadi_interface.h"
}

class TransverseOcp : public Ocp {
public:
    TransverseOcp();
    virtual ~TransverseOcp();

    void Shift() override;

    double ComputeLyapunovBound(double epsN2);
    std::vector<double> ComputeTunnelBound(int Nta, double epsN2);
protected:
    casadi_memory *tunnel_bound_memory;
};

#endif // TRANSVERSE_OCP_H_
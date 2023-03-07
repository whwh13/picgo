#ifndef MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvFF_hpp
#define MEDYAN_Mechanics_ForceField_VolumeConservation_VolConsrvFF_hpp

#include <memory>
#include <vector>

#include "Mechanics/ForceField/VolumeConservation/VolConsrvMembrane.hpp"

namespace medyan {

struct VolumeConservationFFFactory {

    auto operator()(
        const std::string& type
    ) const {
        using namespace std;

        std::vector< std::unique_ptr< ForceField > > res;

        if (type == "MEMBRANE") {
            res.push_back(
                make_unique< VolumeConservationMembrane >());
        }
        else if(type == "") {}
        else {
            log::error("Volume conservation FF type {} is not recognized.", type);
            throw std::runtime_error("Membrane volume conservation FF type not recognized");
        }

        return res;
    }
};

} // namespace medyan

#endif

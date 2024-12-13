/*
 * Volcanic Skies
 * Copyright (C) 2024 P. Cilliers Pretorius, University of Cape Town
 *
 * This file is part of the Volcanic Skies project.
 *
 * Volcanic Skies is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License (GPL) as published by the Free Software 
 * Foundation, either version 2 of the License, or (at your discretion) any later version.
 *
 * Volcanic Skies is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with 
 * Volcanic Skies in the LICENSE file. If not, see <https://www.gnu.org/licenses/>.
 *
 * Additional information and disclaimers regarding liability and third-party 
 * components can be found in the NOTICE file included with this project.
 *
 */
#ifndef PWM_UTILS_CLOUDTYPE_H
#define PWM_UTILS_CLOUDTYPE_H
namespace PWM{
    namespace Utils{
        enum class cloudType{
            EMPTY,
            STRATUS,
            ALTOSTRATUS,
            CIRROSTRATUS,
            NIMBOSTRATUS,
            CIRRUS,
            CUMULUS,
            STRATOCUMULUS,
            ALTOCUMULUS,
            CIRROCUMULUS,
            CUMULONIMBUS,
            CTYPEEND
        };
    }
}
#endif // PWM_UTILS_CLOUDTYPE_H

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

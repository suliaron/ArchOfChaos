#include "model.h"

const char *Model::FormalismToString(Formalism formalism) noexcept
{
    switch (formalism) {
        case Formalism::NEWTONIAN:
            return "NEWTONIAN";

        case Formalism::HAMILTONIAN:
            return "HAMILTONIAN";
    }

    return "UNKNOWN";
}

const char *Model::IndicatorTypeToString(IndicatorType indicator) noexcept
{
    switch (indicator) {
        case IndicatorType::NONE:
            return "NONE";

        case IndicatorType::FLI:
            return "FLI";

        case IndicatorType::LCI:
            return "LCI";

        case IndicatorType::RLI:
            return "RLI";
    }

    return "UNKNOWN";
}
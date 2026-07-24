#include <cstdint>
#include <iostream>

struct InitData {
    // Initial epoch [day]
    double t0 = 0.0;

    // Gaussian gravitational constant
    double mu = 0.0;

    // Semimajor axis grid
    double a0 = 0.0;
    double a1 = 0.0;
    std::uint32_t Na = 0;

    // Eccentricity grid
    double e0 = 0.0;
    double e1 = 0.0;
    std::uint32_t Ne = 0;

    // Fixed orbital elements
    double omega = 0.0;
    double tau = 0.0;

    InitData() = default;
    ~InitData() = default;
};

class GridIterator {
public:
    explicit GridIterator(const InitData& data)
        : data_(data),
        ia_(0),
        ie_(0),
        da_((data.a1 - data.a0) / static_cast<double>(data.Na)),
        de_((data.e1 - data.e0) / static_cast<double>(data.Ne))
    {}

    double a() const
    {
        return data_.a0 + static_cast<double>(ia_) * da_;
    }

    double e() const
    {
        return data_.e0 + static_cast<double>(ie_) * de_;
    }

    bool next()
    {
        if (ia_ < data_.Na) {
            ++ia_;
            return true;
        }

        ia_ = 0;

        if (ie_ < data_.Ne) {
            ++ie_;
            return true;
        }

        return false;
    }

private:
    const InitData& data_;

    std::uint32_t ia_;
    std::uint32_t ie_;

    double da_;
    double de_;
};

int main()
{
    InitData init;

    init.a0 = 1.0;
    init.a1 = 2.0;
    init.Na = 4;

    init.e0 = 0.0;
    init.e1 = 0.2;
    init.Ne = 2;

    GridIterator grid(init);

    do {
        std::cout << "a = " << grid.a()
            << ", e = " << grid.e() << '\n';
    } while (grid.next());

    return 0;
}
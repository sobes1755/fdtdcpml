#include <fstream>
#include <experimental/mdspan>

template <typename T>
requires std::is_trivially_copyable_v <T>
void
writeBinary(std::ofstream& ofs, T val)
{

    ofs.write(reinterpret_cast<const char*>(&val), sizeof(T));

}

template <typename T, typename Extents, typename Layout, typename Accessor>
void
writeSlice(std::ofstream& ofs, std::mdspan<T, Extents, Layout, Accessor> slice)
{

    writeBinary(ofs, slice.rank());

    if constexpr (slice.rank() == 1) {

        writeBinary(ofs, slice.extent(0));

        for (std::size_t p = 0; p < slice.extent(0); ++p)
            writeBinary(ofs, slice[p]);

    } else if constexpr (slice.rank() == 2) {

        writeBinary(ofs, slice.extent(0));
        writeBinary(ofs, slice.extent(1));

        for (std::size_t p = 0; p < slice.extent(0); ++p)
        for (std::size_t q = 0; q < slice.extent(1); ++q)
            writeBinary(ofs, slice[p, q]);

    } else if constexpr (slice.rank() == 3) {

        writeBinary(ofs, slice.extent(0));
        writeBinary(ofs, slice.extent(1));
        writeBinary(ofs, slice.extent(2));

        for (std::size_t p = 0; p < slice.extent(0); ++p)
        for (std::size_t q = 0; q < slice.extent(1); ++q)
        for (std::size_t r = 0; r < slice.extent(2); ++r)
            writeBinary(ofs, slice[p, q, r]);

    }

}

template <typename T>
struct scaling_accessor {

    using offset_policy = scaling_accessor;
    using element_type = T;
    using reference = T;
    using data_handle_type = T*;

    T coefficient;

    constexpr reference access(data_handle_type p, std::size_t i) const noexcept {
        return p[i] * coefficient;
    }

    constexpr data_handle_type offset(data_handle_type p, std::size_t i) const noexcept {
        return p + i;
    }

};

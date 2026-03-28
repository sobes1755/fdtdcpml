#pragma once

#include <fstream>
#include <experimental/mdspan>

namespace fdtdcpml {

template <typename T>
struct ScalingAccessor {

    using element_type = T;
    using data_handle_type = T*;
    using reference = T;
    using offset_policy = ScalingAccessor<T>;

    T scale_ = T{1};

    ScalingAccessor() = default;  // They say it's used by mdspan...
    explicit ScalingAccessor(T scale) : scale_(scale) {}  // They say "explicit" is recommended by C++ Core Guidelines...

    template <typename OtherT>
    ScalingAccessor(const ScalingAccessor<OtherT>& other) noexcept : scale_(other.scale_) {}  // Used by submdspan...

    constexpr reference access(data_handle_type p, std::size_t i) const noexcept {
        return p[i] * scale_;
    }
    constexpr data_handle_type offset(data_handle_type p, std::size_t i) const noexcept {
        return p + i;
    }

};

template <typename T>
requires std::is_trivially_copyable_v <T>
void
write_binary(std::ofstream& ofs, T val)
{

    ofs.write(reinterpret_cast<const char*>(&val), sizeof(T));

}

template <typename T, typename Extents, typename Layout, typename Accessor>
void
write_mdspan(std::ofstream& ofs, std::mdspan<T, Extents, Layout, Accessor> span)
{

    write_binary(ofs, span.rank());

    for (std::size_t i = 0; i < span.rank(); ++i) {
        write_binary(ofs, span.extent(i));
    }

    auto write_recursive = [&](auto& self, auto... indices) mutable {
        if constexpr (sizeof...(indices) == Extents::rank()) {
            write_binary(ofs, span[indices...]);
        } else {
            std::size_t dim = sizeof...(indices);
            for (std::size_t i = 0; i < span.extent(dim); ++i) {
                self(self, indices..., i);
            }
        }
    };

    write_recursive(write_recursive);

}

}  // namespace fdtdcpml

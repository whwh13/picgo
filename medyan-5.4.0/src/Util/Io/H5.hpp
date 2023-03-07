#ifndef MEDYAN_Util_Io_H5_hpp
#define MEDYAN_Util_Io_H5_hpp

#include <algorithm>
#include <string>
#include <string_view>
#include <type_traits>

#include <highfive/H5Attribute.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <xtensor/xarray.hpp>

namespace medyan::h5 {

using File = HighFive::File;
using Group = HighFive::Group;
using DataSpace = HighFive::DataSpace;
using Attribute = HighFive::Attribute;
using DataSet = HighFive::DataSet;

// Some type traits.
//------------------------------------------------------------------------------

// Dense Eigen types.
template< typename T >
constexpr bool isEigenDense = std::is_base_of_v<Eigen::DenseBase<T>, T>;

// General scalar types accepted by the library.
template< typename T >
constexpr bool isGeneralScalar =
    std::is_arithmetic_v<T> ||
    std::is_same_v<T, std::string>;

// Add multi-layer pointer.
template< typename T, Size n >
struct AddPointerN {
    using type = std::add_pointer_t<typename AddPointerN<T, n-1>::type>;
};
template< typename T >
struct AddPointerN<T, 0> {
    using type = T;
};
template< typename T, Size n >
using AddPointerN_t = typename AddPointerN<T, n>::type;

// Auxiliary functions for reading and writing data sets.
//------------------------------------------------------------------------------

// General scalar types.
template< typename T, std::enable_if_t< isGeneralScalar<T> >* = nullptr >
inline void writeDataSet(Group& group, std::string_view name, const T& data) {
    using namespace HighFive;
    const std::string namestr(name);
    DataSet dataset = group.exist(namestr) ? group.getDataSet(namestr) : group.createDataSet<T>(namestr, DataSpace::From(data));
    dataset.write(data);
}
template< typename T, std::enable_if_t< isGeneralScalar<T> >* = nullptr >
inline void readDataSet(T& data, const Group& group, std::string_view name) {
    using namespace HighFive;
    DataSet dataset = group.getDataSet(string(name));
    dataset.read(data);
}

// std::vector types.
template< typename T >
inline void writeDataSet(Group& group, std::string_view name, const std::vector<T>& data) {
    using namespace HighFive;
    const std::string namestr(name);
    // Empty cases are handled differently because of the way HighFive handles empty vectors. https://github.com/BlueBrain/HighFive/issues/172
    DataSpace dataspace = data.empty() ? DataSpace(0) : DataSpace::From(data);
    DataSet dataset = group.exist(namestr) ? group.getDataSet(namestr) : group.createDataSet<T>(namestr, dataspace);
    dataset.write(data);
}

template< typename T >
inline void readDataSet(std::vector<T>& data, const Group& group, std::string_view name) {
    using namespace HighFive;
    DataSet dataset = group.getDataSet(std::string(name));
    DataSpace dataspace = dataset.getSpace();
    data.resize(dataspace.getDimensions()[0]);
    dataset.read(data);
}

// xtensor data types.
template< std::size_t n, typename T >
inline void writeDataSet(Group& group, std::string_view name, const xt::xtensor<T, n, xt::layout_type::column_major>& data) {
    using namespace HighFive;
    const std::string namestr(name);
    const auto shape = data.shape();
    DataSpace dataspace = DataSpace(shape.rbegin(), shape.rend());
    DataSet dataset = group.exist(namestr) ? group.getDataSet(namestr) : group.createDataSet<T>(namestr, dataspace);
    dataset.write((AddPointerN_t<const T, n>) data.data());
}
template< std::size_t n, typename T >
inline void readDataSet(xt::xtensor<T, n, xt::layout_type::column_major>& data, const Group& group, std::string_view name) {
    using namespace HighFive;
    DataSet dataset = group.getDataSet(std::string(name));
    auto dims = dataset.getSpace().getDimensions();
    std::reverse(dims.begin(), dims.end());
    data.resize(dims);
    dataset.read((AddPointerN_t<T, n>) data.data());
}

// Auxiliary function to store Eigen matrices in a HDF5 dataset.
// Currently, HighFive gives the same dimension as in the HDF5 file, by NOT respecting the default column major order of the Eigen matrix.
// Using this function, if an Eigen matrix is (m, n) column major order, then the HDF5 dataset will be (n, m) row major. As a result, the memory layout will be the same.
template< typename EigenMatrix, std::enable_if_t< isEigenDense<EigenMatrix> >* = nullptr >
inline void writeDataSet(Group& grp, std::string_view name, const EigenMatrix& mat) {
    using namespace std;
    using Scalar = typename EigenMatrix::Scalar;
    const std::string namestr(name);
    DataSet dataset = grp.exist(namestr)
        ? grp.getDataSet(namestr)
        : grp.createDataSet<Scalar>(namestr, HighFive::DataSpace{ static_cast<size_t>(mat.cols()), static_cast<size_t>(mat.rows()) });
    dataset.write((const Scalar**)mat.data());
}
template< typename EigenMatrix, std::enable_if_t< isEigenDense<EigenMatrix> >* = nullptr >
inline void readDataSet(EigenMatrix& mat, const Group& grp, std::string_view name) {
    using namespace std;
    using Scalar = typename EigenMatrix::Scalar;
    DataSet dataset = grp.getDataSet(string(name));
    const auto sizes = dataset.getSpace().getDimensions();
    mat.resize(sizes[1], sizes[0]);
    dataset.read((Scalar**)mat.data());
}



// Alternative read dataset function that creates and returns the data.
template< typename T >
inline T readDataSet(const Group& group, std::string_view name) {
    T t;
    readDataSet(t, group, name);
    return t;
}

} // namespace medyan::h5

#endif

#ifndef MEDYAN_Util_StableVector_hpp
#define MEDYAN_Util_StableVector_hpp

#include <cstddef> // nullptr_t, size_t
#include <functional> // hash
#include <optional>
#include <type_traits>
#include <utility> // forward, move
#include <vector>


namespace medyan {


// Notes:
// - The StableVector should allow erasing while iterating.
template< typename T >
struct StableVector {
    struct Index {
        using ValueType = std::ptrdiff_t;
        ValueType value = 0;

        Index() = default;
        Index(ValueType value) : value(value) {}

        Index& operator++() { ++value; return *this; }
        Index& operator--() { --value; return *this; }

        // Replace with <=> in C++20.
        bool operator==(Index rhs) const noexcept { return value == rhs.value; }
        bool operator!=(Index rhs) const noexcept { return !(*this == rhs); }
    };

    // Notes:
    // - Must have "Index index" member, because the index can be used in iterator related functions.
    template< bool isConst >
    struct IteratorBase {
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type        = T;
        using difference_type   = std::ptrdiff_t;
        using pointer           = std::conditional_t< isConst, const value_type*, value_type* >;
        using reference         = std::conditional_t< isConst, const value_type&, value_type& >;

        using ParentPointer = std::conditional_t< isConst, const StableVector*, StableVector* >;
        ParentPointer psv = nullptr;
        Index index {};

        IteratorBase() = default;
        IteratorBase(ParentPointer psv, Index newIndex) :
            psv(psv), index(newIndex)
        {
            // Find the next dereferenceable target, or the end if no dereferenceable target is found.
            while(index.value < psv->value.size() && !psv->value[index.value].has_value()) {
                ++index;
            }
        }
        IteratorBase(const IteratorBase&) = default;
        IteratorBase& operator=(const IteratorBase&) = default;

        // Dereferencing
        //---------------------------------------------
        reference operator*() const { return (*psv)[index]; }
        pointer operator->() const { return &(*psv)[index]; }

        // Modification
        //---------------------------------------------
        // Precondition: *this is dereferenceable
        IteratorBase& operator++() {
            do {
                ++index;
            }
            while(index.value < psv->value.size() && !psv->value[index.value].has_value());
            return *this;
        }
        // Precondition: *this is decrementable
        IteratorBase& operator--() {
            do {
                --index;
            }
            while(index.value >= 0 && !psv->value[index.value].has_value());
            return *this;
        }
        IteratorBase operator++(int) { IteratorBase tmp(*this); ++(*this); return tmp; }
        IteratorBase operator--(int) { IteratorBase tmp(*this); --(*this); return tmp; }

        // Comparison
        //---------------------------------------------
        template< bool isConstRhs >
        bool operator==(const IteratorBase<isConstRhs>& rhs) const {
            return psv == rhs.psv && index.value == rhs.index.value;
        }
        template< bool isConstRhs >
        bool operator!=(const IteratorBase<isConstRhs>& rhs) const {
            return !(*this == rhs);
        }
    }; // IteratorBase

    using Iterator      = IteratorBase< false >;
    using ConstIterator = IteratorBase< true >;

    // All data.
    //---------------------------------
    std::vector< std::optional<T> > value;
    std::vector< Index > deletedIndices;


    // Member functions.
    //---------------------------------

    // Insert an element.
    //
    // Note:
    // - Invalidates iterators.
    Index insert(const T& t) {
        Index res;
        if(deletedIndices.empty()) {
            // Append at the back.
            res.value = value.size();
            value.push_back(t);
        } else {
            // Fill in the last hole.
            res = deletedIndices.back();
            value[res.value] = t;
            deletedIndices.pop_back();
        }
        return res;
    }
    Index insert(T&& t) {
        Index res;
        if(deletedIndices.empty()) {
            // Append at the back.
            res.value = value.size();
            value.push_back(std::move(t));
        } else {
            // Fill in the last hole.
            res = deletedIndices.back();
            value[res.value] = std::move(t);
            deletedIndices.pop_back();
        }
        return res;
    }

    // Insert an element in place.
    template< typename... Args >
    Index emplace(Args&&... args) {
        Index res;
        if(deletedIndices.empty()) {
            // Append at the back.
            res.value = value.size();
            value.emplace_back(std::in_place, std::forward<Args>(args)...);
        } else {
            // Fill in the last hole.
            res = deletedIndices.back();
            value[res.value].emplace(std::forward<Args>(args)...);
            deletedIndices.pop_back();
        }
        return res;
    }

    // Remove an element at index.
    //
    // Note:
    // - If index is within value range but the item is already deleted, this function is no-op.
    // - If index is outside value range, the behavior is undefined.
    // - Does not invalidate iterators.
    void erase(Index index) {
        if(value[index.value].has_value()) {
            value[index.value].reset();
            deletedIndices.push_back(index);
        }
    }

    // Remove an element at iterator.
    //
    // Note:
    // - The iterator must point to a valid element, or the behavior is undefined.
    // - Does not invalidate iterators except for this one.
    void erase(Iterator it) {
        erase(indexat(it));
    }

    // Access an element at index.
    // Note:
    // - If index is outside valid range of value, or item at index is deleted, the behavior is undefined.
    T&       operator[](Index index)       { return *value[index.value]; }
    const T& operator[](Index index) const { return *value[index.value]; }

    // Access an element at index safely.
    // Note:
    // - If index is outside valid range of value, or item at index is deleted, an exception will be thrown.
    T&       at(Index index)       { return value.at(index.value).value(); }
    const T& at(Index index) const { return value.at(index.value).value(); }

    // Get number of elements stored.
    auto size() const { return value.size() - deletedIndices.size(); }

    // Returns true if no elements are stored.
    bool empty() const { return size() == 0; }

    // Iterator start.
    Iterator      begin()        noexcept { return Iterator     (this, Index {0}); }
    ConstIterator begin()  const noexcept { return ConstIterator(this, Index {0}); }
    ConstIterator cbegin() const noexcept { return ConstIterator(this, Index {0}); }

    // Iterator end.
    Iterator      end()        noexcept { return Iterator     (this, Index {static_cast<std::ptrdiff_t>(value.size())}); }
    ConstIterator end()  const noexcept { return ConstIterator(this, Index {static_cast<std::ptrdiff_t>(value.size())}); }
    ConstIterator cend() const noexcept { return ConstIterator(this, Index {static_cast<std::ptrdiff_t>(value.size())}); }

    // Utility function, that converts a valid iterator to an index.
    Index indexat(Iterator      it) const { return it.index; }
    Index indexat(ConstIterator it) const { return it.index; }
};

template< typename T >
using StableVectorIndex = typename StableVector<T>::Index;

// Provides hasing for indices.
template< typename T >
struct StableVectorIndexHash {
    std::size_t operator()(StableVectorIndex<T> ti) const {
        return std::hash<typename StableVectorIndex<T>::ValueType>{}(ti.value);
    }
};

} // namespace medyan

#endif

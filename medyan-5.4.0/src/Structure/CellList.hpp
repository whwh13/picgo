#ifndef MEDYAN_Structure_CellList_Hpp
#define MEDYAN_Structure_CellList_Hpp

#include <cstddef> // ptrdiff_t
#include <iterator> // tags
#include <vector>

#include "common.h"

namespace medyan {

// Forward declarations
template< typename TElement, typename THead > class CellListManager;

// Storing information for user of elements,
// like each molecule.
// Note:
//   - The users must make sure that the manager points to a valid location
template< typename TElement, typename THead >
struct CellListElementUser {
    using ManagerType = CellListManager< TElement, THead >;

    ManagerType* manager = nullptr;
    Index head;
    Index index;
};

// Storing information for user of heads,
// like each compartment storing molecules.
// Note:
//   - The users must make sure that the manager points to a valid location
template< typename TElement, typename THead >
struct CellListHeadUser {
    using ManagerType = CellListManager< TElement, THead >;

    ManagerType* manager = nullptr;
    Index index;
};


// The necessary structure for element list
template< typename T >
struct CellListElement {
    T val;
    Index next;
    Index prev;
    bool hasNext;
    bool hasPrev;
};

// The necessary structure for head list
template< typename T >
struct CellListHead {
    T val;
    Size  size = 0;
    Index first;
    Index last;
};

// The class managing a cell list
template< typename TElement, typename THead >
class CellListManager {
public:
    using ElementList = std::vector< CellListElement<TElement> >;
    using HeadList    = std::vector< CellListHead   <THead   > >;
    using ElementUser = CellListElementUser< TElement, THead >;
    using HeadUser    = CellListHeadUser   < TElement, THead >;

    // The class for viewing and iterating the element pointers in a specific cell.
    // The class acts like a double linked list.
    // Only const iterator is offered.
    class CellView {
    private:

        // The const iterator for traversing cell list
        class ConstIterator_ {
        public:
            using iterator_category = std::bidirectional_iterator_tag;
            using value_type        = TElement;
            using difference_type   = std::ptrdiff_t;
            using pointer           = const value_type*;
            using reference         = const value_type&;

        private:
            const ElementList* el_;
            const CellListHead<THead>* head_;
            Index ei_;
            bool atEnd_;

        public:
            ConstIterator_() = default;
            ConstIterator_(const ElementList* el, const CellListHead<THead>* head, Index elementIndex, bool atEnd) :
                el_(el), head_(head), ei_(elementIndex), atEnd_(atEnd || (head_->size == 0))
            { }
            ConstIterator_(const ConstIterator_&) = default;
            ConstIterator_& operator=(const ConstIterator_&) = default;

            // Dereferencing
            //---------------------------------------------
            reference operator*() const { return (*el_)[ei_].val; }
            pointer operator->() const { return &(*el_)[ei_].val; }

            // Modification
            //---------------------------------------------
            // Precondition: *this is dereferenceable
            ConstIterator_& operator++() {
                if((*el_)[ei_].hasNext) ei_ = (*el_)[ei_].next;
                else atEnd_ = true;
                return *this;
            }
            // Precondition: *this is decrementable
            ConstIterator_& operator--() {
                if(atEnd_) {
                    ei_ = head_->last; // head_->size == 0 is UB
                    atEnd_ = false;
                }
                else ei_ = (*el_)[ei_].prev; // (*el_)[ei_].hasPrev == false is UB
                return *this;
            }
            ConstIterator_ operator++(int) { ConstIterator_ tmp(*this); ++(*this); return tmp; }
            ConstIterator_ operator--(int) { ConstIterator_ tmp(*this); --(*this); return tmp; }

            // Comparison
            //---------------------------------------------
            bool operator==(const ConstIterator_& rhs) const {
                return (atEnd_ && rhs.atEnd_) || (!atEnd_ && !rhs.atEnd_ && ei_ == rhs.ei_);
            }
            bool operator!=(const ConstIterator_& rhs) const {
                return !(*this == rhs);
            }
        }; // class ConstIterator_

    public:
        using const_iterator = ConstIterator_;
        using size_type = Size;

        CellView(const ElementList* el, const CellListHead<THead>* head) :
            el_(el), head_(head)
        { }

        size_type size() const noexcept { return head_->size; }
        bool empty() const noexcept { return !size(); }

        const_iterator begin() const noexcept { return const_iterator(el_, head_, head_->first, false); }
        const_iterator end() const noexcept { return const_iterator(el_, head_, 0, true); }

    private:
        const ElementList* el_;
        const CellListHead<THead>* head_;
    };


private:
    HeadList                   headList_;
    ElementList                elementList_;
    std::vector< Index >       elementDeletedIndices_;

public:
    // Accessors
    //-------------------------------------------------------------------------

    THead getHead(Index headIndex) const { return headList_[headIndex].val; }
    THead getHead(const ElementUser& eu) const { return getHead(eu.head); }
    auto  numHeads() const noexcept { return headList_.size(); }

    CellView getElements(Index headIndex) const { return CellView(&elementList_, &headList_[headIndex]); }
    CellView getElements(const HeadUser& hu) const { return getElements(hu.index); }
    auto     numElements() const noexcept { return elementList_.size(); }

    // Element operations with fixed head users
    //-------------------------------------------------------------------------

    void clearElements() {
        elementList_.clear();
        elementDeletedIndices_.clear();
        for(auto& head : headList_) {
            head.size  = 0;
        }
    }

    // Returns the index of the element.
    Index addElement(TElement element, Index headIndex) {
        // Add the new element
        Index newIndex;
        if(elementDeletedIndices_.empty()) {
            elementList_.push_back({element});
            newIndex = elementList_.size() - 1;
        } else {
            newIndex = elementDeletedIndices_.back();
            elementList_[newIndex].val = element;
            elementDeletedIndices_.pop_back();
        }

        // Connect the new element
        registerElement_(newIndex, headIndex);

        return newIndex;
    }
    void addElement(TElement element, ElementUser& eu, const HeadUser& hu) {
        eu.head = hu.index;
        eu.index = addElement(element, eu.head);
    }

    void updateElement(Index elementIndex, Index oldHeadIndex, Index newHeadIndex) {
        // Remove from current head.
        unregisterElement_(elementIndex, oldHeadIndex);
        // Add to the new head.
        registerElement_(elementIndex, newHeadIndex);
    }
    void updateElement(ElementUser& eu, const HeadUser& newHu) {
        updateElement(eu.index, eu.head, newHu.index);
        eu.head = newHu.index;
    }

    void removeElement(Index elementIndex, Index headIndex) {
        unregisterElement_(elementIndex, headIndex);
        elementDeletedIndices_.push_back(elementIndex);
    }
    void removeElement(const ElementUser& eu) {
        removeElement(eu.index, eu.head);
    }

    // Head operations with fixed elements.
    //-------------------------------------------------------------------------

    Index addHead(THead head) {
        headList_.push_back({head});
        return headList_.size() - 1;
    }
    void addHead(THead head, HeadUser& hu) {
        hu.index = addHead(head);
    }

    // Clear all heads and elements.
    //-------------------------------------------------------------------------
    void clearAll() {
        clearElements();
        headList_.clear();
    }

private:

    // Helper function to register an element at an index to a head
    // Note:
    //   - This function causes the head size count to increase by 1.
    //   - This function does not manage allocating the element.
    void registerElement_(const Index index, const Index head) {
        if(headList_[head].size) {
            // Insert as last
            const auto last = headList_[head].last;
            elementList_[index].hasPrev = true;
            elementList_[index].prev = last;
            elementList_[index].hasNext = false;

            elementList_[last].hasNext = true;
            elementList_[last].next = index;
            headList_[head].last = index;
        }
        else {
            // Was empty
            elementList_[index].hasPrev = false;
            elementList_[index].hasNext = false;
            headList_[head].first = index;
            headList_[head].last  = index;
        }

        ++headList_[head].size;
    } // void registerElement_(...)

    // Helper function to unregister an element at an index from a head
    // Note:
    //   - This function causes the head size count to decrease by 1.
    //   - This function does not manage removing the element.
    void unregisterElement_(const Index index, const Index head) {
        --headList_[head].size;

        // Reconnect linked list
        const auto& e = elementList_[index];

        if(e.hasPrev) {
            elementList_[e.prev].hasNext = e.hasNext;
            if(e.hasNext) elementList_[e.prev].next = e.next;
        } else {
            if(e.hasNext) headList_[head].first = e.next;
        }

        if(e.hasNext) {
            elementList_[e.next].hasPrev = e.hasPrev;
            if(e.hasPrev) elementList_[e.next].prev = e.prev;
        } else {
            if(e.hasPrev) headList_[head].last = e.prev;
        }
    } // void unregisterElement_(...)

}; // template<...> class CellListManager



} // namespace medyan

#endif

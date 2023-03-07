
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_DoubleLinkedList_h
#define MEDYAN_DoubleLinkedList_h
#include <cstddef> // ptrdiff_t, size_t
#include <iterator> // tags
#include "Util/Io/Log.hpp"

namespace medyan {
namespace linkedlist{
	// Forward declarations
	template< typename TElement> class DoubleLinkedList;

	template< typename T >
	struct DLLNode {
		T* elementptr;
		DLLNode* prev = nullptr;
		DLLNode* next = nullptr;
		bool hasPrev = false;
		bool hasNext = false;
	};

	template< typename TElement>
	class DoubleLinkedList {

	public:

		//Checks if a Node is at the beginning of the linked list or at the end.
		bool isHeadNode(DLLNode<TElement>* TNodeptr){ return TNodeptr == headNodePtr;}
		bool isTailNode(DLLNode<TElement>* TNodeptr){ return TNodeptr == tailNodePtr;}

		//Gives the number of Nodes in the list.
		std::size_t size(){return numNodes;}

		void addNode(TElement* newElement, DLLNode<TElement>& dllnoderef){

			dllnoderef.elementptr = newElement;

			if(numNodes){
				dllnoderef.hasPrev = true;
				dllnoderef.prev = tailNodePtr;
				tailNodePtr->hasNext = true;
				tailNodePtr->next = &dllnoderef;
				tailNodePtr = &dllnoderef;
			}
			else {
				headNodePtr = &dllnoderef;
				tailNodePtr = &dllnoderef;
			}
			numNodes++;
		}

		void removeNode(DLLNode<TElement>& dllnoderef){

			if(numNodes > 2 ){
				//Connect the prevnode to nextnode
				if (isTailNode(&dllnoderef)) {
					tailNodePtr = dllnoderef.prev;
					//sever connection with previous Node
					dllnoderef.prev->hasNext = false;
					dllnoderef.prev->next = nullptr;
				}
				else if (isHeadNode(&dllnoderef)) {
					headNodePtr = dllnoderef.next;
					//sever connection with next Node
					dllnoderef.next->hasPrev = false;
					dllnoderef.next->prev = nullptr;
				}
				else {
					DLLNode<TElement> *_prevnode = dllnoderef.prev;
					DLLNode<TElement> *_nextnode = dllnoderef.next;
					_prevnode->next = _nextnode;
					_nextnode->prev = _prevnode;
				}
			}
			else if(numNodes == 2 ) {
				if (isTailNode(&dllnoderef)) {
					tailNodePtr = dllnoderef.prev;
					//sever connection with previous Node
					dllnoderef.prev->hasNext = false;
					dllnoderef.prev->next = nullptr;
				}
				else if (isHeadNode(&dllnoderef)) {
					headNodePtr = dllnoderef.next;
					//sever connection with next Node
					dllnoderef.next->hasPrev = false;
					dllnoderef.next->prev = nullptr;
				}
			}

			if(numNodes > 0) {
				//reset node data
				dllnoderef.elementptr = nullptr;
				dllnoderef.hasPrev = false;
				dllnoderef.prev = nullptr;
				tailNodePtr->hasNext = false;
				tailNodePtr->next = nullptr;
				numNodes--;
			}
			else
				LOG(ERROR) << "Cannot remove element from LinkedList of size 0.";
		}

	private:

		class ConstIterator_{
		public:
			using iterator_category = std::bidirectional_iterator_tag;
			using value_type        = TElement*;
			using difference_type   = std::ptrdiff_t;
			using pointer           = const value_type*;
			using reference         = const value_type&;

		private:
			DLLNode<TElement>* _headnodeptr;
			DLLNode<TElement>* _currentptr;
			DLLNode<TElement>* _tailnodeptr;
			bool atEnd_;

		public:
			ConstIterator_() = default;
			ConstIterator_(size_t numNodes, bool atEnd, DLLNode<TElement>* headNodePtr,
			               DLLNode<TElement>* tailNodePtr) :
					atEnd_(atEnd || numNodes == 0), _headnodeptr(headNodePtr), _tailnodeptr(tailNodePtr)
			{
				if(atEnd_)
					_currentptr = _tailnodeptr;
				else
					_currentptr = _headnodeptr;
			}
			ConstIterator_(const ConstIterator_&) = default;
			ConstIterator_& operator=(const ConstIterator_&) = default;


			// Dereferencing
			//---------------------------------------------
			reference operator*() const { return _currentptr->elementptr; }
			pointer operator->() const { return &_currentptr->elementptr; }

			// Modification
			//---------------------------------------------
			// Precondition: *this is dereferenceable
			ConstIterator_& operator++() {
				if(_currentptr->hasNext) _currentptr = _currentptr->next;
				else atEnd_ = true;
				return *this;
			}
			// Precondition: *this is decrementable
			ConstIterator_& operator--() {
				if(atEnd_) {
					_currentptr = _tailnodeptr; // head_->size == 0 is UB
					atEnd_ = false;
				}
				else if(_currentptr->hasPrev)
					_currentptr = _currentptr->prev;
				else
					_currentptr = _headnodeptr;

				// is UB
				return *this;
			}
			ConstIterator_ operator++(int) { ConstIterator_ tmp(*this); ++(*this); return tmp; }
			ConstIterator_ operator--(int) { ConstIterator_ tmp(*this); --(*this); return tmp; }

			// Comparison
			//---------------------------------------------
			bool operator==(const ConstIterator_& rhs) const {
				return (atEnd_ && rhs.atEnd_) || (!atEnd_ && !rhs.atEnd_);
			}
			bool operator!=(const ConstIterator_& rhs) const {
				return !(*this == rhs);
			}
		};

		std::size_t numNodes = 0;
		DLLNode<TElement>* headNodePtr = nullptr;
		DLLNode<TElement>* tailNodePtr = nullptr;

	public:
		using const_iterator = ConstIterator_;
		const_iterator begin() const noexcept { return const_iterator(numNodes, false,
				headNodePtr, tailNodePtr); }
		const_iterator end() const noexcept { return const_iterator(numNodes, true,
		                                                            headNodePtr, tailNodePtr); }
	};
};

} // namespace medyan

#endif //MEDYAN_DoubleLinkedList_h

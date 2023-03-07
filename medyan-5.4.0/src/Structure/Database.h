
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

#ifndef MEDYAN_Database_h
#define MEDYAN_Database_h

#include <algorithm> // max, swap
#include <cstddef> // nullptr_t, size_t
#include <tuple> // apply
#include <utility> // forward
#include <vector>


namespace medyan {

//-----------------------------------------------------------------------------
// DatabaseData stores the data that has the same indexing with the internal
// Database structure.
// DatabaseData must be a vector-like type that satisfies the following criteria.
//   - has type value_type.
//   - value_type is default constructable.
//   - value_type is copy/move assignable.
//   - has function push_back(value_type or ref).
//   - has function operator[](size_t pos) -> value_type&.
//   - has function back() -> value_type&.
//   - has function pop_back().
//   - has function resize(size_t size).
// Note: with C++20 these constraints should be written as a concept.
//-----------------------------------------------------------------------------
struct DatabaseDataDefault {
    struct value_type {};

    value_type dummy;

    void push_back(value_type) {}
    value_type&       operator[](std::size_t)       { return dummy; }
    const value_type& operator[](std::size_t) const { return dummy; }
    value_type&       back()       { return dummy; }
    const value_type& back() const { return dummy; }
    void pop_back() {}
    void resize(std::size_t) {}
};

template< typename T >
class DatabaseBase {
    
private:
    inline static std::vector<T*> elems_;  // Pointer to the elements in the collection
    inline static std::size_t nextId_ = 0; // Next unique id

    std::size_t id_;
    std::size_t index_;

public:
    static const auto& getElements() { return elems_; }
    static auto numElements() { return elems_.size(); }

    // Add element on construction
    DatabaseBase() : id_(nextId_++), index_(elems_.size()) {
        elems_.push_back(static_cast<T*>(this));
    }

    // Copy construction creates a new identity
    DatabaseBase(const DatabaseBase&) : DatabaseBase() {}
    // Move construction will use the identity of moved-from object, and a new identity is created for the moved-from object.
    DatabaseBase(DatabaseBase&& rhs ) : DatabaseBase() {
        std::swap(id_, rhs.id_);
        std::swap(index_, rhs.index_);
        std::swap(elems_[index_], elems_[rhs.index_]);
    }

    // Remove element on destruction
    ~DatabaseBase() {
        if(index_ + 1 != elems_.size()) {
            // Move the data from the last element to the current position
            elems_[index_] = elems_.back();
            // Updata _index of the original last element
            elems_[index_] -> DatabaseBase< T >::index_ = index_;
        }

        // Pop the last element
        elems_.pop_back();
    }

    // Disallow assignments
    DatabaseBase& operator=(const DatabaseBase&) = delete;
    DatabaseBase& operator=(DatabaseBase&&     ) = delete;

    std::size_t getId() const { return id_; }
    std::size_t getIndex() const { return index_; }

    // This function overrides the current id of the element.
    // One should not use it unless in cases like re-initializing the system.
    void overrideId(std::size_t id) {
        id_ = id;
        nextId_ = std::max(id_ + 1, nextId_);
    }

};

template< typename DatabaseData > class DatabaseDataManager {

private:
    inline static DatabaseData dbData_;

public:
    static auto&       getDbData()      { return dbData_; }
    static const auto& getDbDataConst() { return dbData_; }
};


/// A collection class to hold instances of a given class

/*!
 *  The Database class is a template for holding a collection of objects
 *  It is used by all elements in the SubSystem to track instances of
 *  the given class. The databases are then used in various ways, i.e.
 *  mechanical minimization uses all beads for its Minimizer methods,
 *  ForceField uses the collections to calculate forces and energy, etc.
 *  
 *  @param T - class to point to
 *  @param stableIndexing - when this is set to true, the database will be able
 *    to track an additional indexing which will never be invalidated until
 *    rearrange() is called.
 *  @param DatabaseData - the managed vectorized data. When stableIndexing is
 *    true, one uses getStableIndex() to access the data, but the size of the
 *    data might be bigger than the number of elements. When stableIndexing is
 *    false, one uses getIndex() to access the data, where the index might be
 *    invalidated after any other element is created/destroyed.
 */
template< typename T, bool stableIndexing, typename DatabaseData = DatabaseDataDefault > class Database;

// Specialization for unstable indexing
template< typename T, typename DatabaseData >
class Database< T, false, DatabaseData > : public DatabaseBase<T>,
                                           public DatabaseDataManager<DatabaseData> {

public:

    using DbBaseType = DatabaseBase<T>;
    using DbDataType = DatabaseDataManager<DatabaseData>;

    // Default constructor creates default value
    Database() {
        DbDataType::getDbData().push_back(typename DatabaseData::value_type {});
    }
    // Add element on construction
    Database(typename DatabaseData::value_type val) {
        DbDataType::getDbData().push_back(std::move(val));
    }

    // Copy constructor also copies the data.
    Database(const Database& rhs) : DbBaseType(rhs) {
        DbDataType::getDbData().push_back(DbDataType::getDbData()[rhs.getIndex()]);
    }
    // Move constructor creates default data for the moved-from object.
    Database(Database&& rhs) : DbBaseType(std::move(rhs)) {
        DbDataType::getDbData().push_back(typename DatabaseData::value_type {});
    }

    // Remove element on destruction (by swapping)
    ~Database() {
        if(this->getIndex() + 1 != DbBaseType::getElements().size()) {
            // Move the data from the last element to the current position
            DbDataType::getDbData()[this->getIndex()] = std::move(DbDataType::getDbData().back());
        }

        // Pop the last element
        DbDataType::getDbData().pop_back();
    }
};

// Specialization for stable indexing
template< typename T, typename DatabaseData >
class Database< T, true, DatabaseData > : public DatabaseBase<T>,
                                          public DatabaseDataManager<DatabaseData> {
    inline static std::vector<T*> stableElems_;
    inline static std::vector<std::size_t> deletedIndices_;

    std::size_t stableIndex_ = 0;

    // Auxiliary function to initialize new data
    template< typename CaseBack, typename CaseHole >
    void initialize_(CaseBack&& caseBack, CaseHole&& caseHole) {
        if(deletedIndices_.empty()) {
            // Append at the back
            stableIndex_ = stableElems_.size();
            stableElems_.push_back(static_cast<T*>(this));
            caseBack();
        } else {
            // Fill in the last hole
            stableIndex_ = deletedIndices_.back();
            stableElems_[stableIndex_] = static_cast<T*>(this);
            caseHole();
            deletedIndices_.pop_back();
        }
    }

public:
    using DbBaseType = DatabaseBase<T>;
    using DbDataType = DatabaseDataManager<DatabaseData>;

    static const auto& getStableElement(std::size_t pos) { return stableElems_[pos]; }
    // Get raw number of stable elements (including deleted)
    static auto rawNumStableElements() { return stableElems_.size(); }
    // Getting information for debug purposes
    static const auto& getDeletedIndices() { return deletedIndices_; }

    // Calling this function may change the stable indices.
    static void rearrange() {
        using std::size_t;

        const size_t numDeleted = deletedIndices_.size();
        const size_t currentSize = stableElems_.size();
        const size_t finalSize = currentSize - numDeleted;

        std::vector<char> isDeleted(numDeleted);
        // Mark to-be-deleted items with indices bigger than finalSize as deleted
        for(size_t i = 0; i < numDeleted; ++i)
            if(deletedIndices_[i] >= finalSize)
                isDeleted[deletedIndices_[i] - finalSize] = true;

        // Move the not-to-be-deleted items with bigger indices to the to-be-deleted items with small indices
        for(size_t indAfterFinal = 0, i = 0; indAfterFinal < numDeleted; ++indAfterFinal) {
            if(!isDeleted[indAfterFinal]) {
                while(i < numDeleted && deletedIndices_[i] >= finalSize) ++i; // Find (including current i) the next i with small index
                if(i < numDeleted) {
                    // Found. This should always be satisfied.
                    stableElems_[deletedIndices_[i]] = stableElems_[finalSize + indAfterFinal];
                    DbDataType::getDbData()[deletedIndices_[i]] = std::move(DbDataType::getDbData()[finalSize + indAfterFinal]);
                    stableElems_[deletedIndices_[i]]->stableIndex_ = deletedIndices_[i];
                }
                ++i;
            }
        }

        // Remove garbage
        stableElems_.resize(finalSize);
        DbDataType::getDbData().resize(finalSize);
        deletedIndices_.clear();
    }

    // Default constructor creates default value.
    Database() {
        initialize_(
            // pushing back empty
            [] { DbDataType::getDbData().push_back(typename DatabaseData::value_type {}); },
            // filling hole with default value
            [this] { DbDataType::getDbData()[stableIndex_] = typename DatabaseData::value_type {}; }
        );
    }
    // Constructor adds the value.
    Database(typename DatabaseData::value_type val) {
        initialize_(
            // pushing back
            [&] { DbDataType::getDbData().push_back(std::move(val)); },
            // filling the hole
            [&, this] { DbDataType::getDbData()[stableIndex_] = std::move(val); }
        );
    }

    // Copy constructor also copies the data.
    Database(const Database& rhs) : DbBaseType(rhs) {
        initialize_(
            // pushing back
            [&] { DbDataType::getDbData().push_back(DbDataType::getDbData()[rhs.stableIndex_]); },
            // filling the hole
            [&, this] { DbDataType::getDbData()[stableIndex_] = DbDataType::getDbData()[rhs.stableIndex_]; }
        );
    }
    // Move constructor creates default data for the moved-from object.
    // The moved-to (this) object will use all indices in the moved-from (rhs) object. Moved-from object will have default data.
    Database(Database&& rhs) : DbBaseType(std::move(rhs)) {
        // Create empty data slot, but will be used by rhs later because of index swap.
        initialize_(
            // pushing back
            [] { DbDataType::getDbData().push_back(typename DatabaseData::value_type {}); },
            // filling hole
            [this] { DbDataType::getDbData()[stableIndex_] = typename DatabaseData::value_type {}; }
        );
        std::swap(stableIndex_, rhs.stableIndex_);
        std::swap(stableElems_[stableIndex_], stableElems_[rhs.stableIndex_]);
    }

    ~Database() {
        // Only mark as deleted
        deletedIndices_.push_back(stableIndex_);
    }

    std::size_t getStableIndex() const { return stableIndex_; }

};

} // namespace medyan

#endif

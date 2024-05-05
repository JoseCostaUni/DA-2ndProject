
#ifndef DA_TP_CLASSES_UFDS
#define DA_TP_CLASSES_UFDS

#include <vector>
/**
 * @brief Union-Find Disjoint Set (UFDS) data structure.
 */
class UFDS {
public:
    /**
    * @brief Constructor to create UFDS object with specified number of elements.
    * @param N Number of elements in the UFDS.
    */
    UFDS(unsigned int N);
    /**
    * @brief Find the representative element of the set that contains element i.
    * @param i Element whose representative is to be found.
    * @return The representative element of the set containing i.
    */
    unsigned long findSet(unsigned int i);
    /**
   * @brief Check if two elements belong to the same set.
   * @param i First element.
   * @param j Second element.
   * @return True if i and j belong to the same set, otherwise false.
   */
    bool isSameSet(unsigned int i, unsigned int j);
    /**
    * @brief Link two sets together.
    * @param i Element from the first set.
    * @param j Element from the second set.
    */
    void linkSets(unsigned int i, unsigned int j);
private:
    std::vector<unsigned int> path; // Ancestor of node i (which can be itself). It is used to determine if two nodes are part of the same set.
    std::vector<unsigned int> rank; // Upper bound for the height of a tree whose root is node i.
};


#endif //DA_TP_CLASSES_UFDS
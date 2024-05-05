
#include "UFDS.h"

/**
   * @brief Constructor to create UFDS object with specified number of elements.
   * @param N Number of elements in the UFDS.
   */
UFDS::UFDS(unsigned int N) {
    path.resize(N);
    rank.resize(N);
    for (unsigned long i = 0; i < N; i++) {
        path[i] = i;
        rank[i] = 0;
    }
}
/**
 * @brief Find the representative element of the set that contains element i.
 * @param i Element whose representative is to be found.
 * @return The representative element of the set containing i.
 *
 * Time Complexity: O(log N) on average, where N is the number of elements in the UFDS.
 */
unsigned long UFDS::findSet(unsigned int i) {
    if (path[i] != i) path[i] = findSet(path[i]);
    return path[i];
}
/**
     * @brief Check if two elements belong to the same set.
     * @param i First element.
     * @param j Second element.
     * @return True if i and j belong to the same set, otherwise false.
     *
     * Time Complexity: O(log N) on average, where N is the number of elements in the UFDS.
     */
bool UFDS::isSameSet(unsigned int i, unsigned int j) {
    return findSet(i) == findSet(j);
}
/**
   * @brief Link two sets together.
   * @param i Element from the first set.
   * @param j Element from the second set.
   *
   * Time Complexity: O(log N) on average, where N is the number of elements in the UFDS.
   */
void UFDS::linkSets(unsigned int i, unsigned int j) {
    if (!isSameSet(i, j)) {
        unsigned long x = findSet(i), y = findSet(j);
        if (rank[x] > rank[y]) path[y] = x; // x becomes the root due to having a larger rank
        else {
            path[x] = y; // y becomes the root due to having a larger rank, or ...
            if (rank[x] == rank[y]) rank[y]++; // ... due to both nodes having the same rank (in order to break the tie)
        }
    }
}
#ifndef __PROGTEST__
#include <cassert>
#include <iomanip>
#include <cstdint>
#include <iostream>
#include <memory>
#include <limits>
#include <optional>
#include <algorithm>
#include <bitset>
#include <list>
#include <array>
#include <vector>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <queue>
#include <random>
#include <type_traits>
#include <utility>

#endif

template <typename T>
struct Queue {
  // define the node structure with necessary attributes
  struct Node {
    T data;
    // pointers to child nodes and parent node
    Node *left, *right, *parent;

    // size of the subtree rooted at this node
    size_t size;

    // priority used for treap balancing
    int priority;
    Node(const T &data)
        : data(data), left(nullptr), right(nullptr), parent(nullptr), size(1),
          priority(rand()) {}
  };

  // root pointer for the tree
  Node *root;

  // initialize the queue
  Queue() : root(nullptr) {}

  // destructor to clean up memory
  ~Queue() { clear(root); }

  // helper function to delete all nodes
  void clear(Node *node) {
    if (node) {
      clear(node->left);
      clear(node->right);
      delete node;
    }
  }

  // check if the queue is empty
  bool empty() const { return root == nullptr; }

  // return the number of elements in the queue
  size_t size() const { return get_size(root); }

  // helper function to get size of a node
  static size_t get_size(Node *node) { return node ? node->size : 0; }

  // update the size of a node after modification
  static void update_size(Node *node) {
    if (node) {
      node->size = 1 + get_size(node->left) + get_size(node->right);
    }
  }

  // reference type to represent an element in the queue
  struct Ref {
    Node *node;
    Ref(Node *node = nullptr) : node(node) {}
  };

  // add an element to the end of the queue and return its reference
  Ref push_last(T x) {
    Node *new_node = new Node(x);

    root = merge(root, new_node);

    if (root)
      root->parent = nullptr;

    return Ref(new_node);
  }

  // remove and return the first element from the queue
  T pop_first() {
    if (empty())
      throw std::out_of_range("Queue is empty");

    Node *left = nullptr, *right = nullptr;
    split(root, 1, left, right);

    T ret = left->data;

    delete left;

    root = right;

    if (root)
      root->parent = nullptr;

    return ret;
  }

  // get the position of an element by its reference
  size_t position(const Ref &ref) const { return get_position(ref.node); }

  // move an element forward in the queue by the given number of positions
  void jump_ahead(const Ref &ref, size_t positions) {
    size_t pos = position(ref);
    size_t new_pos = positions > pos ? 0 : pos - positions;

    // remove the node from its current position
    Node *left = nullptr, *middle = nullptr, *right = nullptr;
    split(root, pos, left, right);
    split(right, 1, middle, right); // middle is the node being moved
    root = merge(left, right);      // merge the remaining parts
    if (root)
      root->parent = nullptr;

    // insert the node at the new position
    split(root, new_pos, left, right);
    root = merge(merge(left, middle), right);
    if (root)
      root->parent = nullptr;
  }

private:
  // merge two subtrees into one based on priority
  static Node *merge(Node *left, Node *right) {
    if (!left)
      return right; // if left subtree is empty, return right subtree

    if (!right)
      return left; // if right subtree is empty, return left subtree

    if (left->priority > right->priority) { // if left has higher priority
      left->right = merge(left->right, right); // merge right subtree into left's right

      if (left->right)
        left->right->parent = left; // update parent pointer

      update_size(left);
      return left;
    } else {
      // merge left subtree into right's left
      right->left = merge(left, right->left);

      if (right->left)
        right->left->parent = right; // update parent pointer

      update_size(right); // update size of the current node
      return right;
    }
  }

  // split a tree into two parts based on a key
  static void split(Node *node, size_t key, Node *&left, Node *&right,
                    size_t add = 0) {
    if (!node) {
      left = right = nullptr;
      return;
    }

    size_t curr_key = add + get_size(node->left); // calculate the current key of the node

    // if the key is less than or equal to the current key
    if (key <= curr_key) {
      split(node->left, key, left, node->left, add); // split the left subtree

      if (node->left)
        node->left->parent = node;

      right = node;
    } else {
      split(node->right, key, node->right, right, curr_key + 1);

      if (node->right)
        node->right->parent = node;

      left = node;
    }
    update_size(node);
  }

  // calculate the position of a node in the tree
  size_t get_position(const Node *node) const {
    size_t pos = get_size(node->left);
    const Node *current = node;
    while (current->parent) {
      if (current == current->parent->right) {
        pos += get_size(current->parent->left) + 1;
      }
      current = current->parent;
    }
    return pos;
  }
};

#ifndef __PROGTEST__

////////////////// Dark magic, ignore ////////////////////////

template < typename T >
auto quote(const T& t) { return t; }

std::string quote(const std::string& s) {
  std::string ret = "\"";
  for (char c : s) if (c != '\n') ret += c; else ret += "\\n";
  return ret + "\"";
}

#define STR_(a) #a
#define STR(a) STR_(a)

#define CHECK_(a, b, a_str, b_str) do { \
    auto _a = (a); \
    decltype(a) _b = (b); \
    if (_a != _b) { \
      std::cout << "Line " << __LINE__ << ": Assertion " \
        << a_str << " == " << b_str << " failed!" \
        << " (lhs: " << quote(_a) << ")" << std::endl; \
      fail++; \
    } else ok++; \
  } while (0)

#define CHECK(a, b) CHECK_(a, b, #a, #b)

#define CHECK_EX(expr, ex) do { \
    try { \
      (expr); \
      fail++; \
      std::cout << "Line " << __LINE__ << ": Expected " STR(expr) \
        " to throw " #ex " but no exception was raised." << std::endl; \
    } catch (const ex&) { ok++; \
    } catch (...) { \
      fail++; \
      std::cout << "Line " << __LINE__ << ": Expected " STR(expr) \
        " to throw " #ex " but got different exception." << std::endl; \
    } \
  } while (0)

////////////////// End of dark magic ////////////////////////


void test1(int& ok, int& fail) {
  Queue<int> Q;
  CHECK(Q.empty(), true);
  CHECK(Q.size(), 0);

  constexpr int RUN = 10, TOT = 105;

  for (int i = 0; i < TOT; i++) {
    Q.push_last(i % RUN);
    CHECK(Q.empty(), false);
    CHECK(Q.size(), i + 1);
  }

  for (int i = 0; i < TOT; i++) {
    CHECK(Q.pop_first(), i % RUN);
    CHECK(Q.size(), TOT - 1 - i);
  }

  CHECK(Q.empty(), true);
}

void test2(int& ok, int& fail) {
  Queue<int> Q;
  CHECK(Q.empty(), true);
  CHECK(Q.size(), 0);
  std::vector<decltype(Q.push_last(0))> refs;

  constexpr int RUN = 10, TOT = 105;

  for (int i = 0; i < TOT; i++) {
    refs.push_back(Q.push_last(i % RUN));
    CHECK(Q.size(), i + 1);
  }

  for (int i = 0; i < TOT; i++) CHECK(Q.position(refs[i]), i);

  Q.jump_ahead(refs[0], 15);
  Q.jump_ahead(refs[3], 0);

  CHECK(Q.size(), TOT);
  for (int i = 0; i < TOT; i++) CHECK(Q.position(refs[i]), i);

  Q.jump_ahead(refs[8], 100);
  Q.jump_ahead(refs[9], 100);
  Q.jump_ahead(refs[7], 1);

  static_assert(RUN == 10 && TOT >= 30);
  for (int i : { 9, 8, 0, 1, 2, 3, 4, 5, 7, 6 })
    CHECK(Q.pop_first(), i);

  for (int i = 0; i < TOT*2 / 3; i++) {
    CHECK(Q.pop_first(), i % RUN);
    CHECK(Q.size(), TOT - 11 - i);
  }

  CHECK(Q.empty(), false);
}

template < int TOT >
void test_speed(int& ok, int& fail) {
  Queue<int> Q;
  CHECK(Q.empty(), true);
  CHECK(Q.size(), 0);
  std::vector<decltype(Q.push_last(0))> refs;

  for (int i = 0; i < TOT; i++) {
    refs.push_back(Q.push_last(i));
    CHECK(Q.size(), i + 1);
  }

  for (int i = 0; i < TOT; i++) {
    CHECK(Q.position(refs[i]), i);
    Q.jump_ahead(refs[i], i);
  }

  for (int i = 0; i < TOT; i++) CHECK(Q.position(refs[i]), TOT - 1 - i);

  for (int i = 0; i < TOT; i++) {
    CHECK(Q.pop_first(), TOT - 1 - i);
    CHECK(Q.size(), TOT - 1 - i);
  }

  CHECK(Q.empty(), true);
}

int main() {
  int ok = 0, fail = 0;
  if (!fail) test1(ok, fail);
  if (!fail) test2(ok, fail);
  if (!fail) test_speed<1'000>(ok, fail);

  if (!fail) std::cout << "Passed all " << ok << " tests!" << std::endl;
  else std::cout << "Failed " << fail << " of " << (ok + fail) << " tests." << std::endl;
}

#endif

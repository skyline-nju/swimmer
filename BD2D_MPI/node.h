#ifndef NODE_H
#define NODE_H
#include <rand.h>
#include <list>

template <class TPar>
class UniNode : public TPar {
public:
  UniNode() : TPar(), next(nullptr) {}
  UniNode(double x0, double y0, double vx0, double vy0): TPar(x0, y0, vx0, vy0), next(nullptr) {}
  template <typename TBc> 
  UniNode(Ran &myran, const TBc &bc): TPar(myran, bc), next(nullptr) {}

  void append_at_front(UniNode<TPar> ** head) {
    next = *head;
    *head = this;
  }

  void break_away(UniNode<TPar> *pre_node) const { pre_node->next = next; }
  void break_away(UniNode<TPar> **head, UniNode<TPar> *pre_node) const;
  
  UniNode *next;
};

template<class TPar>
void UniNode<TPar>::break_away(UniNode<TPar>** head,
                               UniNode<TPar>* pre_node) const {
  if (pre_node) {
    pre_node->next = next;
  } else {
    *head = next;
  }
}

template <class TPar>
class BiNode : public TPar {
public:
  BiNode() : TPar(), prev(nullptr), next(nullptr) {}
  BiNode(double x0, double y0, double vx0, double vy0): TPar
  (x0, y0, vx0, vy0), prev(nullptr), next(nullptr) {}
  template <typename TBc>
  BiNode(Ran &myran, const TBc &bc): TPar(myran, bc), prev(nullptr), next(nullptr) {}
  
  void append_at_front(BiNode<TPar> ** head);

  void break_away() const;

  void break_away(BiNode<TPar> **head) const;

  void break_away(BiNode<TPar> **head, BiNode<TPar> *pre_node) const {
    break_away(head);
  }

  BiNode *prev;
  BiNode *next;
};

template <class TPar>
void BiNode<TPar>::append_at_front(BiNode<TPar> ** head) {
  prev = nullptr;
  next = *head;
  if (next) {
    next->prev = this;
  }
  *head = this;
}

template<class TPar>
void BiNode<TPar>::break_away() const {
  prev->next = next;
  if (next) {
    next->prev = prev;
  }
}

template<class TPar>
void BiNode<TPar>::break_away(BiNode<TPar>** head) const {
  if (prev) {
    break_away();
  } else {
    *head = next;
    if (next) {
      next->prev = nullptr;
    }
  }
}

template <class TNode, class BiFunc>
void for_each_node_pair(TNode* head, BiFunc f_ij) {
  TNode *node1 = head;
  while (node1->next) {
    TNode *node2 = node1->next;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

template <class TPar, class BiFunc>
void for_each_node_pair(const std::list<TPar *> &cl, BiFunc f_ij) {
  auto end2 = cl.cend();
  auto end1 = std::prev(end2);
  for (auto it1 = cl.cbegin(); it1 != end1; ++it1) {
    for (auto it2 = std::next(it1); it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template <class TNode, class BiFunc>
void for_each_node_pair(TNode* head1, TNode* head2, BiFunc f_ij) {
  TNode *node1 = head1;
  do {
    TNode *node2 = head2;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  } while (node1);
}

template <class TPar, class BiFunc>
void for_each_node_pair(const std::list<TPar *> &cl1,
                        const std::list<TPar *> &cl2, BiFunc f_ij) {
  auto end1 = cl1.cend();
  auto end2 = cl2.cend();
  auto beg2 = cl2.cbegin();
  for (auto it1 = cl1.cbegin(); it1 != end1; ++it1) {
    for (auto it2 = beg2; it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}
#endif
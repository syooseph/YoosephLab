/* The basic implementation of tree structure */
#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <queue>
#include <stack>

#ifndef _TREE_BASIC_H_
#define _TREE_BASIC_H_

template<typename DataType>
struct TreeNodeType {
  DataType data;
  TreeNodeType<DataType>* parent;
  std::list<TreeNodeType<DataType>*> children;
};

template<typename DataType>
class Tree  {
 public:
  Tree(void);
  ~Tree(void);
  TreeNodeType<DataType>* GetRoot(void);
  TreeNodeType<DataType>* AddChild(const DataType child_data, TreeNodeType<DataType>* parent);
  void TraverseDepthFirst(std::list<DataType>& recorded_data);
  void TraverseBreathFirst(std::list<DataType>& recorded_data);
  void FetchLeaves(std::list<TreeNodeType<DataType>*>& leaves);
  void SpellPath(TreeNodeType<DataType>* path_end, std::list<DataType>& path_data);
  void Destroy(void);
 private:
  TreeNodeType<DataType>* root_;
};

template<typename DataType>
Tree<DataType>::Tree(void)  {
  root_ = NULL;
  return;
}

template<typename DataType>
Tree<DataType>::~Tree(void) {
  return;
}

template<typename DataType>
TreeNodeType<DataType>* Tree<DataType>::AddChild(const DataType child_data, TreeNodeType<DataType>* parent)  {
  TreeNodeType<DataType>* child = new TreeNodeType<DataType>;
  child->data = child_data;
  child->parent = parent;
  if(parent != NULL)  {
    parent->children.push_back(child);
  } else  {
    root_ = child;
  }
  return child;
}

template<typename DataType>
TreeNodeType<DataType>* Tree<DataType>::GetRoot(void) {
  return root_;
}

template<typename DataType>
void Tree<DataType>::TraverseBreathFirst(std::list<DataType>& recorded_data)  {
  if(root_ == NULL)  {
    return;
  }
  std::queue<TreeNodeType<DataType>*> unvisited;
  unvisited.push(root_);
  while(!unvisited.empty()) {
    TreeNodeType<DataType>* current = unvisited.front();
    recorded_data.push_back(current->data);
    for(auto it = current->children.begin(); it != current->children.end(); ++ it) {
      unvisited.push(*it);
    }
    unvisited.pop();
  }
  return;
}

template<typename DataType>
void Tree<DataType>::TraverseDepthFirst(std::list<DataType>& recorded_data) {
  if(root_ == NULL)  {
    return;
  }
  std::stack<TreeNodeType<DataType>*> unvisited;
  unvisited.push(root_);
  while(!unvisited.empty()) {
    TreeNodeType<DataType>* current = unvisited.top();
    recorded_data.push_back(current->data);
    unvisited.pop();
    for(auto it = current->children.rbegin(); it != current->children.rend(); ++ it) {
      unvisited.push(*it);
    }
  }
  return;
}

template<typename DataType>
void Tree<DataType>::Destroy(void)  {
  if(root_ == NULL)  {
    return;
  }
  std::stack<TreeNodeType<DataType>*> unvisited;
  unvisited.push(root_);
  while(!unvisited.empty()) {
    TreeNodeType<DataType>* current = unvisited.top();
    unvisited.pop();
    for(auto it = current->children.rbegin(); it != current->children.rend(); ++ it) {
      unvisited.push(*it);
    }
    delete current;
  }
  return;
}

template<typename DataType>
void Tree<DataType>::FetchLeaves(std::list<TreeNodeType<DataType>*>& leaves)  {
  if(root_ == NULL)  {
    return;
  }
  std::stack<TreeNodeType<DataType>*> unvisited;
  unvisited.push(root_);
  while(!unvisited.empty()) {
    TreeNodeType<DataType>* current = unvisited.top();
    // consider the node as a leaf if it has no child
    if(current->children.size() == 0)  {
      leaves.push_back(current);
    }
    unvisited.pop();
    for(auto it = current->children.rbegin(); it != current->children.rend(); ++ it) {
      unvisited.push(*it);
    }
  }
  return;
}

template<typename DataType>
void Tree<DataType>::SpellPath(TreeNodeType<DataType>* path_end, std::list<DataType>& path_data)  {
  TreeNodeType<DataType>* current = path_end;
  while(current != NULL) {
    path_data.push_back(current->data);
    current = current->parent;
  }
  return;
}

#endif

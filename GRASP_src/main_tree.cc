#include "tree_basic.h"
#include <iostream>
#include <list>
#include <vector>

using namespace std;

int main()  {
  Tree<int> tree_progress;
  list<TreeNodeType<int>* > to_expand;
  list<TreeNodeType<int>* > new_expand;
  // insert the root node
  to_expand.push_back(tree_progress.AddChild(22, NULL));
  cout << "before the loop" << endl;
  unsigned int i;
  for(i = 0; i < 3; ++ i) {
    cout << "stage  " << i << endl;
    cout << "size of expand list: " << to_expand.size() << endl;
    for(auto it = to_expand.begin(); it != to_expand.end(); ++ it) {
      int current_data = (*it)->data;
      cout << "current_data:  " << current_data << endl;
      new_expand.push_back(tree_progress.AddChild(current_data * 2, *it));
      new_expand.push_back(tree_progress.AddChild(current_data * 2 + 1, *it));
      cout << "end of for loop" << endl;
    }
    to_expand = new_expand;
    new_expand.clear();
  }
  cout << tree_progress.GetRoot()->children.size() << endl;
  cout << endl;
  list<int> record_nums;
  tree_progress.TraverseDepthFirst(record_nums);
  list<TreeNodeType<int>*> leaves;
  tree_progress.FetchLeaves(leaves);
  for(auto it = leaves.begin(); it != leaves.end(); ++ it) {
    list<int> path_data;
    tree_progress.SpellPath(*it, path_data);
    for(auto it_list = path_data.begin(); it_list != path_data.end(); ++ it_list) {
      cout << *it_list << " ";
    }
    cout << endl;
  }
  cout << endl;
  return 0;
}

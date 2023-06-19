#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include "utils.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

bool in_pathway(std::vector<int> pathway, int i) {
  auto temp = std::find(pathway.begin(), pathway.end(), i);
  if (temp == pathway.end()) {
    return false;
  } else {
    return true;
  }
}

void remove_duplicates(List &tree, int &changed, std::vector<int> pathway = {}, int current_pos = 0) {
  
  IntegerVector nodes = tree["nodes"];
  List children = tree["children"];
  LogicalVector in_tree = tree["in_tree"];
  
  std::vector<int> current_ch_set = children.at(current_pos);
  int ch_set_size = current_ch_set.size();
  
  // check if the mutation has already occurred in the pathway
  if (in_pathway(pathway, nodes.at(current_pos))) {
    in_tree.at(current_pos) = false;
    std::vector<int> parents = tree["parents"];
    int pa_idx = parents.at(current_pos) - 1;
    std::vector<int> pa_ch_set = children.at(pa_idx);
    auto to_remove = std::find(pa_ch_set.begin(), pa_ch_set.end(), current_pos + 1);
    pa_ch_set.erase(to_remove);
    
    if (ch_set_size > 0) {
      pa_ch_set.insert(pa_ch_set.end(), current_ch_set.begin(), current_ch_set.end());
      for (int c : current_ch_set) {
        parents.at(c - 1) = pa_idx + 1;
      }
    }
    children.at(pa_idx) = pa_ch_set;
    
    tree["parents"] = parents;
    tree["children"] = children;
    changed++;
  }
  
  if (ch_set_size > 0) {
    // check if duplicated mutations are in the same subclone
    for (int k {0}; k < (ch_set_size - 1); ++k) {
      int ch1 = current_ch_set.at(k) - 1; // c++ indexing starts from 0
      if (in_tree.at(ch1)) {
        for (int l {k + 1}; l < ch_set_size; ++l) {
          int ch2 = current_ch_set.at(l) - 1; // c++ indexing starts from 0
          if (nodes.at(ch1) == nodes.at(ch2)) {
            in_tree.at(ch2) = false;
            std::vector<int> ch1_set = children.at(ch1);
            std::vector<int> ch2_set = children.at(ch2);
            ch1_set.insert(ch1_set.end(), ch2_set.begin(), ch2_set.end());
            children.at(ch1) = ch1_set;
            std::vector<int> parents = tree["parents"];
            for (int c : ch2_set) {
              parents.at(c - 1) = ch1 + 1;
            }
            tree["in_tree"] = in_tree;
            tree["children"] = children;
            tree["parents"] = parents;
            changed++;
          }
        }
      }
    }
    if (changed > 0) {
      std::vector<int> new_current_ch_set {};
      for (int k {0}; k < ch_set_size; ++k) {
        int ch1 = current_ch_set.at(k) - 1; // c++ indexing starts from 0
        if (in_tree.at(ch1)) {
          new_current_ch_set.push_back(ch1 + 1);
        }
      }
      children.at(current_pos) = new_current_ch_set;
      current_ch_set = new_current_ch_set;
    }
    
    if (current_pos != 0) {
      pathway.push_back(nodes.at(current_pos));
    }
    
    // go down the tree
    for (int pos : current_ch_set) {
      remove_duplicates(tree, changed, pathway, pos - 1); // c++ indexing starts from 0
    }
  }
}

void get_augmented_tree(std::vector<int> mut_indices, List &tree,
                        std::vector<int> pathway = {}, int current_pos = 0) {
  
  IntegerVector nodes = tree["nodes"];
  List children = tree["children"];
  LogicalVector in_tree = tree["in_tree"];
  IntegerVector parents = tree["parents"];
  
  if (in_tree.at(current_pos)) {
    std::vector<int> current_ch_set = children.at(current_pos);
    int ch_set_size = current_ch_set.size();
    
    std::vector<int> nodes_to_add = my_setdiff(mut_indices, pathway);
    
    if (current_pos != 0) {
      auto to_remove = std::find(nodes_to_add.begin(), nodes_to_add.end(), nodes.at(current_pos));
      if (to_remove != nodes_to_add.end()) {
        nodes_to_add.erase(to_remove);
      }
    }
    
    for (int i {0}; i < ch_set_size; ++i) {
      std::vector<int> node = {nodes.at(current_ch_set.at(i) - 1)};
      nodes_to_add = my_setdiff(nodes_to_add, node);
    }
    // nodes_to_add = my_setdiff(nodes_to_add, current_ch_set);
    
    for (size_t i {0}; i < nodes_to_add.size(); ++i) {
      nodes.push_back(nodes_to_add.at(i));
      in_tree.push_back(false);
      IntegerVector ch {};
      children.push_back(ch);
      parents.push_back(current_pos + 1);
      current_ch_set.push_back(nodes.length());
    }
    children.at(current_pos) = current_ch_set;
    
    tree["nodes"] = nodes;
    tree["in_tree"] = in_tree;
    tree["children"] = children;
    tree["parents"] = parents;
    
    if (ch_set_size > 0) {
      if (current_pos != 0) {
        pathway.push_back(nodes.at(current_pos));
      }
      for (int pos : current_ch_set) {
        get_augmented_tree(mut_indices, tree, pathway, pos - 1); // c++ indexing starts from 0
      }
    }
  }
}

// [[Rcpp::export]]
List get_augmented_trees(int n, const List &trees) {
  
  List augmented_trees {};
  for (int i {0}; i < trees.length(); ++i) {
    List tree = trees.at(i);
    // // Remove duplicated nodes in the tree
    // int changed {0};
    // do {
    //   changed = 0;
    //   remove_duplicates(tree, changed);
    // } while (changed > 0);
    
    // Convert from T to A(T)
    std::vector<int> mut_indices(n);
    std::iota(mut_indices.begin(), mut_indices.end(), 1);
    get_augmented_tree(mut_indices, tree);
    
    augmented_trees.push_back(tree);
  }
  return augmented_trees;
}

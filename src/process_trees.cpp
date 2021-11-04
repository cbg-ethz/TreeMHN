#include <RcppArmadillo.h>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <regex>
#include <map>
#include <algorithm>
#include <iterator>
#include "utils.h"
using namespace Rcpp;
namespace fs = std::filesystem; // C++17
// [[Rcpp::depends(RcppArmadillo)]]

std::string trim_mut_name(std::string mut_name) {
  std::string new_name {};
  for (const char &s : mut_name) {
    if ((s == '_') || (s == '-')) {
      break;
    }
    new_name.push_back(s);
  }
  return new_name;
}

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
      } else {
        Rcout << nodes.at(current_pos) << " is not found in:" << "\n";
        std::cout << "nodes to add = { ";
        for (int n : nodes_to_add) {
          std::cout << n << ", ";
        }
        std::cout << "}; \n";
        std::cout << "pathway = { ";
        for (int n : pathway) {
          std::cout << n << ", ";
        }
        std::cout << "}; \n";
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
    // Remove duplicated nodes in the tree
    int changed {0};
    do {
      changed = 0;
      remove_duplicates(tree, changed);
    } while (changed > 0);

    // Convert from T to A(T)
    std::vector<int> mut_indices(n);
    std::iota(mut_indices.begin(), mut_indices.end(), 1);
    get_augmented_tree(mut_indices, tree);

    augmented_trees.push_back(tree);
  }
  return augmented_trees;
}

// [[Rcpp::export]]
List parse_trees(std::string path) {

  List trees {};
  int n {0}; // total number of mutations
  std::vector<std::string> mutations {};

  // Loop through all trees in the folder
  //std::regex fname_reg("AML-\\d{2,3}_fp_0\\.01-all_AML-\\d{2,3}-001\\.gv");
  std::regex fname_reg("AML-\\d{2,3}_c?_?fp_0\\.01_fn_0\\.\\d{1,6}-all_AML-\\d{2,3}-001_?c?\\.gv");
  std::regex node_reg("\\d\\[style=empty.+");
  std::regex mut_reg("\".+?\"");
  std::regex mut_idx_reg("\\d{1,4}\\[");
  std::regex edge_reg("\\d{1,3} -> \\d{1,3}");
  std::regex edge_reg_pa("\\d{1,3} ->");
  std::regex edge_reg_ch("-> \\d{1,3}");
  std::smatch m;

  for (const auto & entry : fs::directory_iterator(path)){
    std::ifstream in_file {entry.path()};
    std::string fname = entry.path().filename();

    if (std::regex_match(fname, m, fname_reg)) {

      String tree_name {};
      IntegerVector nodes {0};
      IntegerVector parents {1};
      List children {};
      LogicalVector in_tree {true};
      std::map<int, int> mut_map {};

      // Parse lines to obtain node and edge information
      std::string line {};
      while (std::getline(in_file, line)) {
        // Extract tree name
        if (line.substr(0,6) == "label=") {
          tree_name = std::string(&line[7], &line[line.size() - 2]);
        }

        // Get nodes
        if (std::regex_search(line, m, node_reg)) {
          int mut_idx; // mutation index in the tree
          std::string mut; // mutation name
          int idx {0}; // mutation index in all trees
          IntegerVector ch {}; // children set
          children.push_back(ch);
          // extract the mutation index in the tree
          if (std::regex_search(line, m, mut_idx_reg)) {
            mut = m.str();
            mut.pop_back();
            mut_idx = std::stoi(mut);
          }
          // extract the mutation name in the tree
          if (std::regex_search(line, m, mut_reg)) {
            mut = m.str();
            mut = trim_mut_name(std::string(&mut[1], &mut[mut.size() - 1]));
          }
          // determine the mutation index in all trees (except for the root)
          if (mut != "Root") {
            auto mut_it = std::find(mutations.begin(), mutations.end(), mut);
            if (mut_it != mutations.end()) {
              idx = mut_it - mutations.begin() + 1;
            } else {
              idx = n+1;
              n++;
              mutations.push_back(mut);
            }
            nodes.push_back(idx);
            in_tree.push_back(true);
            parents.push_back(0);
          }
          mut_map[mut_idx] = idx;
        }

        // Get edges
        if (std::regex_search(line, m, edge_reg)) {
          std::string edge = m.str();
          int i {0}, j {0};
          if (std::regex_search(edge, m, edge_reg_pa)) {
            std::string temp = m.str();
            i = std::stoi(std::string(&temp[0], &temp[temp.size() - 3]));
          }
          if (std::regex_search(edge, m, edge_reg_ch)) {
            std::string temp = m.str();
            j = std::stoi(std::string(&temp[3], &temp[temp.size()]));
          }

          if (mut_map[i] == 0) {
            parents[j+1] = 1;
            IntegerVector ch = children.at(0);
            ch.push_back(j+2);
            children.at(0) = ch;
          } else {
            parents[j+1] = i+2;
            IntegerVector ch = children.at(i+1);
            ch.push_back(j+2);
            children.at(i+1) = ch;
          }
        }
      }

      // Add the tree to list
      List tree = List::create(Named("tree_ID") = tree_name,
                               Named("nodes") = nodes,
                               Named("parents") = parents,
                               Named("children") = children,
                               Named("in_tree") = in_tree);
      trees.push_back(tree);
    }

    in_file.close();
  }

  // Convert trees (T) to augmented trees (A(T)) and remove repeated mutations
  List augmented_trees = get_augmented_trees(n, trees);

  List Tree_summary = List::create(Named("n") = n,
                                   Named("mutations") = mutations,
                                   Named("trees") = augmented_trees);
  return Tree_summary;

}

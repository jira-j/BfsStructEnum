#pragma once 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <array>
#include <valarray>
#include <string>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>
#include <tuple>

#include <time.h>
#include <curl/curl.h>

typedef int valence_value_type;
typedef int label_value_type; // `signed' is required because `-' is used.
typedef std::vector< std::vector<int> > carbon_position;
typedef std::vector< std::vector<int> > auto_group;

bool compare_function(std::vector<int>,std::vector<int>);

const int num_distinct_atoms = 4;
//std::string atomchar[] = {"C", "N", "O", "H"};
std::vector<std::string> atomchar = {"C", "N", "O", "H"};
const std::string bondchar[] = {"X", "", "=", "#"};
//int valence[] = {4, 3, 2, 1};
std::vector<int> valence = {4, 3, 2, 1};
const std::string input_atomchar[] = {"C", "N", "O", "H"};
const std::string input_bondchar[] = {"X", "", "=", "#"};
const int input_valence[] = {4, 3, 2, 1};
const int max_valence = 10;
int first_atom_valence_one = 3;
int num_special_atom =0;

static bool do_print = false;
static int num_except_H = 0;
static int num_H = 0;
static int num_lack_H = 0;
static int num_naph_bond = 0;

static std::vector<int> t_r, t_v;
static std::vector<valence_value_type> location;
static std::vector<valence_value_type> lack_degree;

typedef std::valarray<char> is_ident_type;
//typedef std::vector<char> is_ident_type;
class ChemTreeCenter {

	struct Node {
		label_value_type label;
	//bool is_identical;
		std::array<int, max_valence> children;
		valence_value_type num_children;
		int parent;
		valence_value_type multi; // subtree id during construction of simple trees
	        valence_value_type nth; // start from 0 not 1
	//int depth;
		std::array<int, max_valence+1> bond_position;
	//      std::vector<int> bond_position;
	};

	std::vector<int> rest_atoms;
	std::valarray<Node> nodes;
	int num_nodes;
#ifdef CUT_BY_NUM_H
	int num_H_to_be_added;
#endif
	//int deepest_head;

	void init() { // initialize nodes
	//nodes.resize(num_except_H);
		num_nodes = 0;
#ifdef CUT_BY_NUM_H
		num_H_to_be_added = 0;
#endif
	//deepest_head = 1;
	}
	inline void printseq() const;
	inline void printseq_single() const;
	inline void printsmi(const int i = 0) const;
	inline void printsmi_single(const int i = 0) const;
	inline void printmol(const std::string& filename) const;
       
public:

	ChemTreeCenter(const std::initializer_list<int>& _rest, const int _num_except_H) : rest_atoms(_rest), nodes(std::valarray<Node>(_num_except_H)) {
		init();
	}
	ChemTreeCenter(const std::vector<int>& _rest, const int _num_except_H) : rest_atoms(_rest), nodes(std::valarray<Node>(_num_except_H)) {
		init();
	}

	inline void update_identical(is_ident_type& is_ident, const int i) const; // when a new node i is added
	inline void update_identical_multi(is_ident_type& is_ident, const int i) const ; 
	inline void update_identical_end(is_ident_type& is_ident, const int i) const ; // when a child of a node i is recognized not to be added more
	inline bool add_root(const label_value_type atom_label);
	inline bool add_root_child(is_ident_type& is_ident, const label_value_type atom_label);
	inline bool add_node(is_ident_type& is_ident, const int parenti, const label_value_type atom_label);
	inline void del_root() {
		num_nodes = 0;
		++(rest_atoms[nodes[0].label]);
	}
	inline void del_last_node() {
		--num_nodes;
		--(nodes[nodes[num_nodes].parent].num_children);
		++(rest_atoms[nodes[num_nodes].label]); 
	}
	inline label_value_type begin_atom_label(const is_ident_type& is_ident, const int i) const;
	inline label_value_type begin_atom_label_root() const
	{
		const int nc = num_nodes - 1;
		if (nc > 0) {
			return nodes[nc].label;
		}
		return 0;
	}
	inline valence_value_type max_multi(const is_ident_type& is_ident, const int i) const;
	inline int is_normal(const int deepest_head) const;
	inline bool is_multi_normal(const int v) const;
	inline bool can_be_added(const int i) const
	{
		return (valence[nodes[i].label] - nodes[i].num_children - 1 > 0);
	}
	inline bool can_be_added_root() const
	{
		return (valence[nodes[0].label] > nodes[0].num_children);
	}

	inline bool remain(const label_value_type atom_label) const
	{
		return (rest_atoms[atom_label] > 0);
	}
	inline int starti() const
	{
		return nodes[num_nodes-1].parent;
	}
	inline int get_num_nodes() const
	{
		return num_nodes;
	}
  inline int get_num_child(int index){
    return nodes[index].num_children;
  }
  inline void reset_bond_position(size_t index){
    for(size_t i=0; i<max_valence; i++){
      nodes[index].bond_position[i] = 0;
    }
  }

	inline bool share_only_root(int i) const // This function is available only during construction of simple trees
	{
		return (nodes[i].multi != nodes[num_nodes - 1].multi);
#if 0
int j = num_nodes - 1;
while (i != j) {
i = nodes[i].parent;
j = nodes[j].parent;
}
return (i == 0);
#endif
}
inline valence_value_type get_center(const int deepest_head) const // This function is available only during construction of simple trees
{
	const valence_value_type subtree = nodes[deepest_head].multi;
	if (subtree != nodes[num_nodes - 1].multi) {
		return 0;
	} else {
		return subtree;
	}
#if 0
int i = deepest_head;
int j = num_nodes - 1;
while (i != j) {
i = nodes[i].parent;
j = nodes[j].parent;
}
if (i == 0) {
return 0;
}
int previ;
do {
previ = i;
i = nodes[i].parent;
} while (i != 0);

return previ;
#endif
}
inline void set_multi_bond(const int i, const valence_value_type multiple) 
{
	nodes[i].multi = multiple;
}
inline void fill_rest_single_bond(const int i)
{
	for (int j = i; j < num_nodes; ++j) {
		nodes[j].multi = 1;
	}
}
inline void calc_lack_degree() const
{
	lack_degree.clear();
	lack_degree.push_back(valence[nodes[0].label] - nodes[0].num_children);
	for(int i = 1; i < num_nodes; ++i){
		lack_degree.push_back(valence[nodes[i].label] - nodes[i].num_children - 1);
	}
}
inline int get_parent(const int i) const
{
	return nodes[i].parent;
}
inline void print() const;
inline void print_single() const;


//===========================================================
  bool adjacent(int index1, int index2){
    if(nodes[index1].parent == index2 || nodes[index2].parent == index1)
      return true;

    return false;
  }

		bool can_be_added_multi(int i){
			if(i>0){
				//if current node is benzene and parent node is pyridine -> true
			        if(nodes[i].label<num_special_atom && nodes[get_parent(i)].label<num_special_atom){
				        return true;
			        }
			        //current node is special atom and parent node is not benzene -> true
				//if(nodes[i].label<num_special_atom && nodes[get_parent(i)].label==0){
				//return true;
				//}
				//if current nodes is not special atom and parent is special atom -> false
				if(nodes[i].label>=num_special_atom && nodes[get_parent(i)].label>=num_special_atom){
					return true;
				}
			}
			return false;
		}
		int get_label(int i){
			return nodes[i].label;
		}
                int get_parent(int i){
		        return nodes[i].parent;
		}
		int get_multi(int i){
			return nodes[i].multi;
		}
                int get_nth(int i){
		  return nodes[i].nth;
		}
		void assign_bond_position(int index,std::vector< std::vector<int> >&valid_pos){
		  size_t start_j = 0;
		  size_t step = 1;
		  if(index != 0){
		    nodes[index].bond_position[0] = valid_pos[0][0];
		    start_j = 1;
		  }else{
		    nodes[index].bond_position[0] = -1;
		  }
		  for(size_t j=start_j; j < valid_pos.size(); j++){
		    for(size_t k=0; k < valid_pos[j].size(); k++){
		      //if( c_index[j][k] !=-1 ){
		      nodes[index].bond_position[step] = valid_pos[j][k];
		      step++;
		
		      /*if( c_index[j][k] > index ){ //c_index[j][k] is child node
			nodes[index].bond_position[nodes[c_index[j][k]].nth+1] = valid_pos[j][k];
			}else{
			  std::cout<<"hello"<<std::endl;
			  nodes[index].bond_position[0] = valid_pos[j][k];
			  }*/
			//}
		    }
		  }
		}
		void clear_bond_position(int i){
			for(size_t j=0;j<nodes[i].num_children;j++)
				nodes[i].bond_position[j] = -1;
		}

		bool is_equal(int,int);
		void generate_adj_list(std::vector< std::vector<int> >&,int);
		
		
                size_t label_substructure(const int index, const std::vector< std::vector<int> > & sympath_collection,  std::vector<carbon_position> &result,  const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi, std::ofstream &outputfile);
                size_t label_child_and_next(std::vector<carbon_position> &result,  int j,  int k,  int index,  const std::vector< std::vector<int> > & sympath_collection,  const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi, std::ofstream &outputfile);

                bool is_normal_sympath(const std::vector< std::vector<int> > & cp_list,  const std::vector<  std::vector<int> > & sympath_collection, const int index, const std::vector< auto_group > &aut_group);
		void find_sympath(std::vector< std::vector<int> > &chain);
                void find_sympath_down(int step,int index_step,std::vector<int> &temp,std::vector< std::vector<int> > &temp2);
		bool chain_middle_symmetry(std::vector<int> chain);
                int chain_symmetry(int index1,int index2);
                int chain_middle_pair_symmetry(const std::vector<int> &path, const carbon_position &cp_list, const int latest_unequal);
                bool is_tsub_equal(int index1, int index2);
                bool is_tsub_cp_equal(int index1, int index2, int path_index1, int path_index2, const std::vector<int> &path);
                bool is_tsub_has_substr(int index);
                bool chain_end_symmetry(int index1,int index2);
                bool is_chain_end_redundant(std::vector<int> chain, const carbon_position &carbon_pos);
                bool is_symmetry(int index);//find if this node is symmetry or not (all its child nodes must be assigned carbon_pos), use before is_normal_benzene
		bool is_trisymmetry(std::vector<int> chain1,std::vector<int> chain2);
                bool is_normal_str(const std::vector< std::vector<int> > &cp_list,  const std::vector< std::vector<int> > &aut_group);
		bool is_trisymmetry_redundant(int index_end1,int index_end2,const carbon_position & temp_pos);
		bool is_updown_symmetry(int index,int dealing_child,int first_unassigned_nth);
		bool is_fused_benzene(std::vector<int> chain){
			for(int i=0;i<chain.size()-1;i++){
				if(this->nodes[chain[i]].parent==chain[i+1]){
					if(this->nodes[chain[i]].multi!=2)
						return false;
				}
				if(this->nodes[chain[i+1]].parent==chain[i]){
					if(this->nodes[chain[i+1]].multi!=2)
						return false;
				}
			}
			return true;
		}

		void get_bond_position(carbon_position &cp_child,int index){
			for(int i=0;i<cp_child.size();i++){
				for(int j=0;j<cp_child[i].size();j++){
					if(cp_child[i][j]==-1){
						cp_child[i][j] = cp_child[i][j-1]+1;
					}else if(cp_child[i][j]>index){
						cp_child[i][j] = this->nodes[index].bond_position[nodes[cp_child[i][j]].nth+1];
					}else{
					  cp_child[i][j] = this->nodes[index].bond_position[0];
					}
				}
			}
		}
                void get_bond_position(std::vector<int> &cp_child,int index){
			for(int i=0; i<cp_child.size(); i++){
			  if(cp_child[i] > index){
			    cp_child[i] = this->nodes[index].bond_position[nodes[cp_child[i]].nth+1];
			  }else{
			    cp_child[i] = this->nodes[index].bond_position[0];
			  }
			}
		}
		void show() const{
			print_tree(0);
		}
void write_smiles(int index, std::ofstream & file, int & num_cycle, const std::vector<std::string> & str_smi) const;

void write(std::ofstream& outputfile, const std::vector<std::string> &str_smi) const{
    int num_cycle = 0;
    
    write_smiles(0, outputfile, num_cycle, str_smi);
  }


bool debug(){
  //return false;
  return true;
  return this->nodes[0].num_children==3 && this->nodes[1].label==0 && this->nodes[3].label==2 && this->nodes[2].multi==2;
}

bool code(){
  return false;//true;
}

		void print_tree(int index) const{
			Node n = nodes[index];
			
			int temp = index;
			while(temp != 0){
			        std::cout<<"        ";
				temp = nodes[temp].parent;
			}
			std::cout<<" label =" << atomchar[n.label] <<" bond = "<<n.multi;//<<(n.is_identical?" (identical) ":" ");
			if(index > 0)
			  std::cout<<" parent = "<<n.parent;
			if(index > 0 && nodes[n.parent].label<num_special_atom && nodes[n.parent].bond_position.size()>0){
				std::cout<<" c_position ="<<nodes[n.parent].bond_position[n.nth+1];	
			}
			if(index > 0 && n.label < num_special_atom){
			  std::cout<<"         parent_position ="<<n.bond_position[0];
			}
			std::cout<<std::endl;

			if(index == this->num_nodes-1){
				return;
			}
			if(n.num_children!=0){
				print_tree(n.children[0]);
			}

			if(nodes[index+1].parent == n.parent){
				print_tree(index+1);   
			}

		}
          void write_tree(int index,std::ofstream& output_file) const{
			Node n = nodes[index];
       
			int temp=index;
			while(temp!=0){
				output_file<< "        ";
				temp = nodes[temp].parent;
			}
			output_file<<" label =" <<atomchar[n.label] <<" bond = "<<n.multi;//<<(n.is_identical?" (identical) ":" ");
			output_file<<" parent = "<<n.parent;
			/*if(index>0 && nodes[n.parent].label == 0 && nodes[n.parent].bond_position.size()>0){
				output_file<<" c_position ="<<nodes[n.parent].bond_position[n.nth+1];
				if(n.multi ==2){
				  output_file<<","<<nodes[n.parent].bond_position[n.nth+1]+1;
				  output_file<<"         parent_position ="<<n.bond_position[0];
				}
				}*/
			output_file<<"\n";

			if(index == this->num_nodes-1){
				return;
			}
			if(n.num_children!=0){
			  write_tree(n.children[0],output_file);
			}

			if(nodes[index+1].parent == n.parent){
			  write_tree(index+1,output_file);   
			}

		}

	};

	inline void ChemTreeCenter::update_identical(is_ident_type& is_ident, const int i) const{
		using namespace std;

		location.clear();
		int j = i;
		int last_ident = -1;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; // because of BFS 
				BOOST_REVERSE_FOREACH(const auto& loc, location) {
					k = nodes[k].children[loc];
				}
				if (nodes[k].label != nodes[i].label) {
					is_ident[j] = (char)0;
				} else {
					last_ident = j;
				}
			}
			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}

		if (last_ident < 0) {
			last_ident = i;
			is_ident[last_ident] = (char)(-1);
		}
		last_ident = nodes[last_ident].parent;
		while (is_ident[last_ident] == (char)0) {
			is_ident[last_ident] = (char)(-1);
			last_ident = nodes[last_ident].parent;
		}

	};



	inline void ChemTreeCenter::update_identical_end(is_ident_type& is_ident, const int i) const
	{
		using namespace std;

		location.clear();
		int j = i;
		int last_ident = -1;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					k = nodes[k].children[loc];
				}
				if(nodes[k].num_children != nodes[i].num_children) {
					is_ident[j] = (char)0;
				} else {
					last_ident = j;
				}
			}

			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}
		
		if (last_ident < 0) {
		  last_ident = i;
		  is_ident[last_ident] = (char)(-1);
		  //is_ident[i]=(char)(-1);
		}
		
		last_ident = nodes[last_ident].parent; 
		
		while (is_ident[last_ident] == (char)0  && last_ident > 0) { // add && last_ident > 0 :1/3/2016
		  
			is_ident[last_ident] = (char)(-1);
			last_ident = nodes[last_ident].parent; 
		}
	}

	inline bool ChemTreeCenter::add_root(const label_value_type atom_label)
	{ 
		auto& root = nodes[0];
		root.label = atom_label;
		root.num_children = 0;
		root.parent = -1;
		root.multi = 1;
//root.is_identical = false;
//is_ident[0] = false;
//root.nth = 0;
//root.depth = 0;
		num_nodes = 1;
		--(rest_atoms[atom_label]); 
		return num_nodes == num_except_H;
	}

	inline bool ChemTreeCenter::add_root_child(is_ident_type& is_ident, const label_value_type atom_label)
	{
		auto& node = nodes[num_nodes];
		auto& parnode = nodes[0];
		node.label = atom_label;
		const valence_value_type nc = parnode.num_children;
		node.parent = 0;
		node.multi = num_nodes;
		node.nth = nc;
		node.num_children = 0;

		parnode.children[nc] = num_nodes;
		++(parnode.num_children);

		if ((nc > 0) and (nodes[nc].label == atom_label)) {
			is_ident[num_nodes] = (char)1;
		} else {
			is_ident[num_nodes] = (char)(-1);
		}

		--(rest_atoms[atom_label]);
		++num_nodes;
		return num_nodes == num_except_H;
	}

	inline bool ChemTreeCenter::add_node(is_ident_type& is_ident, const int parenti, const label_value_type atom_label)
	{
//std::cerr << parenti << " " << num_nodes << " " << atomchar[atom_label] << std::endl;
		auto& node = nodes[num_nodes];
		auto& parnode = nodes[parenti];
		node.label = atom_label;
		const valence_value_type nc = parnode.num_children;
		node.parent = parenti;
		if (parenti == 0) { // subtree id
			node.multi = num_nodes;
		} else {
			node.multi = parnode.multi;
		}
		node.nth = nc;
		node.num_children = 0;

		parnode.children[nc] = num_nodes;
		++(parnode.num_children);

		is_ident[num_nodes] = (char)(nc > 0);
		update_identical(is_ident, num_nodes);

		--(rest_atoms[atom_label]);
		++num_nodes;

		return num_nodes == num_except_H;
	}

	inline label_value_type ChemTreeCenter::begin_atom_label(const is_ident_type& is_ident, const int i) const
	{
		using namespace std;

		int begin_atom = 0;
		const valence_value_type nc = nodes[i].num_children;
		if (nc > 0) {
			begin_atom = nodes[num_nodes-1].label;
		}
		location.clear();
		location.push_back(nc);
		int j = i;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 

				bool flag = false;
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					if(nodes[k].num_children > loc) {
						k = nodes[k].children[loc];
					} else {
						begin_atom = first_atom_valence_one;
						flag = true;
						break;
					}
				}
				if(flag) break;
				else{
					if (nodes[k].label > begin_atom) {
						begin_atom = nodes[k].label;
					}
				}
			}
			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}

		return begin_atom;
	}

	inline int ChemTreeCenter::is_normal(const int deepest_head) const
	{
		using namespace std;

		const valence_value_type v = get_center(deepest_head);

		if (v == 0) {
			return 1;
		}

		{ // root
			const label_value_type labelrv = nodes[0].label - nodes[v].label;
			if (labelrv != 0) {
				return (int)(labelrv > 0);
			}
		}

// children of root
		t_r.clear();
		t_v.clear();
		valence_value_type ncv = nodes[v].num_children + 1;
		const valence_value_type minv = (v < ncv) ? v : ncv;
		for (int i = 1; i < minv; ++i) {
			const int vc = nodes[v].children[i-1];
			const label_value_type labelrv = nodes[i].label - nodes[vc].label;
			if (labelrv != 0) {
				return (int)(labelrv > 0);
			}
			t_r.push_back(i);
			t_v.push_back(vc);
		}
		const valence_value_type nc0 = nodes[0].num_children;
		const valence_value_type minnc = (nc0 < ncv) ? nc0 : ncv;
		for (int i = v + 1; i <= minnc; ++i) {
			const int vc = nodes[v].children[i-2];
			const label_value_type labelrv = nodes[i].label - nodes[vc].label;
			if (labelrv != 0) {
				return (int)(labelrv > 0);
			}
			t_r.push_back(i);
			t_v.push_back(vc);
		}
		ncv -= nc0;
		if (ncv != 0) {
			return (int)(ncv > 0);
		}

		size_t j = 0;
		while (j < t_r.size()) {
			const auto& nr = nodes[t_r[j]];
			const auto& nv = nodes[t_v[j]];
			const valence_value_type ch_r = nr.num_children;
			const valence_value_type ch_v = nv.num_children;
			const valence_value_type minnc = (ch_r < ch_v) ? ch_r : ch_v;
			for(int i = 0;i < minnc;i++){
				const int nri = nr.children[i];
				const int nvi = nv.children[i];
				const label_value_type labelrv = nodes[nri].label - nodes[nvi].label;
				if (labelrv != 0) {
					return (int)(labelrv > 0);
				}
				t_r.push_back(nri);
				t_v.push_back(nvi);
			}
			if (ch_r != ch_v) {
				return (int)(ch_r < ch_v);
			}
			j++;
		}
		return 2+v;
	}

	inline bool ChemTreeCenter::is_multi_normal(const int v) const
	{
		using namespace std;

		t_r.clear();
		t_v.clear();
//cerr << v << endl;
		const int r = 0;
		for(int i = 1; i < v; ++i) {
			t_r.push_back(i);
		}
		for (int i = v + 1; i <= nodes[r].num_children; ++i){
			t_r.push_back(i);
		}
		for(int i =0;i<nodes[v].num_children;++i){
			t_v.push_back(nodes[v].children[i]);
		}
		size_t j = 0;
		while (j < t_r.size()) {
			const auto& nr = nodes[t_r[j]];
			const auto& nv = nodes[t_v[j]];
			const valence_value_type mvr = nv.multi - nr.multi;
			if (mvr != 0) {
				return (mvr > 0);
			} 
			for(int i = 0; i < nr.num_children; ++i){
				t_r.push_back(nr.children[i]);
				t_v.push_back(nv.children[i]);
			}
			++j;
		}

		return true;
	}

	inline void ChemTreeCenter::update_identical_multi(is_ident_type& is_ident, const int i) const
	{
		using namespace std;

		location.clear();
		int j = i;
		int last_ident = -1;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					k = nodes[k].children[loc];
				}
				if(nodes[k].multi != nodes[i].multi) {
					is_ident[j]= (char)0;
				} else {
					last_ident = j;
				}
			}
			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}

		if (last_ident < 0) {
			last_ident = i;
			is_ident[last_ident] = (char)(-1);
		}
		last_ident = nodes[last_ident].parent;
		while (is_ident[last_ident] == (char)0) {
			is_ident[last_ident] = (char)(-1);
			last_ident = nodes[last_ident].parent;
		}
	}

	inline valence_value_type ChemTreeCenter::max_multi(const is_ident_type& is_ident, const int i) const
	{
		valence_value_type maxmulti = max_valence;
		location.clear();
		int j = i;
		while (is_ident[j] >= (char)0) {
			if (is_ident[j] > (char)0) {
				int k = j - 1; 
				BOOST_REVERSE_FOREACH(const valence_value_type loc, location) {
					k = nodes[k].children[loc];
				}
				const valence_value_type kmulti = nodes[k].multi;
				if(kmulti < maxmulti)
				{
					maxmulti = kmulti;
					if (maxmulti == 1) break;
				}
			}

			location.push_back(nodes[j].nth);
			j = nodes[j].parent;
		}
		return maxmulti;
	}

//===aut function=====
void print_cp(const std::vector< std::vector<int> > &cp_list){
  for(size_t i = 0; i<cp_list.size(); i++){
    std::cout<<"[";
    for(size_t j = 0; j<cp_list[i].size(); j++){
      std::cout<<cp_list[i][j];
      if(j+1<cp_list[i].size())
	std::cout<<" ";
    }
    std::cout<<"]";
  }
  std::cout<<std::endl;
}

void print_autgroup(const auto_group &group){
  for(const auto&mapping : group) {
    for(size_t i=0; i < mapping.size(); ++i) {
      std::cout<<mapping[i]<<"\t";
    }
    std::cout<<std::endl;
  }	    
}

float get_degree(int bond_type){
  switch(bond_type){
  case 1:
  case 2:
  case 3:
    return static_cast<float>(bond_type);
    break;
  case 4:
  case 5:
    return 1.5;
    break;
  }
  return -1;
}

class Mol {
  std::string atomlabels;
  std::vector<std::tuple<int,int,int>> edges;
  std::vector<std::vector<std::pair<int,int>>> adjedges;
  boost::multi_array<int, 2> bondmat;

public:
  void read(const std::string& filename);
  void autsub(std::vector<std::vector<int>>& group, std::vector<int>& mapping, 
	      const std::vector<std::pair<int, int>>& searchorder, const std::vector<int>& restedges, 
	      std::vector<bool>& dvisit, const size_t si) const;
  void aut(std::vector<std::vector<int>>& group) const;
  std::string get_atomlabels(){
    return atomlabels;
  }
  void get_edges(std::vector< std::tuple<int,int,int> > &get_edge){
    get_edge = edges;
  }

  int get_numcycles();
  void find_cycle(std::vector< std::vector<int> > &path_collection);
  void remove_redundant_cycle(std::vector< std::vector<int> > &cycle_collection);
  std::valarray<int> get_lack_valence();
};

void Mol::read(const std::string& filename) 
{
  using namespace std;

  ifstream ifs(filename);

  int num_atoms,  num_bonds;
  string buffer;
  size_t num_header = 3;
  for(size_t i = 0; i < num_header; i++){
    getline(ifs, buffer);
  }
  ifs >> num_atoms >> num_bonds;
  getline(ifs, buffer);
  
  atomlabels.resize(num_atoms);
  float x,y,z;
  for (int i = 0; i != num_atoms; ++i) {
    ifs >> x >> y >> z >> atomlabels[i];
    getline(ifs, buffer);
  }
  edges.resize(num_bonds);
  adjedges.resize(num_atoms);
  bondmat.resize(boost::extents[num_atoms][num_atoms]);
  for (int ei = 0; ei != num_bonds; ++ei) {
    int vi, vj, bondtype;
    ifs >> vi >> vj >> bondtype;
    vi--;
    vj--; //vi and vj is in [1,num_atoms] change it to [0,num_atoms-1]
    getline(ifs, buffer);
    edges[ei] = make_tuple(vi, vj, bondtype);
    adjedges[vi].push_back(make_pair(ei, vj));
    adjedges[vj].push_back(make_pair(ei, vi));
    bondmat[vi][vj] = bondtype;
    bondmat[vj][vi] = bondtype;
    if (vi == vj) {
      cerr << "error: same node\n";
    }
  }
}

void Mol::autsub(std::vector<std::vector<int>>& group, std::vector<int>& mapping, 
		 const std::vector<std::pair<int, int>>& searchorder, const std::vector<int>& restedges, 
		 std::vector<bool>& dvisit, const size_t si) const
{
  int svi, svj;
  std::tie(svi, svj) = searchorder[si];
  const int vi = mapping[svi];
  const int bondtype = bondmat[svi][svj];
  for (const auto& ejt : adjedges[vi]) {
    int ej, vj;
    std::tie(ej, vj) = ejt;
    if (dvisit[vj]) {
      if ((bondtype == bondmat[vi][vj]) and (atomlabels[svj] == atomlabels[vj])) {
        mapping[svj] = vj;
        if (si+1 == searchorder.size()) {
          bool valid = true;
          for (const auto& re : restedges) {
            int ri, rj, rb;
	    std::tie(ri, rj, rb) = edges[re];
            const int mi = mapping[ri];
            const int mj = mapping[rj];
            if (bondmat[mi][mj] != rb) {
              valid = false;
              break;
            }
          }
          if (valid) {
            group.push_back(mapping);
          }
        } else {
          dvisit[vj] = false;
          autsub(group, mapping, searchorder, restedges, dvisit, si+1);
          dvisit[vj] = true;
        }
      }
    }
  }
}

void Mol::aut(std::vector<std::vector<int>>& group) const
{
  using namespace std;
  vector<pair<int, int>> searchorder;
  vector<int> restedges;
  vector<bool> visit(atomlabels.size(), true);
  vector<bool> src(edges.size(), true); 
  int vi = 0;
  visit[vi].flip();
  size_t cnt = 1;

  while (cnt < atomlabels.size()) {
    bool deadend = true;
    for (const auto& ae : adjedges[vi]) {
      if (visit[ae.second]) {
        src[ae.first] = false;
        visit[ae.second] = false;
        searchorder.push_back(make_pair(vi, ae.second));
        ++cnt;
        vi = ae.second;
        deadend = false;
        break;
      }
    }
    if (deadend) {
      for (size_t i = 0; i != edges.size(); ++i) {
        int ix, iy, ib;
	std::tie(ix, iy, ib) = edges[i];
        if (visit[ix] and (not visit[iy])) {
          vi = iy;
          break;
        } else if ((not visit[ix]) and visit[iy]) {
          vi = ix;
          break;
        }
      }
    }
  }

  for (size_t i = 0; i != edges.size(); ++i) {
    if (src[i]) {
      restedges.push_back(i);
    }
  }

  vector<int> mapping(atomlabels.size()); 
  for (size_t i = 0; i != atomlabels.size(); ++i) {
    if (atomlabels[0] == atomlabels[i]) {
      mapping[0] = i;
      vector<bool> dvisit(atomlabels.size(), true);
      dvisit[i] = false;
      autsub(group, mapping, searchorder, restedges, dvisit, 0);
    }
  }
}

int Mol::get_numcycles(){
  std::vector< std::vector<int> > cycle_collection;
  
  for(size_t i=0; i<atomlabels.size(); i++){
    std::vector<int> temp_path(1,i);
    cycle_collection.push_back(temp_path);
    find_cycle(cycle_collection);    
  }

  remove_redundant_cycle(cycle_collection);

  return cycle_collection.size();
}

void Mol::find_cycle(std::vector< std::vector<int> > &path_collection){

  for(size_t i=0; i<path_collection.size(); i++){

    int last_index = path_collection[i].size()-1;
    if(path_collection[i][last_index] >= 0){
      //path is not finish yet -> find the next edge
      std::vector<int> temp_path = path_collection[i];
      path_collection.erase(path_collection.begin()+i);
      i--;
      
      bool end_node = true; 

      for(size_t j=0; j<edges.size(); j++){
	int next_node = -1;
	if(std::get<0>(edges[j])==temp_path[last_index]){
	  next_node = std::get<1>(edges[j]);
	}else if(std::get<1>(edges[j])==temp_path[last_index]){
	  next_node = std::get<0>(edges[j]);
	}
	if(next_node!=-1){
	  if(last_index == 0 || next_node != temp_path[last_index-1]){	    	   
	    temp_path.push_back(next_node);
	  
	    std::vector<int>::iterator start_cycle = std::find(temp_path.begin(),temp_path.begin()+last_index,next_node);
	    //start_cycle = iterator point to index of the same node as the last node
	    if(start_cycle != temp_path.begin()+last_index){
	      //if we found start_cycle before last node -> this path contains a cycle	      
	      std::vector<int> full_path = temp_path;
	      temp_path.erase(temp_path.begin(),start_cycle); //remove the part outside a cycle

	      //if it contains a cycle, then the path is finish -> push -1 to the last element to indicate that this path is finish
	      temp_path.push_back(-1);
	      path_collection.push_back(temp_path);
	      find_cycle(path_collection);

	      temp_path = full_path; //restore the same temp_path
	      temp_path.pop_back();
	    }else{	    
	      path_collection.push_back(temp_path);
	      find_cycle(path_collection);
	   
	      temp_path.pop_back();	      
	    }
	    end_node = false;
	  }
	}
      }//end for edges

      if(end_node){
	temp_path.push_back(-2);
      }

    }
  }
}

void Mol::remove_redundant_cycle(std::vector< std::vector<int> > &cycle_collection){
  for(size_t i=0; i<cycle_collection.size(); i++){
    //delete -1 and the start point of cycle at the end of the path (start point still at the begin of the path)
    cycle_collection[i].erase(cycle_collection[i].end()-2,cycle_collection[i].end());
    
    std::vector<int>::iterator min_index = std::min_element(cycle_collection[i].begin(),cycle_collection[i].end()); 

    std::vector<int> temp_cycle = cycle_collection[i];
    if(min_index!=cycle_collection[i].begin()){
      //if min element is not the first element, shift the cycle such that min element become the first
      int range = cycle_collection[i].end()-min_index;

      std::move(min_index,cycle_collection[i].end(),temp_cycle.begin());
      std::move(cycle_collection[i].begin(),min_index,temp_cycle.begin()+range);
    }

    if(temp_cycle[1]>temp_cycle[temp_cycle.size()-1]){
      //if sequence of node in reverse order is less than the original order, reverse it
      std::reverse(temp_cycle.begin()+1,temp_cycle.end());    
    }
    cycle_collection[i] = temp_cycle;

    if(std::find(cycle_collection.begin(),cycle_collection.begin()+i,temp_cycle)!=cycle_collection.begin()+i){
      cycle_collection.erase(cycle_collection.begin()+i);
      i--;
      continue;
    }
  }
}


std::valarray<int> Mol::get_lack_valence(){

  std::valarray<int> result(atomlabels.size());

  for(size_t i=0; i<atomlabels.size(); i++){
    float sum_bond=0; //summation of its current degree

    for(size_t j=0; j<edges.size(); j++){
      if(std::get<0>(edges[j]) == i || std::get<1>(edges[j]) == i ){
	float degree = get_degree(std::get<2>(edges[j]) );
	if(degree!=-1){
	  sum_bond += degree;
	}
      }
    }
    
    std::string temp = std::string(1,atomlabels[i]);
    int index = std::find(input_atomchar,input_atomchar+num_distinct_atoms, temp) - input_atomchar;

    if(index < num_distinct_atoms){ 
      result[i] = input_valence[index]-sum_bond; //lack_valence = its valence - current degree
    }else{
      std::cout<<"Error: atoms in input structure are not appear in input atoms"<<std::endl;
    }

  }
  
  return result;
}

void get_automorphism(std::vector< std::vector<int> > &position, const std::vector<int> & aut_group){
  
  for(size_t i=0; i<position.size(); i++){ 
    for(size_t j=0; j<position[i].size(); j++){
      if(position[i][j]==-1){
	std::sort(position[i].begin(),position[i].begin()+j);
	return;
      }
      position[i][j] = aut_group[ position[i][j] ];
    }
    std::sort(position[i].begin(),position[i].end());
  }        

  /*for(size_t i=0; i<position.size(); i++){ 
    for(size_t j=0; j<position[i].size(); j++){
      std::cout<<position[i][j]<<" ";
    }
  }  
  std::cout<<std::endl;
  */
}
//====================logP value
/*
struct string {
  char *ptr;
  size_t len;
};

void init_string(struct string *s) {
  s->len = 0;
  s->ptr = malloc(s->len+1);
  if (s->ptr == NULL) {
    fprintf(stderr, "malloc() failed\n");
    exit(EXIT_FAILURE);
  }
  s->ptr[0] = '\0';
}

size_t writefunc(void *ptr, size_t size, size_t nmemb, struct string *s)
{
  size_t new_len = s->len + size*nmemb;
  s->ptr = realloc(s->ptr, new_len+1);
  if (s->ptr == NULL) {
    fprintf(stderr, "realloc() failed\n");
    exit(EXIT_FAILURE);
  }
  memcpy(s->ptr+s->len, ptr, size*nmemb);
  s->ptr[new_len] = '\0';
  s->len = new_len;

  return size*nmemb;
}

void get_logp(std::string smiles){
  using namespace std;
  CURL *curl;
  CURLcode res;
  

  curl = curl_easy_init();
  string url = "http://www.vcclab.org/web/alogps/calc?SMILES="+smiles;
  
  if(curl) {
    struct string logp;
    init_string(&logp);

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str() );
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writefunc);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &logp);
    
    res = curl_easy_perform(curl);

    printf("%s\n", s.ptr);
    free(s.ptr);

    curl_easy_cleanup(curl);
  } 
}
*/

//====================

void count_structure(std::vector<int> &str_atom_number, int str_bond_number, std::vector< std::vector<int> > &input_atom_set, std::vector<int> &bond_set, int index_structure){
  int initial_set_number = input_atom_set.size();

  for(size_t set_i=0; set_i < initial_set_number; set_i++){
    std::vector<int> temp_input_atom(input_atom_set[set_i]);
    int temp_bond = bond_set[set_i];

    bool is_atom_enough = true;

    for(size_t i=0; i<temp_input_atom.size(); i++){
      if(temp_input_atom[i] < str_atom_number[i]){
	is_atom_enough = false;
	break;
      }
    }
    if(temp_bond < str_bond_number) std::cout<<"dou not enough    chemical formula:"<<temp_bond<<"  naph"<<str_bond_number<<std::endl;
    while(is_atom_enough && temp_bond >= str_bond_number){ 
      //plus amount of substructure by one
      temp_input_atom[index_structure]++;
      
      for(size_t i=0; i<temp_input_atom.size(); i++){
	temp_input_atom[i] -= str_atom_number[i];
	if(temp_input_atom[i] < str_atom_number[i]){
	  is_atom_enough = false;
	}
      }

      temp_bond -= str_bond_number;

      input_atom_set.push_back(temp_input_atom);
      bond_set.push_back(temp_bond);
    }
  }

}

void ChemTreeCenter::generate_adj_list(std::vector< std::vector<int> > &adj_list, int index){
  //value in c_index is the index of that node ex. [ [0] [3 4] ] means node[0], node[3] and node[4] adjacent with this node and tsub(node[3]) = tsub(node[4])
	
	if(index!=0){
	  std::vector<int> temp_list;
	  temp_list.push_back(nodes[index].parent);
	  if(nodes[index].multi==2){
	    temp_list.push_back(-1);
	  }
	  adj_list.push_back(temp_list);
	}

	std::vector<int> child_list;
	for(size_t j=0; j<nodes[index].num_children; j++){
	  child_list.push_back( nodes[index].children[static_cast<int>(j)] );
	}

	int num_skip=0;
	for(size_t j=0; j < child_list.size(); j++){
	  num_skip=0;
	  std::vector<int> temp(1, child_list[j+num_skip]);
	  if(nodes[child_list[j+num_skip]].multi==2){
	    temp.push_back(-1);
	  }

	  for(size_t k=j+1; k < child_list.size(); k++){
	    if(this->is_equal(child_list[j], child_list[k])){
	      temp.push_back(child_list[k]);
	      num_skip++;
	      if(nodes[child_list[k]].multi==2){
		temp.push_back(-1);
	      }
	      child_list.erase(child_list.begin()+k);
	      k--;
	    }else{
	      break;
	    }
	  }

	  adj_list.push_back(temp);
	}
}

bool ChemTreeCenter::is_equal(int i,int j){
	if(nodes[i].label != nodes[j].label || nodes[i].multi!=nodes[j].multi || nodes[i].num_children != nodes[j].num_children)
		return false;
	if(nodes[i].label >= num_special_atom && nodes[i].num_children == 0)
		return true;
	if(nodes[i].label < num_special_atom && nodes[i].bond_position != nodes[j].bond_position)
	  return false;

	//bool result = is_equal(nodes[i].children[0],nodes[j].children[0]); 
	for(size_t index=0;index<nodes[i].num_children;index++){
	  if(!is_equal(nodes[i].children[index],nodes[j].children[index])){
	    return false;
	  }
	  //result = result & is_equal(nodes[i].children[index],nodes[j].children[index]);
	}
	return true;
	//return result;
}

int max(std::vector<int> container){
	int result=container[0];
	for(int i=1;i<container.size();i++){
		if(container[i]>result){
			result = container[i];
		}
	}
	if(result<0)
	   result = 0;
	return result;
}

bool found(const std::vector<int> &container,int val){
  //return T if val is found in container, otherwise return F
  return std::find(container.begin(), container.end(), val) != container.end();
}

bool is_cp_available(const std::vector< std::vector<int> > &cp_position, const std::valarray<int> &str_lack_valence){
  
  std::valarray<int> cp_valence(str_lack_valence);

  for(size_t i=0; i<cp_position.size(); i++){
    for(size_t j=0; j<cp_position[i].size(); j++){
      
      if(cp_position[i][j] != -1){
	cp_valence[cp_position[i][j]]--;
      }else{
	return true;
      }

      if(cp_valence[cp_position[i][j]]<0){
	return false;
      }
    }
  }

  return true;
}


bool ChemTreeCenter::is_trisymmetry(std::vector<int> chain1,std::vector<int> chain2){
  if(!is_fused_benzene(chain1) || !is_fused_benzene(chain2)){
    return false;
  }
  if(chain1.size()%2==0){
    return false;
  }

  int center_node = (chain1.size())/2;
  
  if(chain1[center_node]!=chain2[center_node]){
    return false;
  }
  if(chain1[center_node]>chain1[center_node+1] || chain2[center_node]>chain2[center_node+1]){
    return false;
  }

  for(int i=1;i<=(chain1.size()/2)-1;i++){
    Node parent_benzene[3] = {this->nodes[chain1[center_node-i]],this->nodes[chain1[center_node+i]],this->nodes[chain2[center_node+i]]};
    Node child_benzene[3] = {this->nodes[chain1[center_node-i-1]],this->nodes[chain1[center_node+1+i]],this->nodes[chain2[center_node+i+1]]};
    //check for carbon position of each node in the chain, if they are not the same one -> no tri radical symmetry
    if(parent_benzene[0].bond_position[child_benzene[0].nth]!=parent_benzene[1].bond_position[child_benzene[1].nth])
      return false;
    if(parent_benzene[0].bond_position[child_benzene[0].nth]!=parent_benzene[2].bond_position[child_benzene[2].nth])
      return false;
  }
  return true;
}

void flip_trisymmetry(carbon_position (&flip_cp)[3],int i,int range){
  for(int k=0;k<3;k++){
    for(int j=0;j<range;j++){
      flip_cp[k][i][j] = 5-flip_cp[k][i][j];
    } 
    std::sort(flip_cp[k][i].begin(),flip_cp[k][i].begin()+range);
  }
}

void rotate_trisymmetry(carbon_position (&flip_cp)[3],int i){
  std::vector<int> temp = flip_cp[0][i];
  flip_cp[0][i] = flip_cp[1][i];
  flip_cp[1][i] = flip_cp[2][i];
  flip_cp[2][i] = temp;
}

bool greater_than(const carbon_position (&original_cp)[3],const carbon_position (&flip_cp)[3],int i,int range){
  for(int k=0;k<3;k++){
    for(int j=0;j<range;j++){
      if(original_cp[k][i][j]>flip_cp[k][i][j])
	return true;
      if(original_cp[k][i][j]<flip_cp[k][i][j])
	return false;
    }
  }
  return false;
}

bool ChemTreeCenter::is_trisymmetry_redundant(int index_end1,int index_end2,const carbon_position & temp_pos){
  carbon_position cp_end1,cp_end2;
  generate_adj_list(cp_end1,index_end1);
  //std::sort(cp_end1.begin(),cp_end1.end(),compare_function);
  generate_adj_list(cp_end2,index_end2);
  //std::sort(cp_end2.begin(),cp_end2.end(),compare_function);

  get_bond_position(cp_end1,index_end1);
  get_bond_position(cp_end2,index_end2);

  carbon_position original_cp[3] = {cp_end1,cp_end2,temp_pos};
  carbon_position flip_cp[3] = {cp_end1,cp_end2,temp_pos};

  for(int i=0;i<temp_pos.size();i++){
    int range = std::find(temp_pos[i].begin(),temp_pos[i].end(),-2)-temp_pos[i].begin();
    if(range >0){
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      
      rotate_trisymmetry(flip_cp,i); //with third rotation, flip_cp is back to original_cp 
      flip_trisymmetry(flip_cp,i,range);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
      rotate_trisymmetry(flip_cp,i);
      if(greater_than(original_cp,flip_cp,i,range)){
	return true;
      }
    }else{
      break;
    }
  }  
  return false;
};

bool ChemTreeCenter::is_updown_symmetry(int index,int dealing_adjacent,int first_unassigned_index){
  //dealing_adjacent is index (not nth) of adjacent node that is filling carbon_position, so we will not consider it in finding symmetry of next benzene
  //first call dealing_adjacent is child node but can be parent node in the future
  
  carbon_position c_index;
  generate_adj_list(c_index,index);

  //std::sort(c_index.begin(),c_index.end(),compare_function);
  

  if(this->get_label(dealing_adjacent)>=num_special_atom || this->nodes[dealing_adjacent].multi==1){
    if(dealing_adjacent>index){
      return true;
    }
  }

  int cp_dealing_adjacent1 = 0;
  //position of carbon in benzene index that bond with dealing_adjacent
  if(nodes[dealing_adjacent].parent==index){
    cp_dealing_adjacent1 = this->nodes[index].bond_position[this->nodes[dealing_adjacent].nth+1];
  }else{
    cp_dealing_adjacent1 = this->nodes[index].bond_position[0];
  }
  int cp_dealing_adjacent2 = cp_dealing_adjacent1+1;

  carbon_position cp(c_index),flip_cp(c_index);
  for(int i=0;i<c_index.size();i++){
    for(int j=0;j<c_index[i].size();j++){
      if(c_index[i][j]==-1){
	cp[i][j] = 0;
      }else{
	if(c_index[i][j]>index){
	  cp[i][j] = this->nodes[index].bond_position[nodes[c_index[i][j]].nth+1];
	}else{
	  cp[i][j] = this->nodes[index].bond_position[0];
	}
      }
      flip_cp[i][j] = (6+cp_dealing_adjacent1+cp_dealing_adjacent2-cp[i][j])%6;
    }
    std::sort(cp[i].begin(),cp[i].end());
    std::sort(flip_cp[i].begin(),flip_cp[i].end());
  }
 
  if(cp!=flip_cp){
    return false;
  }else{
    //check if benzene node 'b' bonding with current node is symmetry or not
    //if b is not symmetry, then current node is not symmetry
    for(int i=0;i<c_index.size();i++){
      for(int j=0;j<c_index[i].size();j++){
	if(this->get_label(c_index[i][j])<num_special_atom){
	  if(this->get_label(c_index[i][j])!=0){
	    return false; //if it is pyridine, it cannot be updown symmetry
	  }
	  if(c_index[i][j]==nodes[index].parent && nodes[dealing_adjacent].parent==index && nodes[index].multi==2){//c_index is parent
	    //int index_of_parent =this->nodes[index].parent;
	    if(!is_updown_symmetry(c_index[i][j],index,first_unassigned_index)){
	      return false;
	    }
	  }else if(c_index[i][j]!=-1 && c_index[i][j]<first_unassigned_index && c_index[i][j]!=this->nodes[dealing_adjacent].nth){//c_index is child that is not reserve carbon for naph_bond
	    //int index_of_child = this->nodes[index].children[c_index[i][j]];
	    if(this->nodes[c_index[i][j]].multi==2){
	      if(!is_updown_symmetry( c_index[i][j], index, first_unassigned_index)){
		return false;
	      }
	    }
	  }
	}
      }
    }
    return true;
  }
}

int ChemTreeCenter::chain_symmetry(int index1,int index2){
  //-1 if cp[index1] > cp[index2], where index1 < index2
  // 1 if cp[index1] < cp[index2], where index1 < index2
  // 0 if cp[index1] = cp[index2], where index1 < index2	     

  if(index1 > index2){
    int temp = index1;
    index1 = index2;
    index2 = temp;
  }

  if(this->get_label(index1)!=this->get_label(index2)){
    std::cout<<"Error chain_symmetry (1)"<<std::endl;
    return -2;
  }

  /*
  std::vector<int> children1,children2;
  for(int i=0;i<this->nodes[index1].num_children;i++){
    children1.push_back(this->nodes[index1].children[i]);
  }
  for(int i=0;i<this->nodes[index2].num_children;i++){
    children2.push_back(this->nodes[index2].children[i]);
  }
  */
  bool is_adjacent = false;
  int index_parent,index_child;
  
  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent==index2){
    std::cout<<"Error nodes cannot be adjacent"<<std::endl;
    is_adjacent = true;
    index_parent = index2;
    index_child = index1;
    //children2.erase(children2.begin()+this->nodes[index1].nth);
  }else if(this->nodes[index2].parent==index1){
    std::cout<<"Error nodes cannot be adjacent"<<std::endl;
    is_adjacent = true;
    index_parent = index1;
    index_child = index2;
    //children1.erase(children1.begin()+this->nodes[index2].nth);
  }

  if(!is_adjacent){
    if(this->get_label(index1) < num_special_atom){
      if(this->nodes[index1].bond_position == this->nodes[index2].bond_position)
	return 0; 
      if(this->nodes[index1].bond_position > this->nodes[index2].bond_position)
	return -1;
      if(this->nodes[index1].bond_position < this->nodes[index2].bond_position)
	return 1;
    }
  }else{
    std::cout<<"error in chain_symmetry"<<std::endl;
    return 0;
  }

  /*
  //parent must be root unless the tree is not left heavy
  if(index_parent!=0){
    return false;
  }
  //if it is not the first child of this label_atom(ex.carbon), it and its parent cannot have the same structure
  for(int i=this->nodes[index_child].nth-1;i>=0;i--){
    if(this->get_label(this->nodes[index_parent].children[i])==this->get_label(index_child)){
      return false;
    }
    }*/

  //std::array<int,max_valence> cp_child=this->nodes[index_child].bond_position;
  //std::array<int,max_valence> cp_parent=this->nodes[index_parent].bond_position;
  carbon_position cp_child,cp_parent;
  generate_adj_list(cp_child,index_child);
  generate_adj_list(cp_parent,index_parent);
  //delete information of child benzene from cp_parent
  for(int i=0;i<cp_parent.size();i++){
    for(int j=0;j<cp_parent[i].size();j++){
      if(cp_parent[i][j]==this->nodes[index_child].nth){
	if(j==0 && cp_parent[i].size()==this->nodes[index_child].multi){
	  cp_parent.erase(cp_parent.begin()+i);
	}else{
	  cp_parent[i].erase(cp_parent[i].begin()+j);
	  if(this->nodes[index_child].multi==2){
	    cp_parent[i].erase(cp_parent[i].begin()+j+1);
	  }
	}
	break;
      }
    }
  }
  //std::sort(cp_child.begin(),cp_child.end(),compare_function);
  //std::sort(cp_parent.begin(),cp_parent.end(),compare_function);
  get_bond_position(cp_child,index_child);
  get_bond_position(cp_parent,index_parent);  
  
  if(this->nodes[index_child].multi==1){
    int nth = this->nodes[index_child].nth;
    //compare carbon position of child nodes (excluding index1 in cp_parent)
    //parent must be root so can compare by cp_child==cp_parent
    /*for(int i=0;i<this->nodes[index_child].num_children;i++){
      int shift=0;
      if(i>=nth){
	shift=1;
      }
      if(cp_child[i]!=cp_parent[i+shift]){
	return false;
      }
      }*/
    if(cp_child.size()!=cp_parent.size()){
      std::cout<<"error!!"<<std::endl;
    }else{
      if(cp_child!=cp_parent){
	return false;
      }
    }
    return true;
  }else if(this->nodes[index1].multi==2){
    //used to shift carbon position of child (in parent bond_position) to 0
    //not consider the case that more than 2 benzenes merge together ================================
    //if consider that case, cp of naph bond must be included (now only one cp per one atom from cp_child/cp_parent)
    int cp_shift = this->nodes[index_parent].bond_position[this->nodes[index_child].nth];
    
    /*for(int i=0;i<this->nodes[index_parent].num_children;i++){
      cp_parent[i] = (cp_shift-cp_parent[i]+6)%6;
      }*/
    for(int i=0;i<cp_parent.size();i++){
      for(int j=0;j<cp_parent[i].size();j++){
	cp_parent[i][j] = (cp_shift-cp_parent[i][j]+6)%6;
      }
      std::sort(cp_parent[i].begin(),cp_parent[i].end());
    }
    
    bool leftright_sym = true;
    //left-right symmetry between 2 benzenes
    /*
    for(int i=0;i<this->nodes[index_child].num_children;i++){
      int shift=0;
      if(i>=this->nodes[index_child].nth){
	shift=1;
      }
      if(cp_child[i]!=cp_parent[i+shift]){
	leftright_sym = false;
      }
      }*/
    if(cp_child.size()!=cp_parent.size()){
      std::cout<<"error!!"<<std::endl;
    }else{
      if(cp_child!=cp_parent){
	leftright_sym = false;
      }
    }
    if(leftright_sym){
      return true;
    }

    //left-right and up-down symmetry between 2 benzenes
    /*for(int i=0;i<this->nodes[index_child].num_children;i++){
      int shift=0;
      if(i>=this->nodes[index_child].nth){
	shift=1;
      }
      if(cp_child[i]!=5-cp_parent[i+shift]){
	return false;
      }
      }*/
    for(int i=0;i<cp_parent.size();i++){
      for(int j=0;j<cp_parent[i].size();j++){
	cp_parent[i][j] = 5-cp_parent[i][j];
      }
      std::sort(cp_parent[i].begin(),cp_parent[i].end());
    }

    if(cp_child.size()!=cp_parent.size()){
      std::cout<<"error!!"<<std::endl;
    }else{
      if(cp_child!=cp_parent){
	return false;
      }
    }
    return true;
  }
  std::cout<<"error in chain_symmetry case:1,2"<<std::endl;
  return true; 
}

bool ChemTreeCenter::is_tsub_equal(int index1, int index2){
  if(get_label(index1) != get_label(index2)){
    return false;
  }
  
  std::vector<int> child_index1, child_index2;
  for(int i=0; i<this->nodes[index1].num_children; i++){
    child_index1.push_back(nodes[index1].children[i]);
  }
  for(int i=0; i<this->nodes[index2].num_children; i++){
    child_index2.push_back(nodes[index2].children[i]);
  }

  bool is_adjacent = false;

  if(nodes[index1].parent == index2){
    is_adjacent = true;
    child_index2.erase(child_index2.begin() + nodes[index1].nth);
  }else if(nodes[index2].parent == index1){
    is_adjacent = true;
    child_index1.erase(child_index1.begin() + nodes[index2].nth);
  }

  if(child_index1.size() != child_index2.size()){
    return false;
  }

  if(!is_adjacent){
    if(nodes[index1].multi != nodes[index2].multi){
      return false;
    }
  
    if(nodes[index1].parent != nodes[index2].parent){
      if(!adjacent(nodes[index1].parent, nodes[index2].parent)){
	//if parent of both nodes are not adjacent, can check nth directly 
	if(nodes[index1].nth != nodes[index2].nth){
	  return false;
	}
      }else{
	//if parent of both nodes are adjacent, must care about nth of one parent in nth 
	int parent1 = nodes[index1].parent;
	int parent2 = nodes[index2].parent;
	if(nodes[ parent1 ].parent == parent2){
	  int nth = nodes[parent1].nth;
	  int shift = (nth>nodes[index2].nth ? 0:-1);
	  if(nodes[index1].nth != nodes[index2].nth+shift){
	    return false;
	  }
	}else{
	  int nth = nodes[parent2].nth;
	  int shift = (nth>nodes[index1].nth ? 0:-1);
	  if(nodes[index1].nth+shift != nodes[index2].nth){
	    return false;
	  }
	}
      }
    }
  }else{ //index1 and index2 are adjacent
    if(index2!=0){
      return false;
    }
  }

  //check if all child nodes has the same structure or not
  for(int i=0; i<child_index1.size(); i++){
    if(!is_equal(child_index1[i], child_index2[i])){
      return false;
    }
  }
  return true;
}

bool ChemTreeCenter::is_tsub_cp_equal(int index1, int index2, int path_index1, int path_index2, const std::vector<int> &path){  
  std::vector<int> child_index1, child_index2;
  for(size_t i=0; i<nodes[index1].num_children; i++){
    child_index1.push_back(nodes[index1].children[i]);
  }
  for(size_t i=0; i<nodes[index2].num_children; i++){
    child_index2.push_back(nodes[index2].children[i]);
  }

  if(nodes[index1].parent == index2){
    child_index2.erase(child_index2.begin() + nodes[index1].nth);
  }else if(nodes[index2].parent == index1){
    child_index1.erase(child_index1.begin() + nodes[index2].nth);
  }

  if(get_label(index1) < num_special_atom && get_label(index2) < num_special_atom){
    if(!found(path, index1) && !found(path, index2) ){
      //check only equal or not -> no need to use carbon_position_list, just bond_position
      std::vector<int>  cp_index1, cp_index2;
      cp_index1 = child_index1;
      cp_index2 = child_index2;
      
      if(index1 != 0) cp_index1.insert(cp_index1.begin(), nodes[index1].parent);
      if(index2 != 0) cp_index2.insert(cp_index2.begin(), nodes[index2].parent);
      
      get_bond_position(cp_index1, index1);
      get_bond_position(cp_index2, index2);

      if(cp_index1 != cp_index2) return false;
    }
  }
  
  bool result = true;
  for(size_t i=0; i<child_index1.size(); i++){
    if(path_index1==0 || (child_index1[i] > path[path_index1-1] && child_index2[i] > path[path_index2+1]) ){
    result = result && is_tsub_cp_equal(child_index1[i], child_index2[i], path_index1-1, path_index2+1, path);
    }
  }
  return result;
}

bool ChemTreeCenter::chain_end_symmetry(int index1,int index2){
  
  std::vector<int> children1,children2;
  for(int i=0;i<this->nodes[index1].num_children;i++){
    children1.push_back(this->nodes[index1].children[i]);
  }
  
  for(int i=0;i<this->nodes[index2].num_children;i++){
    children2.push_back(this->nodes[index2].children[i]);
  }
  
  bool is_adjacent = false; 
  //is two end-nodes adjacent with each other or not

  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent==index2){
    is_adjacent = true;
    children2.erase(children2.begin()+this->nodes[index1].nth);
  }else if(this->nodes[index2].parent==index1){
    is_adjacent = true;
    children1.erase(children1.begin()+this->nodes[index2].nth);
  }

  if(!is_adjacent){
    if(this->nodes[index1].multi!=this->nodes[index2].multi){
      return false;
    }
    if(children1.size()!=children2.size()){
      return false;
    }
  }else{
    if(index2!=0){
      return false;
    }
    if(children1.size()!=children2.size()){
      return false;
    }
  }

  //check if all adjacent nodes have the same structure or not
  for(int i=0;i<children1.size();i++){
    if(!is_equal(children1[i],children2[i])){
      return false;
    }
  }
  return true;
}

bool ChemTreeCenter::is_tsub_has_substr(int index){
  bool result = false;
  /*if(get_label(index) < num_special_atom){
    return true;
    }*/
  for(int i=0; i<nodes[index].num_children; i++){
    if(nodes[nodes[index].children[i]].num_children == 0){
      return false;
    }
    result = result || is_tsub_has_substr(nodes[index].children[i]);
  }
  return result;
}

bool ChemTreeCenter::is_chain_end_redundant(std::vector<int> chain,const carbon_position &carbon_pos){
  //get index of benzene chain and check if both end is redundant or not
  int index1 = chain[0],index2 = chain[chain.size()-1];
  
  bool is_adjacent = false; 
  //is two end-nodes adjacent with each other or not

  //remove index1/index2 if it is another node's child (consider other adjacent nodes only)
  if(this->nodes[index1].parent==index2 || this->nodes[index2].parent==index1){
    is_adjacent = true;
  }	
     
  if(this->nodes[index1].multi==1 && this->nodes[index2].multi==1){
    carbon_position co_carbon_pos;
    generate_adj_list(co_carbon_pos,index2);
    this->get_bond_position(co_carbon_pos,index2);
    
    int shift = 0;
    size_t step = 0;
    //step uses to track the nth of current node corresponding to carbon_pos[i][j]
    //if index2 is not root it have parent so nth starts from -1(parent)
    if(index2 != 0){
      step = -1;
    }
    
    for(size_t i=0; i<carbon_pos.size(); i++){    
      for(size_t j=0; j<carbon_pos[i].size(); j++){
	if(carbon_pos[i][j]!=-1){
	  if(is_adjacent && step >= nodes[index1].nth){ 
	    shift = 1;
	  }    

	  if(carbon_pos[i][j] > co_carbon_pos[i][j]){
	    return true;
	  }else if(carbon_pos[i][j] < co_carbon_pos[i][j]){
	    return false;
	  }
	}else{
	  break;
	}
	step++;
      }
    }
    return false;
  }/*else if(this->nodes[index1].multi==2){
    
    carbon_position c_index1;
    generate_adj_list(c_index1,index1);
    std::sort(c_index1.begin(),c_index1.end(),compare_function);
    int j_parent =-1;
    int k_parent =-1;

    for(size_t i=0;i<c_index1.size();i++){
      std::vector<int>::iterator it= std::find(c_index1[i].begin(),c_index1[i].end(),nodes[index1].parent);
      if(it!=c_index1[i].end()){
	j_parent = i;
	k_parent = it-c_index1[i].begin();
	break;
      }
    }
    if(j_parent ==-1 || k_parent==-1){
      std::cout<<"error"<<std::endl;
    }
    int carbon_pos_parent = carbon_pos[j_parent][k_parent];
    
    //parent is not assigned carbon position yet
    if(carbon_pos_parent==-1){
      return false;
    }
    //convert carbon position of index1 to be in the range of 1-2-3-4 based on carbon_pos_parent    	  
    std::vector< std::vector<int> > new_carbon_pos1;
    new_carbon_pos1 = carbon_pos;

    //erase parent from new_carbon_pos1
    new_carbon_pos1.erase(new_carbon_pos1.begin()+j_parent);

    if(carbon_pos_parent!=0){ // if carbon_pos_parent !=0 no need to change value of new_carbon_pos1
      for(size_t i=0;i<new_carbon_pos1.size();i++){
	for(size_t j=0;j<new_carbon_pos1[i].size();j++){
	  if(new_carbon_pos1[i][j]!=-1){
	    new_carbon_pos1[i][j] = (new_carbon_pos1[i][j]+5-carbon_pos_parent)%6; 
	  }
	}
	std::sort(new_carbon_pos1[i].begin(),new_carbon_pos1[i].end());
      }
    }
    std::vector< std::vector<int> > new_carbon_pos2;
    new_carbon_pos2 = new_carbon_pos1;

    if(is_fused_benzene(chain)){
      //rotate new_carbon_pos2's carbon position in reverse direction (1-2-3-4) to (4-3-2-1)
      //necessary only when two benzene is in the same cyclic structure
      for(int i=0;i<new_carbon_pos2.size();i++){
	for(int j=0;j<new_carbon_pos2[i].size();j++){
	  if(new_carbon_pos2[i][j]!=-1){
	    new_carbon_pos2[i][j] = 5-new_carbon_pos2[i][j];
	  }
	}
	std::sort(new_carbon_pos2[i].begin(),new_carbon_pos2[i].end());
      }
    }
	 
    bool check_newcarbon = true;
    bool check_carbon = true;
    bool leftright_updown_sym = false;

    if(is_adjacent){
      leftright_updown_sym = true;
    }else if(is_fused_benzene(chain)){
      int root_index; //position of the top node (node with lowest depth) in the chain
      for(int i=0;i<chain.size();i++){
	if(chain[i]<chain[i+1]){
	  root_index = i;
	  break;
	}
      }
      carbon_position temp_c_index;
      generate_adj_list(temp_c_index,chain[root_index]);
      std::sort(temp_c_index.begin(),temp_c_index.end(),compare_function);
      //check if chain is a straight chain of benzene rings or not and can be flipped in both left-right and up-down direction
      if(abs(this->nodes[chain[root_index]].bond_position[this->nodes[chain[root_index-1]].nth]-this->nodes[chain[root_index]].bond_position[this->nodes[chain[root_index+1]].nth] )==3 && temp_c_index[0].size()%2==0  && chain[root_index]==0){
	//root_index has even child and no parent
	  
	bool no_odd_child = true;
	for(int i=root_index-1;i>0;i--){
	  temp_c_index.clear();
	  generate_adj_list(temp_c_index,chain[i]);
	  if(this->nodes[chain[i]].bond_position[this->nodes[chain[i-1]].nth]!=2){
	    no_odd_child = false;	     
	    break;
	  }
	  for(int j=0;j<temp_c_index.size();j++){
	    if(temp_c_index[j].size()%2!=0){
	      no_odd_child=false;
	      break;
	    }
	  }
	}
	if(no_odd_child){
	  leftright_updown_sym = true;
	}
	
	// check only first half because the first half and the latter half are the same (from chain_symmetry)
	//if(no_odd_child){
	  //for(int i=root_index+1;i<chain.size()-1;i++){
	  //  temp_c_index.clear();
	  //  generate_adj_list(temp_c_index,chain[i]);
	  //  if(temp_c_index[0].size()%2!=0 || this->nodes[chain[i]].bond_position[this->nodes[chain[i+1]].nth]!=2){
	  //    no_odd_child = false;
	  //    break;
	  //}
	  //}
	  //if(no_odd_child){
	  //  leftright_updown_sym = true;
	  //}
	  //}
      }
    }
    carbon_position index2_c_index; //store nth of benzene's child
    generate_adj_list(index2_c_index,index2);
    std::sort(index2_c_index.begin(),index2_c_index.end(),compare_function);
    
    if(is_adjacent){
      //in this case index2 is parent of index1

      int nth = this->nodes[index1].nth;
      int cp_shift = 6+this->nodes[index2].bond_position[nth+1];
      //delete information of index1 (child benzene) from  index2_c_index
      //delete information of parent of index2 from index2_c_index
      //and convert carbon position of nodes adjacent to index2 to be between 1-4 based on cp_shift
      for(int i=0;i<index2_c_index.size();i++){
	for(int j=0;j<index2_c_index[i].size();j++){
	  if(index2_c_index[i][j]==index1){
	    if(j==0 && index2_c_index[i].size()==2){
	      index2_c_index.erase(index2_c_index.begin()+i);
	      j=-1;//make j at new iteration =0 while i is still the same
	    }else{
	      index2_c_index[i].erase(index2_c_index[i].begin()+j);
	      index2_c_index[i].erase(index2_c_index[i].begin()+j+1);
	      j--;
	    }
	  }else if(index2_c_index[i][j]==this->nodes[index2].parent){
	    index2_c_index.erase(index2_c_index.begin()+i);
	    j=-1;//make j at new iteration =0 while i is still the same
	  }else{
	    if(carbon_pos[i][j]!=-1){
		int nth2 = nodes[index2_c_index[i][j]].nth;
		if(this->nodes[index2].bond_position[nth+1]!=0 && this->nodes[index2].bond_position[nth+1]<3){
		  index2_c_index[i][j] = (cp_shift-this->nodes[index2].bond_position[nth2+1])%6;
		}else{
		  index2_c_index[i][j] = (this->nodes[index2].bond_position[nth2+1]+5-this->nodes[index2].bond_position[nth+1])%6;
		}
	    }else{
	      if(j!=0){
		return false;
	      }
	    }
	  }
	}
	std::sort(index2_c_index[i].begin(),index2_c_index[i].end());
      }
    }else{ //not adjacent
     
      //not adjacent -> remove parent from index2_c_index and adjust it to 1-2-3-4 based on carbon position of parent
      int j_parent=-1;
      int k_parent=-1;
      for(size_t i=0;i<index2_c_index.size();i++){
	std::vector<int>::iterator it= std::find(index2_c_index[i].begin(),index2_c_index[i].end(),nodes[index2].parent);
	if(it!=index2_c_index[i].end()){
	  j_parent = i;
	  k_parent = it-index2_c_index[i].begin();
	  break;
	}
      }
      
      if(j_parent==-1 || k_parent==-1){
	std::cout<<"error line 1524"<<std::endl;
      }

      int carbon_pos_parent = nodes[index2].bond_position[0];

      //erase parent from new_carbon_pos1
      index2_c_index.erase(index2_c_index.begin()+j_parent);

      //convert carbon position of index1 to be in the range of 1-2-3-4 based on carbon_pos_parent    	  
      if(carbon_pos_parent!=0){ // if carbon_pos_parent !=0 no need to change value of new_carbon_pos1
	for(size_t i=0;i<index2_c_index.size();i++){
	  for(size_t j=0;j<index2_c_index[i].size();j++){
	    if(index2_c_index[i][j]!=-1){
	      index2_c_index[i][j] = (nodes[index2].bond_position[nodes[index2_c_index[i][j]].nth+1]+5-carbon_pos_parent)%6; 
	    }else{
	      index2_c_index[i][j] = (nodes[index2].bond_position[nodes[index2_c_index[i][j]].nth+1]+1+5-carbon_pos_parent)%6;
	    }
	  }
	  std::sort(index2_c_index[i].begin(),index2_c_index[i].end());
	}
      } 
    }
    
    for(int i=0;i<new_carbon_pos1.size();i++){    
      for(int j=0;j<new_carbon_pos1[i].size();j++){
	if(new_carbon_pos1[i][j]!=-1){
	   if(leftright_updown_sym){  
	    //their common parent does not have parent, has only even unique childs
	    //if(debug()) std::cout<<"left right sym   "<<carbon_pos[i][j]<<"   "<<new_carbon_pos[i][j]<<"    "<<index2_c_index[i][j]<<std::endl;
	    if(check_newcarbon && check_carbon){
	      if(new_carbon_pos2[i][j]>index2_c_index[i][j] && new_carbon_pos1[i][j]>index2_c_index[i][j]){
		return false;
	      }
	      if(new_carbon_pos2[i][j]<index2_c_index[i][j] || new_carbon_pos1[i][j]<index2_c_index[i][j]){
		return true;
	      }
		
	      //when reach here either new_carbon_pos or carbon_pos (but not both) must equal to index2.bond_position[....]
	      if(new_carbon_pos2[i][j]==index2_c_index[i][j]){
		//if new_carbon_pos is equal, we will compare only new carbon_pos in the next round 
		check_carbon = false;
	      }else{
		check_newcarbon = false;
	      }
	    }else if(check_carbon && !check_newcarbon){
	      if(new_carbon_pos1[i][j]>index2_c_index[i][j]){
		return false;
	      }
	      if(new_carbon_pos1[i][j]<index2_c_index[i][j]){
		return true;
	      } 
	    }else if(!check_carbon && check_newcarbon){
	      if(new_carbon_pos2[i][j]>index2_c_index[i][j]){
		return false;
	      }
	      if(new_carbon_pos2[i][j]<index2_c_index[i][j]){
		return true;
	      }
	    }
		  
	  }else{
	    if(debug() && code()) std::cout<<"not sym  cp "<<new_carbon_pos2[i][j]<<" cp2 "<<nodes[index2].bond_position[c_index[i][j]]<<std::endl;
	    if(new_carbon_pos2[i][j]>index2_c_index[i][j]){
	      return false;
	    }
	    if(new_carbon_pos2[i][j]<index2_c_index[i][j]){
	      return true;
	    }
	  }
    
	}//carbon_pos[i][j]!=-2	  
      }//for j of carbon_pos[i][j]
    }//for i of carbon_pos[i][j]
    return false;
  }//multi==2  */

  return false;
}

/*bool ChemTreeCenter::chain_middle_symmetry(std::vector<int> chain){
	int middle = chain.size()/2;
	
	if(chain[middle]>chain[middle+1] || chain[middle]>chain[middle-1]){
	  //if middle of chain is not parent of its adjacent node in chain then it cannot be symmetry
	  return false;
	}
	if(this->nodes[chain[middle+1]].multi!=this->nodes[chain[middle-1]].multi){
	  return false;
	}
	if(this->get_label(chain[middle])!=0){
	  return true;
	}

	std::vector< std::vector<int> > c_index,carbon_pos,flip_carbon_pos; //order of children having the same structure are grouped together 
	generate_adj_list(c_index,chain[middle]);
	//std::sort(c_index.begin(),c_index.end(),compare_function);
	
	if(chain[middle]==0 && c_index.size()==1 && c_index[0].size()==2 && chain[middle]==0){
	  //if its adjacent node has only two benzene with single bond and not have parent, it is always symmetry
	  return true;
	}

	carbon_pos = c_index;
	flip_carbon_pos = c_index;
	//initialize carbon_position of corresponding child node
	for(int i=0;i<c_index.size();i++){
	  for(int j=0;j<c_index[i].size();j++){
	    if(c_index[i][j]==nodes[chain[middle]].parent){
	      carbon_pos[i][j] = this->nodes[chain[middle]].bond_position[0];
	    }else if(c_index[i][j]!=-1){
	      carbon_pos[i][j] = this->nodes[chain[middle]].bond_position[c_index[i][j]];
	    }else{
	      carbon_pos[i][j] = carbon_pos[i][j-1]+1;
	    }
	  }
	}

	if(c_index[0].size()==1 && c_index[0][0]!=this->nodes[chain[middle+1]].nth && c_index[0][0]!=this->nodes[chain[middle-1]].nth){
	  for(int i=0;i<c_index.size();i++){
	    for(int j=0;j<c_index[i].size();j++){
	      flip_carbon_pos[i][j] = (6-carbon_pos[i][j])%6;
	    }
	    std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	  }
	  if(flip_carbon_pos == carbon_pos){
	    return true;
	  }else{ 
	    return false; 
	  }    
	}

	float c_pos1=this->nodes[chain[middle]].bond_position[this->nodes[chain[middle-1]].nth];
	float c_pos2=this->nodes[chain[middle]].bond_position[this->nodes[chain[middle+1]].nth];

	//check symmetry of benzene based on position of two adjacent benzene
	switch(abs(c_pos1-c_pos2)){
	case 1:
	case 2:
	case 4:
	case 5:
	  {  
	    int sum_pos = static_cast<int>(c_pos1+c_pos2)%6;
	    if(this->nodes[chain[middle+1]].multi==2 && this->nodes[chain[middle-1]].multi==2){
	      sum_pos+=1;
	    }
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		flip_carbon_pos[i][j] = (sum_pos+6-carbon_pos[i][j])%6;
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }
	    
	    if(flip_carbon_pos == carbon_pos)
	      return true;
	    else 
	      return false;     
	    break;
	  }
	case 3:
	  {
	    //flip by left-right reflection
	    float min = std::min(c_pos1,c_pos2);
	    if(this->nodes[chain[middle+1]].multi==2 && this->nodes[chain[middle-1]].multi==2){
	      min+=0.5;
	    }
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		flip_carbon_pos[i][j] = static_cast<int>(2*min+6-carbon_pos[i][j])%6;
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }
	    if(flip_carbon_pos == carbon_pos){
	      return true;
	    }
	    //flip by up-down reflection
	    int sum_pos = static_cast<int>(c_pos1+c_pos2)%6;
	    if(this->nodes[chain[middle+1]].multi==2 && this->nodes[chain[middle-1]].multi==2){
	      sum_pos+=1;
	    }
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		if(carbon_pos[i][j]>sum_pos){
		  flip_carbon_pos[i][j] = sum_pos+6-carbon_pos[i][j];
		}else{
		  flip_carbon_pos[i][j] = sum_pos-carbon_pos[i][j];
		}
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }
	    if(flip_carbon_pos == carbon_pos){
	      return true;
	    }
	    //flip by both left-right and up-down reflection
	    for(int i=0;i<c_index.size();i++){
	      for(int j=0;j<c_index[i].size();j++){
		flip_carbon_pos[i][j] = (carbon_pos[i][j]+3)%6;
	      }
	      std::sort(flip_carbon_pos[i].begin(),flip_carbon_pos[i].end());
	    }
	    if(flip_carbon_pos == carbon_pos){
	      return true;
	    }else{
	      return false;
	    }
	    break;
	  }  
	}
	std::cout<<"error in middle_chain_symmetry"<<std::endl;
	return false;
}*/

int ChemTreeCenter::chain_middle_pair_symmetry(const std::vector<int> &path, const carbon_position &cp_list, const int latest_unequal){
  int index1 = path[path.size()/2-1];
  int index2 = path[path.size()/2];

  if(index1 > index2){
    int temp = index1;
    index1 = index2;
    index2 = temp;
  }

  bool is_adjacent = false;
  int index_parent, index_child;

  if(nodes[index1].parent == index2){
    index_parent = index2;
    index_child = index1;
    is_adjacent = true;
  }else if(nodes[index2].parent == index1){
    index_parent = index1;
    index_child = index2;
    is_adjacent = true;
  }
  
  if(!is_adjacent){
    carbon_position cp_index1, cp_index2;

    cp_index1 = cp_list;
    generate_adj_list(cp_index2, index2);
    //std::sort(cp_index2.begin(), cp_index2.end(), compare_function);
    get_bond_position(cp_index2, index2);

    if(cp_index1 < cp_index2)   return  1;
    if(cp_index1 > cp_index2)   return -1;
    return 0;   //if(cp_index1 == cp_index2)
  }else{//is_adjacent
    if(index_parent != 0) return false;

    carbon_position cp_parent, adj_list_parent;
    cp_parent = cp_list;
    generate_adj_list(adj_list_parent, index_parent);

    carbon_position cp_child;
    generate_adj_list(cp_child, index_child);
    get_bond_position(cp_child, index_child);

    int nth = nodes[index_child].nth;

    if(cp_child.size() != cp_parent.size()){
      std::cout<<"Error"<<std::endl;
      return 0;
    }else{
      //because index_child > index_parent always holds 
      if(cp_child < cp_parent)	return -1; 
      if(cp_child > cp_parent)  return  1; 
      return 0; //if(cp_child == cp_parent) 
    }
  }// end if-else is_adjacent
}

void ChemTreeCenter::find_sympath_down(int step,int index_step,std::vector<int> &chain,std::vector< std::vector<int> > &sympath_collection){
  //index_step is index of step in temp ex.if step is the first element, index_step = 0
  std::array<int, max_valence> child = this->nodes[step].children;
  
  for(int i=0;i<this->nodes[step].num_children;i++){
    int index_child = this->nodes[step].children[i];
    //if(this->nodes[index_child].num_children > 0 && 
    if(!found(chain,index_child) && index_child<chain[0]){ 
      //index_child < step to assure that one chain is checked only once
      chain.push_back(index_child);
      
      find_sympath_down(index_child,index_step-1,chain,sympath_collection);
      chain.pop_back();
    }
  }
  if(this->get_label(step) < num_special_atom){
    sympath_collection.push_back(chain);
  }		
}

void ChemTreeCenter::find_sympath(std::vector< std::vector<int> > &sympath_collection){
  for(size_t index=0; index<this->get_num_nodes(); index++){
    if(get_label(index) < num_special_atom){
      //find list of index of benzene chain 
      std::vector<int> upward_chain;
      upward_chain.push_back(index);
      int step = this->nodes[index].parent;

      while(step>=0){
	upward_chain.push_back(step);

	this->find_sympath_down(step,upward_chain.size()-1,upward_chain, sympath_collection);
	step = this->nodes[step].parent;
      }
    }
  }
}

bool ChemTreeCenter::is_normal_sympath(const std::vector< std::vector<int> > & cp_list, 
const std::vector< std::vector<int> > & sympath_collection, 
const int index, 
const std::vector< auto_group > &aut_group){
	/*if(k!=position[j].size()-1){
	  //check for normal chain if all carbon in set j is known (k is the last element of j) only.
	  return true;
	  }*/

	for(size_t i=0; i<sympath_collection.size(); i++){		
	  //remove normal node from sympath to make it simpler to check (only substructure nodes remain)
	  std::vector<int> str_path;
	  std::vector<int> sympath = sympath_collection[i];
	  for(size_t j=0; j<sympath.size(); j++){
	    if(get_label(sympath[j]) < num_special_atom){
	      str_path.push_back(sympath[j]);
	    }
	  }   
	  
	  //if it is the most middle special node in the chain -> check for normal chain; otherwise, do nothing
	  if(*std::min_element(str_path.begin(), str_path.end()) == index){
	    
	    /*
	    if(str_path.size()%2 != 0){
	      bool is_sym = false;
	      
	      //if the middle node is structure node, check if its structure (it and its adjacent nodes) is symmetric or not	     
	      //but consider subtrees of two child nodes in sympath as if they have same structures 
	      //although atom position of some of their nodes are different
	      //so create a copy of cp_list that always treats subtrees of two child nodes as the same structure (let them be in the same set) store it in newcp_list
	      
	      std::vector< std::vector<int> > adj_list;
	      generate_adj_list(adj_list,index);
	      //print_cp(adj_list);

	      int index_in_sympath = std::min_element(sympath.begin(),sympath.end()) - sympath.begin(); //get index of current node in sympath
	      //std::cout<<"index two child nodes:" << sympath[index_in_sympath-1]<<" and "<< sympath[index_in_sympath +1] <<std::endl;

	      std::vector< std::vector<int> > newcp_list = cp_list;
	      for(int adj_i=0; adj_i < adj_list.size(); adj_i++){
		//std::cout<<"adj_i:"<<adj_i<<std::endl;
		//print_cp(adj_list);
		std::vector<int> adj_set = adj_list[adj_i];

		if( !found(adj_set, sympath[index_in_sympath-1] ) && !found(adj_set, sympath[index_in_sympath+1]) ){ 
		  //if those two child nodes are not in this set of nodes -> do nothing
		  //std::cout<<"both are not in"<<adj_i<<std::endl;
		}else{

		  if( found(adj_set, sympath[index_in_sympath-1]) && found(adj_set, sympath[index_in_sympath+1])){
		    //if both child nodes in this set -> do nothing
		    //std::cout<<"both are in"<<adj_i<<std::endl;
		  }else{		    
		    //otherwise, merge another child node in the set before add it in newcp_list

		    for(int adj_j=adj_i+1; adj_j < adj_list.size(); adj_j++){
		      //std::cout<<"adj_j:"<<adj_j<<std::endl;
		      
		      int index_in_adj = -1;

		      if( found( adj_list[adj_j], sympath[index_in_sympath-1]) ){ 
			index_in_adj = std::find(adj_list[adj_j].begin(), adj_list[adj_j].end(), sympath[index_in_sympath-1] ) - adj_list[adj_j].begin();
		      }else if( found(adj_list[adj_j], sympath[index_in_sympath+1])){
			index_in_adj = std::find(adj_list[adj_j].begin(), adj_list[adj_j].end(), sympath[index_in_sympath+1] ) - adj_list[adj_j].begin();
		      }else{
			

			print_cp(adj_list);
			std::cout<<"Error in normal chain"<<std::endl;
			exit(0);
			return false;
		      }
		      
		      //temp.push_back(cp_list[j][index_in_adj]);
		      newcp_list[adj_i].push_back(newcp_list[adj_j][index_in_adj]);
		      adj_list[adj_i].push_back(adj_list[adj_j][index_in_adj]);
	
		      //if set containing another child (of those two chil nodes) has one member remove that set; otherwise remove onlythat member
		      if(adj_list[adj_j].size() == 1){
			newcp_list.erase( newcp_list.begin() + adj_j);
			adj_list.erase( adj_list.begin() + adj_j);
		      }else{
			newcp_list[adj_j].erase( newcp_list[adj_j].begin() + index_in_adj);
			adj_list[adj_j].erase( adj_list[adj_j].begin() + index_in_adj);
		      }
		      if(index_in_adj != -1){
			break;
		      }
					
		    }
		    break;
		  }
		}
		}

	      for(int aut_i=1; aut_i < aut_group[get_label(index)].size(); aut_i++){
		std::vector< std::vector<int> > flip_cp_list = newcp_list;
		get_automorphism(flip_cp_list, aut_group[get_label(index)][aut_i]);
		if(flip_cp_list == newcp_list){
		  is_sym = true;
		}
	      }
	      if(!is_sym){
		//if it is not symmetric -> no need to check for symmetric path
		return true;
	      }
	    }*/

	    if(sympath_collection[i].size() > 3){
	      size_t middle_left = sympath_collection[i].size()/2-1;
	      size_t middle_right = (sympath_collection[i].size()%2 == 0? sympath_collection[i].size()/2: sympath_collection[i].size()/2+1);

	      //check if carbon position of nodes in tsub of two middle nodes is same or not
	      //if it is not the same, always normal
	      if(!is_tsub_cp_equal(sympath_collection[i][middle_left], sympath_collection[i][middle_right], middle_left, middle_right, sympath_collection[i]) ){
		return true;
	      }
	    }

	    int latest_unequal = 0;
	    //keep the latest unequal carbon position of a pair of special atom
	    
	    int max = (str_path.size()%2 == 0? str_path.size()/2-1: str_path.size()/2 );
	    //check chain_symmetry for  all pairs of node except the middle one or the most middle pair
	    for(size_t path_i=0; path_i < max; path_i++){
	      int compare_result = chain_symmetry(str_path[path_i], str_path[str_path.size()-1-path_i]);
	      //-1 if cp[index1] > cp[index2], where index1 < index2
	      // 1 if cp[index1] < cp[index2], where index1 < index2
	      // 0 if cp[index1] = cp[index2], where index1 < index2
	      if(compare_result != 0){
		latest_unequal = compare_result;
	      }
	    }
	    if(str_path.size()%2==0){
	      int compare_result = chain_middle_pair_symmetry(str_path, cp_list, latest_unequal);

	      if(compare_result != 0){
		latest_unequal = compare_result;
	      }
	    }

	    if(latest_unequal == -1){
	      return false;
	    }//otherwise go to check the next sympath
	  }
	  /*if(is_chain_end_redundant(sympath_collection[i], cp_list)){
	    return false;
	    }*/
	}
	//if no symmetry occur for every chain -> it is normal_benzene
	//if every symmetry chain is not redundant -> it is normal_benzene
	return true;
}

bool ChemTreeCenter::is_normal_str(const std::vector< std::vector<int> > &cp_list, const std::vector< std::vector<int> > &aut_group){

  for(size_t aut_i=1; aut_i < aut_group.size(); aut_i++){
    std::vector< std::vector<int> > flip_cp_list = cp_list;

    get_automorphism(flip_cp_list, aut_group[aut_i]);
    if(flip_cp_list < cp_list){
      return false;
    }
    
    /*for(int a=0;a<j;a++){
      for(int b=0;b<position[a].size();b++){
      flip_position.push_back(flip_carbon(position[a][b],label,position[0],flip_mode));
      }
      std::sort(flip_position.begin(),flip_position.end());
      for(int b=0;b<flip_position.size();b++){
      if(flip_position[b]<position[a][b]){
      return false;
      }
      if(flip_position[b]>position[a][b]){
      if(k==position[j].size()-1 && a>0){
      normal = true;
      }
      return true;
      
      }
      }
      flip_position.clear();
      }*/

    //check normal form of child group j
    /*for(int b=0;b<k;b++){
      flip_position.push_back(flip_carbon(position[j][b],label,position[0],flip_mode));
    }
    flip_position.push_back(flip_carbon(carbon_pos,label,position[0],flip_mode)); //add new carbon position to flip vector
    std::sort(flip_position.begin(),flip_position.end());
    for(int b=0;b<k;b++){
      if(flip_position[b]<position[j][b]){
	return false;
      }
      if(flip_position[b]>position[j][b]){
	if(k==position[j].size()-1){
	  if(j>0 || get_label(index)!=0){
	    normal = true;
	  }
	}
	return true;
      }
    }
 
    //compare last carbon position
    if(flip_position[k]<carbon_pos){
      return false;
    }else if(flip_position[k]>carbon_pos){
      if(k==position[j].size()-1){
	if(j>0 || get_label(index)!=0){
	  normal = true;
	}
      }
      return true;
      } */
  }//end for loop aut_group

  return true;
}

size_t ChemTreeCenter::label_child_and_next(std::vector<carbon_position> &result,  int j,  int k,  int index,  const std::vector< std::vector<int> > & sympath_collection,  const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi, std::ofstream & outputfile){
	//j and k are indices of position that will be labeled by carbon position
	//index is index of benzene (used to determine if it is root or not)
 
        carbon_position current_cp = result[result.size()-1];
	if(j+1==current_cp.size() && k+1==current_cp[j].size()){	  

	  assign_bond_position(index, current_cp);
	  if( this->is_normal_sympath(current_cp, sympath_collection, index, aut_group) ){
	    return this->label_substructure(index-1, sympath_collection, result, aut_group, str_lack_valence, str_smi, outputfile);
	  }else{
	    return 0;
	  }
	}

	size_t num=0;
	int max_pos = aut_group[nodes[index].label][0].size();
	//int max_pos = str_lack_valence[get_label(index)].sum();  //max_pos = maximum carbon position

	//bool init_normal = normal;
	for(int cp=0; cp<max_pos; cp++){
	  //normal = init_normal;
	  int next_j = j;
	  int next_k = k+1;

	  if(k+1==current_cp[j].size()){
	    next_j = j+1;
	    next_k = 0;
	  }

	  if( max(current_cp[next_j] ) <= cp){
	    current_cp[next_j][next_k] = cp;
	    
	    if(!is_cp_available(current_cp, str_lack_valence[get_label(index)])){
	      //check if the node "cp" still has free valence or not
	      //need because one node "cp" can have more than 1 single bond
	      current_cp[next_j][next_k] = -1;
	      continue;
	    }

	    if(is_normal_str(current_cp, aut_group[get_label(index)])){ 
	      //if(k != current_cp[nj].size()-1 || this->is_normal_chain(current_cp, sympath_collection))
	      result[result.size()-1][next_j][next_k] = cp;
	      num+=label_child_and_next(result, next_j, next_k, index, sympath_collection, aut_group, str_lack_valence, str_smi, outputfile);
	      result[result.size()-1][next_j][next_k] = -1;
	    }	  

	    current_cp[next_j][next_k] = -1;
	  }
	}
	
	return num;	
}

bool compare_function(std::vector<int> i,std::vector<int> j){
  return i.size()<j.size();
}

bool path_compare_function(std::vector<int> path1, std::vector<int> path2){
  if(path1[0] != path2[0]){
    return path1[0] > path2[0];
  }else{
    return path1.size() < path2.size();
  }
}

size_t ChemTreeCenter::label_substructure(const int index, const std::vector< std::vector<int> > &sympath_collection, std::vector<carbon_position> &result,  const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi, std::ofstream &outputfile){
	using namespace std;

	if(aut_group.size()==0 || index == -1){
	  
	  if(do_print){
	    
	    //if(debug()){
	    std::cout<<"--------------------"<<std::endl;
	    this->show();

	    this->write(outputfile, str_smi);
	    outputfile << "\n";
	    return 1;
	  }
	  //return 0;
	  return 1;
	}
	
	size_t num = 0;

	if(this->get_label(index) < num_special_atom){

	  vector< vector<int> > adj_list; 
	  //group child and parent of node[index] with same structure together
	  //value in c_index is the index of that node ex. [ [0] [3 4] ] means node[0], node[3] and node[4] adjacent with this node and tsub(node[3]) = tsub(node[4])
	  generate_adj_list(adj_list,index);
	 
	  //sort c_index according to size and valence value
	  //std::sort(c_index.begin(),c_index.end(),compare_function);		
	  if(adj_list.size() > 0){
	    
	    carbon_position cp_list = adj_list;
	    for(size_t i=0; i < cp_list.size(); i++){
	      for(size_t j=0; j < cp_list[i].size(); j++){
		cp_list[i][j] = -1; // -1 is initialization
	      }
	    }

	    result.push_back(cp_list);
	  
	    int max_pos = aut_group[nodes[index].label][0].size();
	    //max_pos = maximum carbon position
	    //bool normal = false;

	    for(int cp=0; cp < max_pos; cp++){ 
	      //normal = false;
	      cp_list[0][0] = cp;
	      if(is_cp_available(cp_list, str_lack_valence[get_label(index)])){
		if( is_normal_str(cp_list, aut_group[get_label(index)]) ){
		  //if( 0 != cp_list[j].size()-1 || is_normal_chain(cp_list, sympath_collection) )
		  result[result.size()-1][0][0] = cp;
		  num+= label_child_and_next(result, 0, 0, index, sympath_collection, aut_group, str_lack_valence, str_smi, outputfile);
		}
	      }

	      cp_list[0][0] = -1; 
	    }
		
	    result.pop_back(); //erase(result.begin()+index);
	  }else{
	    //if it has no adjacent nodes (a tree contains only one node)
	    num = this->label_substructure(index-1, sympath_collection, result, aut_group, str_lack_valence, str_smi, outputfile);
	  }
	}else{ // else of if(index==0 && num_children>0)
	  carbon_position temp;
	  result.push_back(temp);
	  
	  num = this->label_substructure(index-1, sympath_collection, result, aut_group, str_lack_valence, str_smi, outputfile);
	  
	  result.pop_back(); //erase(result.begin()+index);
	}

	return num;
}

size_t enum_substructure(ChemTreeCenter& tree,std::ofstream& outputfile,const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi){
	size_t num=0; 
	std::vector<int> multi;
	
	if(num_lack_H == 0){ 
	//store previous information of multi in multi
		for(int index=0;index<tree.get_num_nodes();index++){
			multi.push_back(tree.get_multi(index));
			tree.set_multi_bond(index,1);
		}
	}
	std::vector<carbon_position> result;

	std::vector< std::vector<int> > sympath_collection;
	tree.find_sympath(sympath_collection);

	for(size_t i=0; i<sympath_collection.size(); i++){
	  bool skip = false;
	  for(size_t a=0; a<sympath_collection[i].size()/2; a++){
	    //check for subtree of each pair
	    int right_pair = sympath_collection[i].size()-1-a;
	    if(!tree.is_tsub_equal(sympath_collection[i][a], sympath_collection[i][right_pair]) ){
	      sympath_collection.erase(sympath_collection.begin()+i);
	      i--;
	      skip = true;
	      break;
	    }
	  }
	  if(!skip && tree.is_tsub_has_substr(sympath_collection[i][0]) ){
	    sympath_collection.erase(sympath_collection.begin()+i);
	    i--;
	    continue;
	  }
	  /*if(!skip && !chain_end_symmetry(sympath_collection[i][0], sympath_collection[i][sympath_collection[i].size()-1])){
	    sympath_collection.erase(sympath_collection.begin()+i);
	    i--;
	    continue;
	    } */ 
	}
	
	std::sort(sympath_collection.begin(), sympath_collection.end(), path_compare_function);

	if(tree.debug())
	  num = tree.label_substructure(tree.get_num_nodes()-1, sympath_collection, result, aut_group, str_lack_valence, str_smi, outputfile);

	result.clear();
	
	//outputfile<<"  num="<<num<<"\n";
	//tree.write(outputfile);
	//outputfile<<"=============\n";
	//if(sympath_collection.size()==0){
	//std::cout<<"====="<<num<<"======"<<std::endl;
	//tree.show();
	//std::cout<<"===================="<<std::endl;
	
	for(size_t index=0; index<tree.get_num_nodes(); index++){
	  tree.reset_bond_position(index);
	}

	if(num_lack_H == 0){
	//restore previous information of multi
		for(int index=0;index<multi.size();index++){
			tree.set_multi_bond(index,multi[index]);
		}
	}
	return num;
}

inline size_t set_multi_and_next(ChemTreeCenter& tree, is_ident_type is_ident, const int i, valence_value_type multiple, const int double_bond, const int triple_bond, const int normal_type,std::ofstream& outputfile,  const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi) 
{
	using namespace std;

	tree.set_multi_bond(i, multiple);
	if(double_bond==0 && triple_bond==0){
		tree.fill_rest_single_bond(i+1);
		if(normal_type >= 2){
			if(not tree.is_multi_normal(normal_type - 2)) {
				return 0;
			}
		}
		return enum_substructure(tree, outputfile, aut_group, str_lack_valence, str_smi);
	}

	--multiple;
	lack_degree[i] -= multiple;
	const int pi = tree.get_parent(i);
	lack_degree[pi] -= multiple;

	tree.update_identical_multi(is_ident, i);

	const int num_nodes = tree.get_num_nodes();

	size_t num = 0;
	for (int k = i + 1; k + triple_bond + double_bond  <= num_nodes; ++k) {
		if(tree.can_be_added_multi(k)){ 

			const valence_value_type c_degree = lack_degree[k];
			const int pk = tree.get_parent(k);
			const valence_value_type p_degree = lack_degree[pk];
			const valence_value_type maxmulti = tree.max_multi(is_ident, k);

			if(tree.get_label(k)!=0){

				if ((c_degree >= 1) and (p_degree >= 1) and (maxmulti >= 2)) {
					if (double_bond > 0) {
					  num += set_multi_and_next(tree, is_ident, k, 2, double_bond - 1, triple_bond, normal_type, outputfile, aut_group, str_lack_valence, str_smi);
					}
					if ((triple_bond > 0) and (c_degree >= 2) and (p_degree >= 2) and (maxmulti >= 3)) {
					  num += set_multi_and_next(tree, is_ident, k, 3, double_bond, triple_bond - 1, normal_type, outputfile, aut_group, str_lack_valence, str_smi);
					}
				}
			}
		}//end if(tree.can_be_added_multi(k))

		tree.set_multi_bond(k, 1);
		tree.update_identical_multi(is_ident, k);
	}

	lack_degree[i] += multiple;
	lack_degree[pi] += multiple;

	return num;
}

inline size_t generate_multi(ChemTreeCenter tree, // copy call because replaceing multi
			     is_ident_type& is_ident, const int normal_type,std::ofstream& outputfile,const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi)
{
	using namespace std;
	
	for (int i = tree.get_parent(num_except_H - 1); i != num_except_H - 1; ++i) {
		tree.update_identical_end(is_ident, i);
	}
	tree.calc_lack_degree();

	const int num_nodes = tree.get_num_nodes();
	
	size_t num = 0;
	for(int triple_bond = 0; triple_bond <= (num_lack_H >> 2); ++triple_bond){
		const int double_bond = (num_lack_H >> 1) - (triple_bond << 1);

		for (int k = 1; k + triple_bond + double_bond <= num_nodes; ++k) {
			if(tree.can_be_added_multi(k)){

				const valence_value_type c_degree = lack_degree[k];
				const int pk = tree.get_parent(k);
				const valence_value_type p_degree = lack_degree[pk];
				const valence_value_type maxmulti = tree.max_multi(is_ident, k);

				if(tree.get_label(k)>=num_special_atom){

					if ((c_degree >= 1) and (p_degree >= 1) and (maxmulti >= 2)) {
						if (double_bond > 0) {
						  num += set_multi_and_next(tree, is_ident, k, 2, double_bond - 1, triple_bond, normal_type, outputfile, aut_group, str_lack_valence, str_smi);
						}
						if ((triple_bond > 0) and (c_degree >= 2) and (p_degree >= 2) and (maxmulti >= 3)) {
						  num += set_multi_and_next(tree, is_ident, k, 3, double_bond, triple_bond - 1, normal_type, outputfile, aut_group, str_lack_valence, str_smi);
						}
					}

				}
			}//end if(tree.can_be_added_multi(k))

			tree.set_multi_bond(k, 1);  
			tree.update_identical_multi(is_ident, k);
		}
	}

	return num;
}

inline size_t add_node_and_next(ChemTreeCenter& tree, is_ident_type is_ident, const int deepest_head, const int parenti, const label_value_type atom_label,int &adjacent_benzene,std::ofstream& outputfile,const std::vector< auto_group > &aut_group,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi){
        using namespace std;
	
	if (tree.add_node(is_ident, parenti, atom_label)) {
		const int type = tree.is_normal(deepest_head);
		if (type > 0 ){
		  if (num_lack_H > 0){ 
		    return generate_multi(tree, is_ident, type, outputfile, aut_group, str_lack_valence, str_smi);
		  } else {
		    return enum_substructure(tree, outputfile, aut_group, str_lack_valence, str_smi);
		  }
		} else {
		  return 0;
		}
	}
	size_t num = 0;
	for (int i = tree.starti(); i < deepest_head; ++i){
		if (tree.can_be_added(i)) {
			const label_value_type begin_atom = tree.begin_atom_label(is_ident, i);
			for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
				if (tree.remain(atomi)) {				  
				  if(atomi==0 && tree.get_label(i)==0){
				    adjacent_benzene++;
				  }
				  num += add_node_and_next(tree, is_ident, deepest_head, i, atomi,adjacent_benzene,outputfile,aut_group,str_lack_valence, str_smi);
				  tree.del_last_node();
				  if(atomi==0 && tree.get_label(i)==0){
				    adjacent_benzene--;
				  }
				}
			}
		}
		tree.update_identical_end(is_ident, i);
	}
	if (tree.share_only_root(deepest_head)) {
		const int num_nodes = tree.get_num_nodes();
		for (int i = deepest_head; i < num_nodes; ++i) {
			if (tree.can_be_added(i)) {
				const label_value_type begin_atom = tree.begin_atom_label(is_ident, i);
				for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
					if (tree.remain(atomi)) {					  
					  if(atomi==0 && tree.get_label(i)==0){
					    adjacent_benzene++;
					  }
					  num += add_node_and_next(tree, is_ident, num_nodes, i, atomi,adjacent_benzene,outputfile,aut_group,str_lack_valence, str_smi);
					  tree.del_last_node();
					  if(atomi==0 && tree.get_label(i)==0){
					    adjacent_benzene--;
					  }
					}
				}
			}
			tree.update_identical_end(is_ident, i);
		}
	}
	return num;
}

inline size_t add_root_child_and_next(ChemTreeCenter& tree, is_ident_type is_ident, const label_value_type atom_label,int & adjacent_benzene, std::ofstream& outputfile, const std::vector< auto_group > &aut_group, const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_smi)
{	using namespace std;

	if (tree.add_root_child(is_ident, atom_label)) {
		const int type = tree.is_normal(1);
	
		if (type > 0) {
			if (num_lack_H > 0) {
			  return generate_multi(tree, is_ident, type, outputfile, aut_group, str_lack_valence, str_smi);
			} else {
			  return enum_substructure(tree, outputfile, aut_group, str_lack_valence, str_smi);
			}
		} else {
			return 0;
		}
	}
	size_t num = 0;
	
	if (tree.can_be_added_root()) {
		const label_value_type begin_atom = tree.begin_atom_label_root();
		for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
			if (tree.remain(atomi)) {
			  if(atomi<num_special_atom && tree.get_label(0)<num_special_atom){
			    adjacent_benzene++;
			  }

			  num += add_root_child_and_next(tree, is_ident, atomi, adjacent_benzene, outputfile, aut_group, str_lack_valence, str_smi);
			  tree.del_last_node();
			  if(atomi==0 && tree.get_label(0)==0){adjacent_benzene--;}
			}
		}
	}

	if (tree.share_only_root(1)) { 	
		const int num_nodes = tree.get_num_nodes();
		for (int i = 1; i < num_nodes; ++i) {
			if (tree.can_be_added(i)) {
				const label_value_type begin_atom = tree.begin_atom_label(is_ident, i);
				for (label_value_type atomi = begin_atom; atomi != first_atom_valence_one; ++atomi) {
				  	
					if (tree.remain(atomi)) {
					  if(atomi<num_special_atom && tree.get_label(0)<num_special_atom){adjacent_benzene++;}
					  num += add_node_and_next(tree, is_ident, num_nodes, i, atomi,adjacent_benzene,outputfile,aut_group,str_lack_valence, str_smi);
						tree.del_last_node();
						if(atomi==0 && tree.get_label(0)==0){adjacent_benzene--;}
					}
				}
			}
			tree.update_identical_end(is_ident, i);

		}
	}
	return num;
}

inline size_t add_root_and_next(ChemTreeCenter& tree, const label_value_type atom_label,std::ofstream& outputfile,const std::vector< auto_group > &aut_group, const std::vector< std::valarray<int> > &str_lack_valence, const std::vector< std::string > &str_smi)
{
	if (tree.add_root(atom_label)) {
		if (do_print) {
		  tree.write(outputfile, str_smi);
		  outputfile << "\n";
		}
		return 1;
	} 

	is_ident_type is_ident(num_except_H);
	is_ident[0] = (char)(-1);

	size_t num = 0;

	for (label_value_type atomi = 0; atomi != first_atom_valence_one; ++atomi) {
		if(tree.remain(atomi)){
		  int temp = 0;
		  if(atomi < num_special_atom && tree.get_label(0) < num_special_atom){ temp=1; }
		  num += add_root_child_and_next(tree, is_ident, atomi, temp, outputfile, aut_group, str_lack_valence, str_smi);
		  tree.del_last_node(); 
		}
	}
	return num;
}

inline static void chk_H(const std::vector<int>& atom_numbers,const int &bond_numbers, const int num_str)
{
	using namespace std;

	int withoutH_deg = 0;
	int withoutH_size = 0;
	int H_size = 0;
	for(int i = 0;i != first_atom_valence_one;i++){
		withoutH_deg += atom_numbers[i]*valence[i];
		withoutH_size += atom_numbers[i];
	}
	const int needH = withoutH_deg - 2*(withoutH_size-1);
	
	for(int i = first_atom_valence_one;i != num_distinct_atoms + num_str; i++){
		H_size += atom_numbers[i];
	}

	num_except_H = withoutH_size;
//num_lack_H = needH-H_size;//-(bond_numbers[1]*2);
	num_lack_H = bond_numbers;
	
	num_H = H_size;
	if (num_lack_H % 2 != 0) {
		cerr << "error # hydrogen atoms" << endl;
		exit(EXIT_FAILURE);
	}
}

inline size_t generate(const std::vector<int>& atom_numbers,const std::vector< auto_group > &group, 
const int num_str, const std::vector< std::string > str_atomlabels, const std::vector< std::vector< std::tuple<int,int,int> > > &str_bonds, 
const std::vector<int> str_numcycles,  const std::vector< std::valarray<int> > &str_lack_valence, const std::vector<std::string> &str_name , 		      
const std::vector<std::string> &str_smi ,const bool _do_print = false, const int round_num = 0)
{	
        using namespace std;
        //group = automorphism group of substructures
	//str_atomlabels = string, each string is atom label of substructure
	//str_bonds = vector of edges, each edge = {node node type_of_bond}
	//str_numcycles = #cycles 
	//str_lack_valence = valarray of available electron for each atom, in the same order as str_atomlabels
        do_print = _do_print;
	//num_special_atom is a global variable 
	num_special_atom = num_str;
	first_atom_valence_one += num_special_atom;

	ofstream outputfile;
	outputfile.open("output.smi");

	int lack_valence = 2;
	for(size_t i=0;i<atom_numbers.size();i++){
	  lack_valence+=input_valence[i]*atom_numbers[i];
	  lack_valence-=2*atom_numbers[i];
	}
	
	//initialize input for enumerating tree
	vector< vector<int> > atom_set; //input including special atom
	vector< int > bond_set; //bond_set is a vector of lack valence corresponding to input_atom_set

	vector<int> temp_number(atom_numbers);
	temp_number.insert( temp_number.begin(), num_str , 0 ); //add special atom at the beginning

	//insert into vector
	atom_set.push_back(temp_number);
	bond_set.push_back(lack_valence);

	if(num_str > 0){
	  std::vector<int> temp ( atom_numbers.size(), 0 );
	  //str_atom_number = vector of number of atom in substructure corresponding to atom label in atom_numbers
	  std::valarray< std::vector<int> > str_atom_number(temp, num_str);
	  //str_bond_number = vector of number of bond required for each substructure
	  std::vector<float> str_bond_number;
	  int label_name = 97;

	  for(size_t str_i = 0; str_i < num_str ; str_i++, label_name++){
	    //for each substructure -----

	    //add valence of new structure 
	    //valence is a global variable declared at the beginning of program
	    valence.insert( valence.begin()+str_i, str_lack_valence[str_i].sum() );
	    //add name of new structure to atomchar
	    wchar_t temp_name = label_name;
	    string name(1,temp_name);
	    atomchar.insert( atomchar.begin()+str_i, name );

	    float temp_bond_number = 2;

	    for(size_t i=0; i < str_atomlabels[str_i].size(); ++i){
	      //for each atom in substructure ------
	      string temp = string(1, str_atomlabels[str_i][i]);
	      int index = std::find(input_atomchar, input_atomchar+num_distinct_atoms, temp) - input_atomchar;
	      //index = index of atom [i] in substructure [str_i]

	      if(index == num_distinct_atoms){
		cout<<"Error: atoms in input structure are not appear in input atoms"<<endl;
		return 0;
	      }else{
		str_atom_number[str_i][index]++;

		if(str_atom_number[str_i][index] > atom_numbers[index]){
		  cout<<"Error: number of atoms is not enough for input structure"<<endl;
		  return 0;
		}

		temp_bond_number += input_valence[index]-2;
	      }
	    }

	    for(size_t i=0; i<str_lack_valence[str_i].size(); i++){
	      temp_bond_number -= str_lack_valence[str_i][i];
	    }
	    str_atom_number[str_i].insert(str_atom_number[str_i].begin(), num_str, 0 ); 
	    //add number of special atom
	    
	    /*for(size_t i=0; i<str_bonds[str_i].size(); i++) {
	      temp_bond_number += 2 * ( get_degree(std::get<2>(str_bonds[str_i][i]) )-1);
	      std::cout<<temp_bond_number<<std::endl;
	    }
	    std::cout<<"-----"<<std::endl;
	    temp_bond_number += 2 * str_numcycles[str_i];
	    std::cout<<temp_bond_number<<std::endl;*/
	    str_bond_number.insert( str_bond_number.begin(), temp_bond_number );
	    
	    count_structure(str_atom_number[str_i], str_bond_number[str_i], atom_set, bond_set, str_i);
	  }//for each substructure
	  
	}
	
	size_t num = 0;
	
	for(size_t i=0; i < atom_set.size(); i++){
	  
	  if(round_num > (signed)atom_set.size() ){
	     cout<<"number of round is too large"<<endl;
	    return 0;
	  }
	  if(round_num > 0){
	    i = round_num - 1;
	  }
	        cout<<"==================="<<endl;
		for(size_t j=0; j<num_distinct_atoms+num_str; j++){
		  cout<<"amount of "<<atomchar[j];
		  cout<<": "<<atom_set[i][j];
		  if(j < num_str)
		    cout<<"  ("<<str_name[j]<<") ";
		  cout<<endl;
		  /*if(j < num_str){
		    print_autgroup(group[j]);
		    for(int k = 0;k<str_lack_valence[j].size(); k++)
		      std::cout<<str_lack_valence[j][k]<<" ";
		    std::cout<<std::endl;
		    std::cout<<valence[j]<<std::endl;
		  }*/
		}
		cout<<"-----------------"<<endl;
		cout<<"#lack valence:"<<bond_set[i]<<endl;
		cout<<"-----------------"<<endl;
		t_r.clear();
		t_v.clear();
		location.clear();
		lack_degree.clear();

		chk_H(atom_set[i], bond_set[i], num_str);
		t_r.reserve(num_except_H>>1);
		t_v.reserve(num_except_H>>1);
		location.reserve(num_except_H>>1);
		lack_degree.reserve(num_except_H);

		ChemTreeCenter tree(atom_set[i], num_except_H);
		
		for (label_value_type atomi = 0; atomi != first_atom_valence_one; atomi++) {
			if (tree.remain(atomi)){
			  num += add_root_and_next(tree, atomi, outputfile, group, str_lack_valence, str_smi);
			  tree.del_root();
			}
		}

		std::cout<<"accumulated #enumerated structure until "<<i+1<<" round is "<<num<<std::endl;
		if(round_num > 0){
		  break;
		}
	}

	outputfile.close();
	return num;
}

#if 0
void ChemTreeCenter::printmol(const std::string& filename) const
{
using namespace std;

ofstream ofs((filename+".mol").c_str(),ios_base::out | ios_base::app);
ofs << "\n";
ofs<<"   RDKit          2D\n";
ofs << "\n";
ofs << boost::format("%3d%3d") % num_nodes % (num_nodes-1);
ofs<<"  0  0  0  0  0  0  0  0999 V2000\n";
for(int i = 0;i<num_nodes;i++){
ofs<<"    0.0000    0.0000    0.0000 "<<atomchar[nodes[i].label]<<"   ";
ofs<<"0  0  0  0  0  0  0  0  0  0  0  0\n";
}

for (int i = 0; i != num_nodes; ++i) {
for (int j = 0; j != nodes[i].num_children; ++j) {
ofs << boost::format("%3d%3d%3d") % (i+1) % (nodes[i].children[j]+1) % nodes[nodes[i].children[j]].multi;
ofs <<"  0\n";
}
}
ofs << "M END\n$$$$\n";
ofs.close();
}
#endif

inline void ChemTreeCenter::printsmi(const int i) const
{
	if (i != 0) {
		std::cout << bondchar[nodes[i].multi];
	}
	std::cout << atomchar[nodes[i].label];
	const valence_value_type nc = nodes[i].num_children;
	if (nc > 1) {
		for (valence_value_type v = 0; v != nc; ++v) {
			std::cout << "(";
			printsmi(nodes[i].children[v]);
			std::cout << ")";
		}
	} else if (nc == 1) {
		printsmi(nodes[i].children[0]);
	}
}

inline void ChemTreeCenter::printsmi_single(const int i) const
{
	std::cout << atomchar[nodes[i].label];
	const valence_value_type nc = nodes[i].num_children;
	if (nc > 1) {
		for (valence_value_type v = 0; v != nc; ++v) {
			std::cout << "(";
			printsmi_single(nodes[i].children[v]);
			std::cout << ")";
		}
	} else if (nc == 1) {
		printsmi_single(nodes[i].children[0]);
	}
}

inline void ChemTreeCenter::printseq() const
{
	using namespace std;

	cout << nodes[0].label << nodes[0].num_children;
	for (int i = 1; i!= num_nodes; ++i) {
		cout << nodes[i].label << nodes[i].multi << nodes[i].num_children;
	}
}

inline void ChemTreeCenter::printseq_single() const
{
	using namespace std;

	cout << nodes[0].label << nodes[0].num_children;
	for (int i = 1; i!= num_nodes; ++i) {
		cout << nodes[i].label << nodes[i].num_children;
	}
}

inline void ChemTreeCenter::print() const
{
	printsmi();
	std::cout << std::endl;
}

inline void ChemTreeCenter::print_single() const
{
	printsmi_single(0);
	std::cout << std::endl;
}


class smiles{
  std::vector<char> atom;
  std::vector<int> cycle,position;

public:
  smiles(std::vector<char> atom, std::vector<int> cycle, std::vector<int> position){
    this->atom = atom;
    this->cycle = cycle;
    this->position = position;
  }
  smiles(){}
    int get_size(){
    return atom.size();
  }

  char get_atom(size_t index){
    if(index < atom.size()){
      return atom[index];
    }else{
      std::cout << "get_atom: index out of bound" <<std::endl;
      exit(0);
      return '\0';
    }
  }

  int get_cycle(size_t index){
    if(index < cycle.size()){
      return cycle[index];
    }else{
      std::cout << "get_cycle: index out of bound" <<std::endl;
      exit(0);
      return -1;
    }
  }

  int get_position(size_t index){
    if(index < position.size()){
      return position[index];
    }else{
      std::cout << "get_position: index out of bound" << std::endl;
      exit(0);
      return -1;
    }

  }

  void push(char a, int c, int p){
    atom.push_back(a);
    cycle.push_back(c);
    position.push_back(p);
  }

  void push(smiles smi, size_t index){
    if(index < smi.get_size()){
      this->push( smi.get_atom(index), smi.get_cycle(index), smi.get_position(index) );
    }else{
      std::cout << "push_smiles(1): index out of bound" << std::endl;
      exit(0);
    }
  }

  void push(smiles smi, size_t start, size_t end){
    if(start > smi.get_size() || end > smi.get_size() ){
      std::cout << "push_smiles(2): index out of bound" << std::endl;
      exit(0);
    }

    if(start > end){
      //reverse order
      for(size_t index = start; index !=end; index--){
	this->push( smi.get_atom(index), smi.get_cycle(index), smi.get_position(index) );
      }
      this->push( smi.get_atom(end), smi.get_cycle(end), smi.get_position(end) );
    }else{
      //forward order
      for(size_t index = start; index !=end; index++){
	this->push( smi.get_atom(index), smi.get_cycle(index), smi.get_position(index) );
      }
      this->push( smi.get_atom(end), smi.get_cycle(end), smi.get_position(end) );
    }
  }

  void push_openbracket(){
    this->push( '(', 0, -1);
  }

  void push_closebracket(){
    this->push( ')', 0, -1);
  }

  size_t leftmost_concat_bracket(std::vector< std::pair<int,int> > & bracket_set, size_t start){
    //start = index of the open bracket -> return index of the leftmost open bracket concatenated to start
    if(atom[start] != '('){
      std::cout<<"concatenate bracket: error"<<std::endl;
      exit(0);
    }

    if(start == 0 || atom[start-1] != ')'){
      return start;
    }

    size_t result = start;
    for(size_t i = bracket_set.size()-1; i < bracket_set.size(); i--){
      if(bracket_set[i].second == result-1){
	result = bracket_set[i].first;
      }
    }
    return result;
  }

  void remove_empty_brackets(){
    for(size_t i =0; i < atom.size()-1; i++){
      if( atom[i] == '(' && atom[i+1] == ')' ){
	atom.erase( atom.begin()+i, atom.begin()+i+2 );
	cycle.erase( cycle.begin()+i, cycle.begin()+i+2 );
	position.erase( position.begin()+i, position.begin()+i+2 );

	i--; 
      }
    }
  }

  void print(){
    using namespace std;
    cout<<"atom     = ";
    for(size_t i=0; i< atom.size(); i++){
      cout<<atom[i]<<" ";
    }
    cout<<endl;
    cout<<"position = ";
    for(size_t i=0; i< position.size(); i++){
      if(position[i] >= 0)
	cout<<position[i]<<" ";
      else
	cout<<"/";
    }
    cout<<endl;
    cout<<"cycle    = ";
    for(size_t i=0; i< cycle.size(); i++){
      cout<<cycle[i]<<" ";
    }
    cout<<endl;
  }
};

void reverse_smiles( int start, int end, smiles &input_smi, smiles &new_smi, std::vector< std::pair<int,int> > &bracket_set){
  using namespace std;

  pair<int,int> innermost(-1,-1);
  for(size_t i=0; i<bracket_set.size() ; i++){
    if(bracket_set[i].first >= start && bracket_set[i].second <= end){
      innermost = bracket_set[i];
    }
  }
  if(innermost.first == -1){
    if(end >= 0){
      new_smi.push( input_smi, end, start);
    }
  }else{
    if(innermost.second+1 <= end)
      reverse_smiles( innermost.second+1, end, input_smi, new_smi, bracket_set );
    size_t end_part1 = input_smi.leftmost_concat_bracket(bracket_set, innermost.first) - 1;
    //if there are concatenated brackets -> end_part1 = the part before the leftmost concatenated bracket
    /*if( input_smi.get_atom(end_part1) == ')' ){
      if(bracket_set.size() > 0){
	for(size_t i = bracket_set.size()-1; i < bracket_set.size(); i--){
	  if(bracket_set[i].second == end_part1){
	    end_part1 = bracket_set[i].first-1;
	  }
	}
      }
      }*/
    new_smi.push( input_smi, end_part1);
    
    new_smi.push_openbracket();
    
    reverse_smiles( start, end_part1-1, input_smi, new_smi, bracket_set);
    new_smi.push_closebracket();

    if(end_part1 != innermost.first-1){
      //if there are concatenated brackets -> add them after reverse of part1
      new_smi.push( input_smi, end_part1+1, innermost.first-1);
    }

    new_smi.push( input_smi, innermost.first+1, innermost.second-1);
  }

}

void arrange_smiles(int index, smiles &input_smi, smiles &new_smi){

  using namespace std;

  vector< pair<int,int> > bracket_set;
  for(size_t i = 0; i < input_smi.get_size(); i++){
    if( input_smi.get_atom(i) == '(' ){
      if(i <= index) index++;

      pair<int,int> temp;
      temp = std::make_pair(i,-1);
      bracket_set.push_back(temp);
    }else if( input_smi.get_atom(i) == ')' ){
      if(i <= index) index++;

      for(int j = bracket_set.size()-1; j>=0 ; j-- ){
	if(bracket_set[j].second == -1){
	  bracket_set[j].second = i;
	  break;
	}
      }
    }
  }
  
  //find index innermost brackets where index is between them (-1 if no such brackets)
  int index_innermost_unpair = -1;
  for(size_t i=0; i<bracket_set.size(); i++){
    if(bracket_set[i].second > index && bracket_set[i].first < index){
      index_innermost_unpair = i;
    }
  }

  if(index_innermost_unpair == -1){
    new_smi.push( input_smi, index);
    new_smi.push_openbracket();
    
    if(index-1 >= 0)
      reverse_smiles(0, index-1, input_smi, new_smi, bracket_set);
    new_smi.push_closebracket();
    
    if(index+1 <= input_smi.get_size()-1)
      new_smi.push( input_smi, index+1, input_smi.get_size()-1);
  }else{ //index_innermost_unpair != -1
    
    pair<int,int> innermost = bracket_set[index_innermost_unpair];

    new_smi.push( input_smi, index);

    new_smi.push_openbracket();
    if(index+1 <= innermost.second-1)
      new_smi.push( input_smi, index+1, innermost.second-1);
    new_smi.push_closebracket();

    if(index-1 >= innermost.first+1)
      reverse_smiles(innermost.first+1, index-1, input_smi, new_smi, bracket_set);
    
    size_t end_part1 = input_smi.leftmost_concat_bracket(bracket_set, innermost.first) - 1;
    /*//if there are concatenated brackets -> end_part1 = the part before the leftmost concatenated bracket
    if( input_smi.get_atom(end_part1) == ')' ){
      if(bracket_set.size() > 0){
	for(size_t i = bracket_set.size()-1; i < bracket_set.size(); i--){
	  if(bracket_set[i].second == end_part1){
	    end_part1 = bracket_set[i].first-1;
	  }
	}
      }
      }*/
    new_smi.push(input_smi, end_part1);

    new_smi.push_openbracket();
    
    //find unpair open bracket from start to end_part1
    vector<int> index_unpair_bracket;
    for(size_t i = 0; i < index_innermost_unpair; i++){
      if(bracket_set[i].first < index && bracket_set[i].second > index){
	index_unpair_bracket.push_back(i);
      }
    }

    size_t leftmost_unpair = innermost.first;
    if(index_unpair_bracket.size() > 0){
      leftmost_unpair = bracket_set[ index_unpair_bracket[0] ].first;

      for(int i = index_unpair_bracket.size()-1; i >= 0 && i < index_unpair_bracket.size(); i--){
	size_t start = bracket_set[ index_unpair_bracket[i] ].first + 1;
	size_t end;
	if( i == index_unpair_bracket.size()-1 ){
	  end = innermost.first;
	}else{
	  end = bracket_set[ index_unpair_bracket[i+1] ].first;
	}
	end = input_smi.leftmost_concat_bracket(bracket_set, end);
	/*if( input_smi.get_atom( end ) == ')' ){
	  if(bracket_set.size() > 0){
	    for(size_t i = bracket_set.size()-1; i < bracket_set.size(); i--){
	      if(bracket_set[i].second == end ){
		end = bracket_set[i].first-1;
	      }
	    }
	  }
	  }*/
	end-=2;	  
	
	if(start <= end)
	  reverse_smiles(start, end, input_smi, new_smi, bracket_set);		
	//find concatenate bracket
	size_t begin_concatenate = input_smi.leftmost_concat_bracket( bracket_set, bracket_set[ index_unpair_bracket[i] ].first);
	if( input_smi.get_atom( bracket_set[ index_unpair_bracket[i] ].first - 1) == ')' ){
	  /*if(bracket_set.size() > 0){
	    for(size_t j = bracket_set.size()-1; j < bracket_set.size(); j--){
	      if(bracket_set[j].second == begin_concatenate - 1){
		begin_concatenate = bracket_set[j].first;
	      }
	    }
	    }*/
	  if(i==0){
	    leftmost_unpair = begin_concatenate;
	  }
	}
	//push the atom left of begin_concatenate
	new_smi.push( input_smi, begin_concatenate - 1);
	//push concatenate with brackets
	if(begin_concatenate != bracket_set[ index_unpair_bracket[i] ].first ){
	  new_smi.push( input_smi, begin_concatenate, bracket_set[ index_unpair_bracket[i] ].first -1 );
	}

	//if(i != 0){
	new_smi.push_openbracket();
	//}
      }
    }

    //new_smi.push_openbracket();

    //push reverse of [start: leftmost unpair bracket]
    /*//find concatenate bracket
    size_t begin_concatenate = leftmost_unpair;
    if( input_smi.get_atom(begin_concatenate-1) == ')' ){
      if(bracket_set.size() > 0){
	for(size_t i = bracket_set.size()-1; i < bracket_set.size(); i--){
	  if(bracket_set[i].second == begin_concatenate - 1){
	    begin_concatenate = bracket_set[i].first;
	  }
	}
      }
    }
    //push the atom left of begin_concatenate
    new_smi.push( input_smi, begin_concatenate - 1);
    //push concatenate with brackets
    if(begin_concatenate != bracket_set[ index_unpair_bracket[i] ].first ){
      new_smi.push( input_smi, begin_concatenate, bracket_set[ index_unpair_bracket[i] ].first -1 );
      }*/
    if(static_cast<int>(leftmost_unpair)-2 >= 0)
      reverse_smiles( 0, leftmost_unpair-2, input_smi, new_smi, bracket_set);

    new_smi.push_closebracket();

    size_t end_innermost = input_smi.get_size();
    if(index_unpair_bracket.size() > 0){
      end_innermost = bracket_set[ index_unpair_bracket[index_unpair_bracket.size()-1] ].second;

      for(size_t i = 0; i < index_unpair_bracket.size(); i++){
	size_t start = bracket_set[ index_unpair_bracket[i] ].second + 1;
	size_t end;
	if( i == 0 ){
	  end = input_smi.get_size()-1;
	}else{
	  end = bracket_set[ index_unpair_bracket[i-1] ].second - 1;
	}	
	
	if(start <= end)
	  new_smi.push( input_smi, start, end);
	//if(i != index_unpair_bracket.size() - 1 ){
	new_smi.push_closebracket();      
	//}
      }
    }
    
    if(end_part1 != innermost.first-1){
      new_smi.push_openbracket();
    }

    if(innermost.second+1 <= end_innermost-1)
      new_smi.push( input_smi, innermost.second+1, end_innermost-1);    

    if(end_part1 != innermost.first-1){
      new_smi.push_closebracket();
    }

    if(end_part1 != innermost.first-1){
      //if there are concatenated brackets put them at the end
      new_smi.push( input_smi, end_part1+1, innermost.first-1);
    }
    
  }

  new_smi.remove_empty_brackets();

}

inline void ChemTreeCenter::write_smiles(int index, std::ofstream & file, int & num_cycle, const std::vector<std::string> & str_smi) const{
  
  if(index !=0){
    file << bondchar[ nodes[index].multi ];
   }

  if(nodes[index].label < num_special_atom){
    //current node is an input structure
    std::string smi = str_smi[ nodes[index].label ];
    std::vector<char> atom;
    std::vector<int> cycle;
    std::vector<int> position;
    std::vector< std::pair<int,int> > bracket_set;
    int step = 0;

    //initialize three vectors
    for(size_t i=0; i < smi.size(); i++){
      if( !isdigit(smi[i]) ){
	cycle.push_back( 0 );
	atom.push_back( smi[i] );
	if( !ispunct(smi[i]) ){
	  position.push_back( step );
	  step++;
	}else{
	  position.push_back( -1 );
	}
      }else{
	cycle[ cycle.size()-1 ] = atoi( &smi[i] );
      }
    }        

    smiles input_smi(atom, cycle, position);

    int index_cycle = num_cycle; //num_cycle start with 0
    num_cycle+= cycle[ std::max_element(cycle.begin(), cycle.end()) - cycle.begin() ];
    
    int begin_cp = 0;
    if(index != 0){
      begin_cp = nodes[index].bond_position[0];
    }

    //begin_cp = position bonding with parent
    //if begin_cp != 0 then need to rearrange smiles so that corresponding atom is the first char in smiles
    if(begin_cp != 0){
      smiles new_smi;
      arrange_smiles(begin_cp, input_smi, new_smi);
      input_smi = new_smi;
    }

    //print from begin_cp to end_of_smiles
    for(size_t i = 0; i < input_smi.get_size() ; i++){
      
      file << input_smi.get_atom(i);

      if( input_smi.get_position(i) != -1 ){
	
	int nth = std::find( nodes[index].bond_position.begin(), nodes[index].bond_position.end(), input_smi.get_position(i) ) - nodes[index].bond_position.begin();
	nth--;

	if(nth < nodes[index].num_children && nth >= 0){
	  file << "(";
	  write_smiles( nodes[index].children[nth], file, num_cycle, str_smi);
	  file << ")";
	}
      }
      if( input_smi.get_cycle(i) != 0 ){
	file << ( index_cycle + input_smi.get_cycle(i) );
      }
    }

    /*
    //print from begin_of_smiles to begin_cp
    for(int i=0; i < begin_cp; i++){
      file << atom[i];

      int nth = std::find( nodes[index].bond_position.begin(), nodes[index].bond_position.end(), i) - nodes[index].bond_position.begin();
      nth--;
      
      if(nth < nodes[index].num_children && nth >= 0){
	file << "(";
	write_smiles( nodes[index].children[nth], file, num_cycle, str_smi);
	file << ")";
      }

      if(cycle[i] != 0){
	file << (index_cycle + cycle[i]);
      }
    }
    */

  }else{
    //current node is a normal atom (not input structure)
    file << atomchar[ nodes[index].label ];

    if(nodes[index].num_children > 0){
      for(int i = 0; i < nodes[index].num_children-1; i++){
	file << "(";
	write_smiles( nodes[index].children[i], file, num_cycle, str_smi);
	file << ")";
      }

      write_smiles( nodes[index].children[ nodes[index].num_children-1 ], file, num_cycle, str_smi);
    }
  }

}


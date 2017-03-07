#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include "bfsenum.hpp"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

#include <time.h>
#include <curl/curl.h>

using namespace std;

bool MolToSmiles(string imol, string &smi){

  ifstream ifs(imol);
  ostringstream oss;
  if(!ifs){
    cout << "Cannot open input file\n";
    return false;
  }

  OpenBabel::OBConversion conv(&ifs, &oss);
  conv.SetInAndOutFormats("MOL", "CAN");
  
  conv.Convert();
  smi = oss.str();
  std::cout<<" smiles:"<<smi<<std::endl;
  istringstream iss(smi);
  ofstream ofs(imol);

  OpenBabel::OBConversion rev(&iss, &ofs);
  rev.SetInAndOutFormats("SMI", "MOL");
  rev.Convert();

  size_t pos = smi.find("\t");
  smi = smi.substr(0, pos);


  OpenBabel::OBConversion ob;
  ob.SetInFormat("MOL");
  OpenBabel::OBMol mol;

  ob.ReadFile(&mol, imol);
  
  // convert kekule representation to aromatic bond
  OpenBabel::OBBond *bond;
  vector<OpenBabel::OBEdgeBase*>::iterator it;
  for( bond = mol.BeginBond(it); bond; bond = mol.NextBond(it) ){
    if( bond->IsAromatic() ){
      bond->SetBO(4);
    }
  }

  ob.SetOutFormat("MOL");
  ob.WriteFile(&mol, imol);
    
  return true;
}

int main(int argc, char **argv)
{
  namespace po = boost::program_options;

  po::options_description opt("Allowed options");

  opt.add_options()
    ("c,c", po::value<int>(), ": input # C atoms")
    ("n,n", po::value<int>(), ": input # N atoms")
    ("o,o", po::value<int>(), ": input # O atoms")
    ("h,h", po::value<int>(), ": input # H atoms")
    ("r,r", po::value<int>(), ": input round number")
    ("input,s", po::value< vector<string> >()->multitoken(), ": input file")
    ("output,t", ": output results")
    ("help,p", ": show this help message");


  po::variables_map argmap;
  po::store(po::parse_command_line(argc, argv, opt), argmap);
  po::notify(argmap);

  if (argmap.count("help") 
      or ((not argmap.count("c")) and(not argmap.count("n")) and(not argmap.count("o")) and (not argmap.count("h")))) {
    std::cerr << "Usage: " << argv[0] << " [option]" << std::endl << opt << std::endl;

    return EXIT_FAILURE;
  }

  int num_C = 0, num_N = 0, num_O = 0, num_H = 0;
  int round_num = 0;

  if (argmap.count("c")) {
    num_C = argmap["c"].as<int>();
  }
  if (argmap.count("n")) {
    num_N = argmap["n"].as<int>();
  }
  if (argmap.count("o")) {
    num_O = argmap["o"].as<int>();
  }
  if (argmap.count("h")) {
    num_H = argmap["h"].as<int>();
  }
  if (argmap.count("r")) {
    round_num = argmap["r"].as<int>();
  }

  std::vector< std::vector< std::vector<int> > > group;
  std::vector<string> str_atomlabels;
  std::vector< vector< std::tuple<int,int,int> > > str_bonds;
  std::vector<int> str_numcycles;
  std::vector< valarray<int> > str_lack_valence;
  std::vector<string> str_name;
  std::vector<string> str_smi;
  int num_str = 0;

  if(argmap.count("input")) {
    vector<string> input_files = argmap["input"].as<vector<string> >();
    num_str = input_files.size();
    std::vector< std::vector<int> > temp;
    group.insert(group.begin(),num_str, temp);
    
    for(size_t i = 0; i<num_str; i++) {
      string temp_smi;
      if( MolToSmiles(input_files[i], temp_smi) ){
	str_smi.push_back(temp_smi);
      }else{
	std::cout << "Cannot open input file" <<std::endl;
	return 1;
      }

      Mol mol;
      str_name.push_back(input_files[i]);
      mol.read(input_files[i]);
      mol.aut(group[i]);

      str_atomlabels.push_back( mol.get_atomlabels() );
      
      std::vector< std::tuple<int,int,int> > temp_bonds;
      mol.get_edges( temp_bonds );
      str_bonds.push_back( temp_bonds );

      str_numcycles.push_back( mol.get_numcycles() );
      str_lack_valence.push_back( mol.get_lack_valence() );
    }  
  }
  
  vector<int> atom_numbers(4);
  atom_numbers[0] = num_C;
  atom_numbers[1] = num_N;
  atom_numbers[2] = num_O;
  atom_numbers[3] = num_H;
  
  /*
    std::vector<char> a{'c','(','c','c','c',')','(','c','c',')','C','c','c'};    
  a = {'c','(','c','c',')','(','c','c','(','c','c','c',')','c',')','c'};
  std::vector<int> b{0, -1, 1, 2, 3, -1, -1, 4, 5, -1, 6, 7, 8};     
  b= {0, -1, 1, 2, -1, -1, 3, 4, -1, 5, 6, 7, -1, 8, -1, 9};
  std::vector<int> c{0,  0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
  smiles smi(a,c,b);
  smiles result;
    arrange_smiles(10, smi, result);

  for(size_t i=0; i<result.get_size(); i++){
    std::cout<<result.get_atom(i)<<"\t"<<result.get_cycle(i)<<"\t"<<result.get_position(i)<<std::endl;
  }
  */


  size_t num = generate(atom_numbers, group, num_str, str_atomlabels, str_bonds, str_numcycles, str_lack_valence, str_name, str_smi, argmap.count("output"), round_num);
  cout<<"======Result======="<<endl;
  cout<<"Total #enumerated structures:" << num << endl;
  return EXIT_SUCCESS;
}


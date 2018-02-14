#include <Rcpp.h>
using namespace Rcpp ;



///////////////////////////////////////////
//
//  Program: PedKin version
const char* version = "1.5";
//   Author: Bonnie Kirkpatrick
//  Company: Intrepid Net Computing
//     Date: May 2, 2017
//  License: GNU GPL v2.0
//
//
//  Given a pedigree of individuals, founder kinship coefficients, 
//  and the individuals of interest compute the kinship on all
//  pairs of individuals individuals of interest in the pedigree.
//  Uses the founder kinship coefficients to initialize the kinship.
//
//  Outputs kinship matrix with the diagonal transformed into the 
//  inbreeding coefficients.
//
//  Algorithm properties:
//    Exact Algorithm
//      O(n^2) time for n individuals in a pedigree
//      O(n^2) space
//    Approximate Algorithm
//      O(n) time for n individuals and sqrt(n) individuals of interest
//      O(n) space
//
///////////////////////////////////////////


#include <vector>
#include <algorithm>
#include <string>
#include <set>
#include <queue>
#include <stack>

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>


using namespace std;



#define FALSE 0
#define TRUE 1


#define UNKNOWN 0

#define MALE 1
#define FEMALE 2
#define SAME_AS 10

#define AFFECTED 2
#define UNAFFECTED 1

#define FIRST_ALLELE 1
#define SECOND_ALLELE 2


#define INVALID -1
#define ERROR -55


#define DOUBLEDIGIT -11
#define INVALID_ALLELE -22
#define OUTOFBOUNDS -33


//#define DEBUG 1  // comment out to remove debugging







typedef int Gender;
typedef int Affection;



class Marker
{
 public:
  Marker(){
    allele1 = UNKNOWN;
    allele2 = UNKNOWN;
  }

  // make sure alleles are sorted when storing them as markers
  Marker(int a, int b){
    if (a <= b)
      {
	allele1 = a;
	allele2 = b;
	if (allele1 == UNKNOWN && allele2 == FIRST_ALLELE)
	  {
	    allele1 = b;
	    allele2 = a;
	  }
      }
    else
      {
	allele1 = b;
	allele2 = a;	
	if (allele1 == UNKNOWN && allele2 == FIRST_ALLELE)
	  {
	    allele1 = a;
	    allele2 = b;
	  }
      }
  }
  Marker(const Marker& m){
    allele1 = m.allele1;
    allele2 = m.allele2;
  }
  const Marker& operator= (const Marker& m) {
    allele1 = m.allele1;
    allele2 = m.allele2;
    return *this;
  }

  // insert an allele in the case that there can be at most one copy of each type of allele
  int insertAllele(int a)
  {
    if (a <= UNKNOWN)
      {
	return INVALID_ALLELE;
      }
    if (allele1 != UNKNOWN && allele2 != UNKNOWN 
	&& (a != allele1 || a != allele2))
      {
	return INVALID_ALLELE;
      }
    if (a == allele1 || a == allele2)
      return 0;
    
    if (allele1 == UNKNOWN && allele2 == UNKNOWN)
      {
	if (a == FIRST_ALLELE)
	  allele1 = a;
	if (a == SECOND_ALLELE)
	  allele2 = a;
	return 1;
      }
    else
      {
	if (allele1 == UNKNOWN)
	  allele1 = a;
	else
	  allele2 = a;
	if (allele1 > allele2)
	  return ERROR;
	return 1;
      }
    return 0;
  }

  int allele1;
  int allele2;
};



class Individual
{
public:
  Individual(){
    mother_id = INVALID;
    father_id = INVALID;
    individual_index = INVALID;
    mother = NULL;
    father = NULL;
  }
  Individual(int fid, int id, int f, int m, Gender s, Affection a, long p, int l){
    family_id = fid;
    individual_id = id;
    individual_index = INVALID;
    mother_id = m;
    father_id = f;
    mother = NULL;
    father = NULL;
    sex = s;
    affection = a;
    file_position = p;
    line_number = l;
  }
  Individual(const Individual& i){
    family_id = i.family_id;
    individual_id = i.individual_id;
    individual_index = i.individual_index;
    father_id = i.father_id;
    mother_id = i.mother_id;
    mother = i.mother;
    father = i.father;
    sex = i.sex;
    affection = i.affection;
    isTyped = i.isTyped;
    file_position = i.file_position;
    line_number = i.line_number;
    disease_marker = disease_marker;
    for(unsigned int c = 0;  c < i.genotypes.size(); c++)
      {
	genotypes.push_back(i.genotypes[c]);
      }
    for(unsigned int c = 0;  c < i.children.size(); c++)
      {
	children.push_back(i.children[c]);
      }
  }
  const Individual& operator= (const Individual& i){
    family_id = i.family_id;
    individual_id = i.individual_id;
    individual_index = i.individual_index;
    father_id = i.father_id;
    mother_id = i.mother_id;
    mother = i.mother;
    father = i.father;
    sex = i.sex;
    affection = i.affection;
    isTyped = i.isTyped;
    file_position = i.file_position;
    line_number = i.line_number;
    disease_marker = disease_marker;
    for(unsigned int c = 0;  c < i.genotypes.size(); c++)
      {
	genotypes.push_back(i.genotypes[c]);
      }
    for(unsigned int c = 0;  c < i.children.size(); c++)
      {
	children.push_back(i.children[c]);
      }
    return *this;
  }
  bool operator== (const Individual& indiv) {
    if (family_id == indiv.family_id)
      return this->individual_id == indiv.individual_id;
    else
      return 0;
  }

  void setIndivIndex(int i){
    individual_index = i;
  }



  int family_id;
  int individual_id;
  int individual_index;
  int father_id;
  int mother_id;
  Individual* mother;
  Individual* father;
  Gender sex;
  Affection affection;
  int isTyped;

  long file_position;
  int line_number;


  vector<Individual*> children; 
  vector<Marker> genotypes;
  vector<Marker> identity;  // name of ancestor alleles that these were inherited from (i.e. IBD state)

  Marker disease_marker;
  Marker disease_identity; // disease SNP IBD state


  vector<bool> segregation_indicator;
  vector<int> cc_membership;
};



class OrderItem
{
 public:
  OrderItem(){
    individual_id = UNKNOWN;
    individual_index = INVALID;
    individual = NULL;
    depth = INVALID;

  }
  OrderItem(const OrderItem& o){
    individual_id = o.individual_id;
    individual_index = o.individual_index;
    individual = o.individual;
    depth = o.depth;
  }
  const OrderItem& operator= (const OrderItem& o) {
    individual_id = o.individual_id;
    individual_index = o.individual_index;
    individual = o.individual;
    depth = o.depth;
    return *this;
  }

  int individual_id;
  int individual_index;
  Individual* individual;
  int depth;
};
struct order_lt
{
  bool operator()(const OrderItem i, const OrderItem j) const
  {
    return i.depth > j.depth;
  }
};


class Family
{
public:
  Family(int id){ family_id = id; }

  bool operator== (const Family& fam) {
    return this->family_id == fam.family_id;
  } 
  bool operator== (vector<Family>::iterator fam) {
    return this->family_id == fam->family_id;
  }

  void setIndices(){
    for(unsigned int i = 0;  i < members.size();  i++){
      members[i]->setIndivIndex(i);
    }
  }


  vector<Individual*>::iterator find_individual(vector<Individual*>::iterator curr, 
						vector<Individual*>::iterator end, 
						int individual_id)
  {
    while (curr != end)
      {
	// check the id for this person
	if ((*curr)->individual_id == individual_id)
	  {
	    return (curr);
	  }
	curr++;
      }
    return end;
  }


  int family_id;
  vector<Individual*> members;

  vector<Individual*> founder; // index of the child in the vector of individuals
  vector<Individual*> nonfounder; // index of the child in the vector of individuals
  vector<Individual*> typed; // index of the child in the vector of individuals
  vector<Individual*> untyped; // index of the child in the vector of individuals

  vector<Individual*> indiv_interest; // list of individuals of interest

  vector<OrderItem> order;  // bottom up order


  vector<vector<double> > founder_kinship;
};



vector<Family> pedigrees;
vector<vector<int> > founder_haplotypes;





int read_allele (char* line, int length, int& index)
{
  while (index < length && isspace(*(line + index)) ){
    index++;
  }
  //printf("read_allele(len %i, idx %i) line[idx] = %c\n", length, index, *(line + index));

  if (index >= length){ 
    //printf("Returning OUTOFBOUNDS due to %i >= %i", index, length);
    return OUTOFBOUNDS;
  } 

  string a;
  a = a + *(line+index);
  int len = 1;
  while (!isspace(*(line+index+len)) && *(line+index+len) != '\0'){
    a = a + *(line+index+len);
    len++;
  }
  index+=len;



  int allele = atoi(a.c_str());


  // For valid SNP alleles
  if (len > 1){
    printf("double digit allele: |%s|\n", a.c_str());
    return DOUBLEDIGIT;
  }
  // check for valid SNP allele
  if (allele != UNKNOWN && allele != FIRST_ALLELE && allele != SECOND_ALLELE){
    return INVALID_ALLELE;
  }


  // check for arbitrary marker allele (i.e. non-SNP markers)
  //if (allele < UNKNOWN){
  //  return INVALID_ALLELE;
  //}

  while (index < length && isspace(*(line + index)) ){
    index++;
  }


  return allele;
}



void read_pedigree_file(const char* in_file)
{


  // read the graph  
   FILE* fp = fopen (in_file, "r");
   if (fp == NULL)
     {
      printf("\nERROR: Cannot open %s.\n\n", in_file);
      exit(-1);
   }



  char line [1024];



  char column[256];
  int line_number = 1;
  long file_position = ftell(fp);
  while ( fscanf (fp, "%s ", column) != EOF)
      {
	//printf("%i]", line_number);
	int famID = atoi(column);
	//printf(" %i", famID);

	fscanf (fp, "%s ", column);
	int indivID = atoi(column);
	//printf(" %i", indivID);

	fscanf (fp, "%s ", column);
	int faID = atoi(column);
	//printf(" %i", faID);

	fscanf (fp, "%s ", column);
	int moID = atoi(column);
	//printf(" %i", moID);

	int isFounder = FALSE;
	if (faID == UNKNOWN && moID == UNKNOWN)
	  {
	    isFounder = TRUE;
	  }


	fscanf (fp, "%s ", column);
	int sex = atoi(column);
	//printf(" %i", sex);

	fscanf (fp, "%s", column);
	int aff = atoi(column);
	//printf(" %i", aff);
	//printf("\n");
	//fflush(stdout);



	// Find the family or make it
	Family newfam(famID);
	vector<Family>::iterator fam;
	fam = find(pedigrees.begin(), pedigrees.end(), newfam);
	if (fam == pedigrees.end()){
	  // make family
	  pedigrees.push_back(newfam);
	  fam = pedigrees.end() - 1;
	}

	//printf(" newfam  ");
	//fflush(stdout);


	// insert individual into pedigree
	Individual* newperson = new Individual(famID, indivID, faID, moID, sex, aff, file_position, line_number);
	vector<Individual*>::iterator person;
	(fam->members).push_back(newperson);
	person = fam->members.end() - 1;
	if (isFounder){
	  (fam->founder).push_back(newperson);
	} else {
	  (fam->nonfounder).push_back(newperson);
	}




	//printf(" newper  ");
	//fflush(stdout);


	// Read all the alleles
	int isTyped = FALSE;
	
	
	char* notNull = fgets(line, 1024, fp);
	int length = strlen(line);
	int odd_allele = INVALID_ALLELE;


	// in case it didn't read an end-of-line character
	while (length > 0 && notNull != NULL)
	  {
	    // for each SNP
	    int index = 0;
	    while (index  < length) // index++ occurs in read_allele
	      {
		int allele1 = INVALID_ALLELE;
		if (odd_allele != INVALID_ALLELE){
		  allele1 = odd_allele;
		  odd_allele = INVALID_ALLELE;
		} else {
		  //printf("read_allele(len %i, idx %i) line[idx] = |%c|\n", length, index, line[index]);
		  allele1 = read_allele(&line[0], length, index);
		}
		if (index >= length && allele1 != OUTOFBOUNDS){
		  //printf("Odd allele: allele1 %i\n", allele1);
		  odd_allele = allele1;
		  break;
		}
		int allele2 = read_allele(&line[0], length, index);


		if (allele1 == DOUBLEDIGIT || allele2 == DOUBLEDIGIT){
		  printf(" INPUT ERROR: Double digit allele read on line %i\n", line_number);
		  exit(-1);
		}
		if (allele1 == INVALID_ALLELE || allele2 == INVALID_ALLELE){
		  printf(" INPUT ERROR: invalid allele (not %i or %i) read on line %i\n", FIRST_ALLELE, SECOND_ALLELE, line_number);
		  exit(-1);
		}
		if ((allele1 == OUTOFBOUNDS && allele2 != OUTOFBOUNDS)
		     || (allele2 == OUTOFBOUNDS && allele1 != OUTOFBOUNDS)){
		  printf(" FATAL ERROR: odd number of alleles in fgets() on line %i\n", line_number);
		  printf("   allele1 %i   allele2 %i\n", allele1, allele2);
		  printf("   at locus # %i\n", ((*person)->genotypes.size()+1) );
		  printf("   index %i, length %i\n", index, length);
		  printf("   odd_allele %i\n", odd_allele);
		  exit(-1);
		}

		if (allele1 != OUTOFBOUNDS && allele2 != OUTOFBOUNDS)
		  {
		    if (!isTyped && (allele1 != UNKNOWN || allele2 != UNKNOWN))
		      {
			isTyped = TRUE;
		      }


		    // We don't want to read the data, so do nothing
		  }
	      }


	    // read the next part of the buffer
	    if (line[length-1] != '\n' && notNull != NULL){
	      //printf("fgets\n");
	      notNull = fgets(line, 1024, fp);
	      length = strlen(line);
	    } else {
	      length = 0;
	    }
	  }
	//printf("  typed? %i\n", isTyped); fflush(stdout);



	// add individual to typed and untyped arrays
	if (isTyped){
	  (fam->typed).push_back( *person );
	} else {
	  (fam->untyped).push_back( *person );
	}
	(*person)->isTyped = isTyped;





	line_number++;

	if (feof(fp) != 0)
	  break;


	file_position = ftell(fp);	
      } // end while not EOF



  fclose(fp); // close in file 

}


// Temp fix - CE
void set_founders()
{

  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++) {
    fam->setIndices();  //


    
    // create founder kinship table in family
    while (fam->founder_kinship.size() < fam->founder.size()) {  // CE family -> fam
      fam->founder_kinship.push_back(vector<double>(fam->founder.size()));
    }
    
    
    
    // For each person in the family
    vector<Individual*>::iterator person;
    vector<Individual*>::iterator found;
    for (person = fam->members.begin();  person != fam->members.end();  person++) {
      //printf("PERSON: %i, f %i, m %i, s %i\n", person->individual_id, person->father_id, person->mother_id, person->sex);
      
      // mother
      if ((*person)->mother == NULL && (*person)->mother_id != UNKNOWN)	{
	int mother_id = (*person)->mother_id;
	found = fam->find_individual(fam->members.begin(), 
				     fam->members.end(), mother_id);
	if (found != fam->members.end())  {
	  (*person)->mother = (*found);
	}
      }
      
      // father
      if ((*person)->father == NULL && (*person)->father_id != UNKNOWN)
	{
	  int father_id = (*person)->father_id;
	  found = fam->find_individual(fam->members.begin(), 
				       fam->members.end(), father_id);
	  if (found != fam->members.end())
	    {
	      (*person)->father = (*found);
	    }
	}
      
      // children
	  for (found = fam->members.begin();  found != fam->members.end(); found++)
	    {
	      //printf("  FOUND: %i, f %i, m %i, s %i\n", (*person)->individual_id, (*person)->father_id, (*person)->mother_id, (*person)->sex);
	      if ((*found)->mother_id == (*person)->individual_id)
		{
		  if ((*person)->sex == FEMALE)
		    {
		      // set mother's child ptr and child's mother pointer
		      (*person)->children.push_back(*found);
		      if ((*found) != NULL){
			(*found)->mother = (*person);
		      } 
		    
		    } else {
		      printf("PED %i ERROR: Sex disagreement between non-female individual %i who is mother of %i\n", fam->family_id, (*person)->individual_id, (*found)->individual_id);
		    }
		} // end found's mother is person

	      if ((*found)->father_id == (*person)->individual_id)
		{
		  if ((*person)->sex == MALE)
		    {
		      // set father's child ptr and child's father pointer
		      (*person)->children.push_back(*found);
		      if ((*found) != NULL)
			{
			  (*found)->father = (*person);
			}
		    }
		  else
		    {
		      printf("PED %i ERROR: Sex disagreement between non-male individual %i who is father of %i\n", fam->family_id, (*person)->individual_id, (*found)->individual_id);
		    }
		} // end found's father is person
	    }
	}
    }


  
  
}


void read_founder_file(const char* in_file)
{

  // read the file
   FILE* fp = fopen (in_file, "r");
   if (fp == NULL)
     {
      printf("\nERROR: Cannot open %s.\n\n", in_file);
      exit(-1);
   }


   char column[256];

   //printf("Founder Kinship Coefficients:  ");
   while ( fscanf (fp, "%s ", column) != EOF)
     {
       int family_id = atoi(column);

       // Scan families for family pointer
       vector<Family>::iterator family = pedigrees.end();
       for (vector<Family>::iterator fam = pedigrees.begin();  
	    fam != pedigrees.end(); 
	    fam++) {
	 if (fam->family_id == family_id) {
	   family = fam;
	 }
       }

       if (family == pedigrees.end()){
	   printf("FATAL ERROR: (kinship file) family %i does not exist\n", family_id);
	   exit(-643);
	 }

       
       // create founder kinship table in family
       while (family->founder_kinship.size() < family->founder.size()) {
	 family->founder_kinship.push_back(vector<double>(family->founder.size()));
       }


       // read the first individual id
       if (fscanf(fp, "%s ", column) != EOF)
	 {
	   int indiv_id = atoi(column);
	   int founder_index = -1;

	   // scan the list of family founders for the individual pointer
	   bool found = false;
	   for (int i = 0;  i < family->founder.size();  i++)
	     {
	       if (family->founder[i]->individual_id == indiv_id)
		 {
		   founder_index = i;
		   found = true;
		 }
	     }
	   if (!found)
	     {
	       printf("FATAL ERROR: (kinship file) founder %i does not exist in family %i\n", indiv_id, family_id);
	       exit(-644);
	     }

	   
	   // read the second individual id
	   if (fscanf(fp, "%s ", column) != EOF)
	     {
	       int indiv_id2 = atoi(column);
	       int founder_index2 = -1;

	       // scan the list of family founders for the individual pointer
	       bool found = false;
	       for (int i = 0;  i < family->founder.size();  i++)
		 {
		   if (family->founder[i]->individual_id == indiv_id2)
		     {
		       founder_index2 = i;
		       found = true;
		     }
		 }
	       if (!found)
		 {
		   printf("FATAL ERROR: (kinship file) founder %i does not exist in family %i\n", indiv_id, family_id);
		   exit(-644);
		 }

	       
	       // read the kinship coefficient
	       if (fscanf(fp, "%s ", column) != EOF)
		 {
		   double kinship = atof(column);
	       
		   // Have two founder indices
		   family->founder_kinship[founder_index][founder_index2] = kinship;
		   family->founder_kinship[founder_index2][founder_index] = kinship;
		 }
	     }
	   else
	     {
	       printf("FATAL ERROR: kinship file must have (indiv id, indiv id, kinship) following every fam id.\n");
	       fflush(stdout);
	       exit(-45345);
	     }
	 }
     }
  fclose(fp);



  // print out the founder kinship
  printf("\nFOUNDER KINSHIP COEFFICIENTS:\n");
  for (vector<Family>::iterator fam = pedigrees.begin();  
       fam != pedigrees.end(); 
       fam++)
    {
      printf("  Family %i:\n", fam->family_id);
      for (int i=0;  i < fam->founder_kinship.size();  i++)
	{
	  for (int j=i;  j < fam->founder_kinship[i].size();  j++)
	    {
	      printf("    Founders %i,%i, %f\n", fam->founder[i]->individual_id, fam->founder[j]->individual_id, fam->founder_kinship[i][j]);
	    }
	}
    }
  printf("\n");   

}




double read_indiv_interest_file(const char* in_file)
{
  // read the graph  
   FILE* fp = fopen (in_file, "r");
   if (fp == NULL)
     {
      printf("\nERROR: Cannot open %s.\n\n", in_file);
      exit(-1);
   }


   char column[256];
   //printf("Recombination rates:  ");
   while ( fscanf (fp, "%s ", column) != EOF)
     {
       int family_id = atoi(column);

       // Scan families for family pointer
       vector<Family>::iterator family = pedigrees.end();
       for (vector<Family>::iterator fam = pedigrees.begin();  
	    fam != pedigrees.end(); 
	    fam++)
	 {
	   if (fam->family_id == family_id)
	     {
	       family = fam;
	     }
	 }

       if (family == pedigrees.end())
	 {
	   printf("FATAL ERROR: family of indiv of interest %i does not exist\n", family_id);
	   exit(-643);
	 }

       // read the individual id
       if (fscanf(fp, "%s ", column) != EOF)
	 {
	   int indiv_id = atoi(column);

	   // scan the list of family members for the individual pointer
	   // push the individual pointer into the indiv of interest list
	   bool found = false;
	   for (int i = 0;  i < family->members.size();  i++)
	     {
	       if (family->members[i]->individual_id == indiv_id)
		 {
		   family->indiv_interest.push_back(family->members[i]);
		   found = true;
		 }
	     }
	   if (!found)
	     {
	       printf("FATAL ERROR: indiv of interest %i does not exist in family %i\n", indiv_id, family_id);
	       exit(-644);
	     }
	 }
       else
	 {
	   printf("FATAL ERROR: indiv. of interest file must have indiv id following every fam id.\n");
	   fflush(stdout);
	   exit(-45345);
	 }
     }
   //printf("\n");

  fclose(fp);

  // print out the individuals of interest
  printf("\nINDIVIDUALS OF INTEREST:\n");
  for (vector<Family>::iterator fam = pedigrees.begin();  
       fam != pedigrees.end(); 
       fam++)
    {
      printf("  Family %i:\n", fam->family_id);
      for (int i=0;  i < fam->indiv_interest.size();  i++)
	{
	  printf("    Indiv %i\n", fam->indiv_interest[i]->individual_id);
	}
    }
  printf("\n");
}




// Sets pointers between individuals in the same pedigree
//  (mother, father, and children)
void set_pedigree_pointers()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      fam->setIndices();

      // For each person in the family
      vector<Individual*>::iterator person;
      vector<Individual*>::iterator found;
      for (person = fam->members.begin();  person != fam->members.end();  person++)
	{
	  //printf("PERSON: %i, f %i, m %i, s %i\n", person->individual_id, person->father_id, person->mother_id, person->sex);

	  // mother
	  if ((*person)->mother == NULL && (*person)->mother_id != UNKNOWN)
	    {
	      int mother_id = (*person)->mother_id;
	      found = fam->find_individual(fam->members.begin(), 
					   fam->members.end(), mother_id);
	      if (found != fam->members.end())
		{
		  (*person)->mother = (*found);
		}
	    }

	  // father
	  if ((*person)->father == NULL && (*person)->father_id != UNKNOWN)
	    {
	      int father_id = (*person)->father_id;
	      found = fam->find_individual(fam->members.begin(), 
					   fam->members.end(), father_id);
	      if (found != fam->members.end())
		{
		  (*person)->father = (*found);
		}
	    }

	  // children
	  for (found = fam->members.begin();  found != fam->members.end(); found++)
	    {
	      //printf("  FOUND: %i, f %i, m %i, s %i\n", (*person)->individual_id, (*person)->father_id, (*person)->mother_id, (*person)->sex);
	      if ((*found)->mother_id == (*person)->individual_id)
		{
		  if ((*person)->sex == FEMALE)
		    {
		      // set mother's child ptr and child's mother pointer
		      (*person)->children.push_back(*found);
		      if ((*found) != NULL){
			(*found)->mother = (*person);
		      } 
		    
		    } else {
		      printf("PED %i ERROR: Sex disagreement between non-female individual %i who is mother of %i\n", fam->family_id, (*person)->individual_id, (*found)->individual_id);
		    }
		} // end found's mother is person

	      if ((*found)->father_id == (*person)->individual_id)
		{
		  if ((*person)->sex == MALE)
		    {
		      // set father's child ptr and child's father pointer
		      (*person)->children.push_back(*found);
		      if ((*found) != NULL)
			{
			  (*found)->father = (*person);
			}
		    }
		  else
		    {
		      printf("PED %i ERROR: Sex disagreement between non-male individual %i who is father of %i\n", fam->family_id, (*person)->individual_id, (*found)->individual_id);
		    }
		} // end found's father is person
	    }
	}
    }


}




void allocate_markers(int number_of_snps)
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      //printf("FAMILY %i\n", fam->family_id);

      // For each person in the family
      vector<Individual*>::iterator person;
      for (person = fam->members.begin();  person != fam->members.end();  person++)
	{
	  // Allocate memory for all the markers we will simulate
	  for (int marker = 0;  marker < number_of_snps;  marker++)
	    {
	      Marker snp(UNKNOWN, UNKNOWN);
	      (*person)->genotypes.push_back(snp);
	      (*person)->identity.push_back(snp);
	    }
	}
    }
}












// returns a particular ordering on the members below the founder in the pedigree
// While updating the order vector, the individuals are stored in the same order that they are
// in the family's member list.
void pedigree_dfs(vector<Family>::iterator fam, 
		  Individual* parent, 
		  int depth)
{
  vector<Individual*>::iterator c;
  for (c = parent->children.begin();  c != parent->children.end();  c++)
    {
      Individual* child = (*c);
      int child_index = child->individual_index;

      // check that fam->order has defined key for this person
      if (fam->order[child_index].individual_id == UNKNOWN)
	{
	  fam->order[child_index].individual_id = child->individual_id;
	  fam->order[child_index].individual_index = child->individual_index;
	  fam->order[child_index].individual = child;
	}
      if (depth > fam->order[child_index].depth)
	{
	  fam->order[child_index].depth = depth;
	  //printf("    child id %i, depth %i\n", child->individual_id, depth);
	  //fflush(stdout);
	}

      pedigree_dfs(fam, child, depth+1);
    }
}


void set_traversal_orders()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      //printf("FAMILY %i\n", fam->family_id);

      // Allocate the order vector
      for (unsigned int c = 0;  c < fam->members.size();  c++)
	{
	  fam->order.push_back(OrderItem());
	}


      // For each founder in the family, determine the top-down order
      // a depth-first search on the pedigree DAG
      for (unsigned int f = 0;  f < fam->founder.size();  f++)
	{
	  Individual* founder = fam->founder[f];

	  fam->order[founder->individual_index].individual_index = founder->individual_index;
	  fam->order[founder->individual_index].individual_id = founder->individual_id;
	  fam->order[founder->individual_index].individual = founder;
	  fam->order[founder->individual_index].depth = 0;
	  pedigree_dfs(fam, founder, 1);
	}

      // sort the order vector by decreasing depth
      sort((fam->order).begin(), (fam->order).end(), order_lt());


      //vector<OrderItem>::iterator bottom_up;
      //for (bottom_up = fam->order.begin();  bottom_up != fam->order.end();  bottom_up++)
      //{
      //  printf("      Person %i:%i depth %i\n", fam->family_id, bottom_up->individual_id, bottom_up->depth);
      //  fflush(stdout);
      //}
    }
}






//
// top-down dynamic programming recursion to compute kinship
//
// founders phi_ff is initialized with kinship coefficients
// every other phi_fj = 0
// phi_ij = (phi_mj + phi_pj)/2 where i is not ancestor of j and i != j
//
int compute_kinship()
{

  // print out all genotypes (ordered by parental source)
  /*  FILE* fp = fopen (out_kinship_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", out_kinship_file);
      exit(-1);
    }
  */

  // For each family
  int hap_index = 0;
  int k = 0; // founder counter to read inbreeding coeffients

  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      vector<vector<double> > phi(fam->members.size(), vector<double>(fam->members.size(), -1.0));

      // for individuals get the top-down order
      vector<Individual*> top_down;
      vector<Individual*> bottom_up;
      stack<Individual*> top_down_queue;
      vector<bool> visited = vector<bool>(fam->members.size(), 0);

      // Find all founders
      for (int f = 0;  f < fam->members.size();  f++)
	{
	  if (fam->members[f]->father_id == UNKNOWN && fam->members[f]->mother_id == UNKNOWN)
	    {
	      top_down_queue.push(fam->members[f]);
	    }
	}
      // populate queue in top-down order
      while(!top_down_queue.empty())
	{
	  Individual* i = top_down_queue.top();

	  // BK added the third condition, 4/23/17
	  bool all_visited = true;
          if (i->children.size() > 0 && visited[i->individual_index] == 0)
	    {
	      for (int c = 0;  c < i->children.size();  c++)
		{
		  Individual* child = i->children[c];
      		  if (visited[child->individual_index] == 0)
      		    {
      		      top_down_queue.push(child);
		      all_visited = false;
      		    }
		}
	    }
	  if (all_visited && visited[i->individual_index] == 0)
	    {
	      visited[i->individual_index] = 1;
	      if (i != top_down_queue.top())
		{
		  printf("ERROR in topo sort\n");
		}
	      top_down_queue.pop();
	      bottom_up.push_back(i);
	    }
	}

      // reverse the list
      for (int c=bottom_up.size()-1;  c >= 0;  c--)
	{
	  top_down.push_back(bottom_up[c]);
	}

      // printf("\n");
      // printf("top down size: %i\n", top_down.size());
      // for (int c=0;  c < top_down.size();  c++)
      // 	{
      // 	  printf("%i] person %i\n", c, top_down[c]->individual_id);
      // 	}
      // printf("\n");
      

      Rcout << "kjhkjhk" << std::endl;

      // Initialize the ancestor sets by top-down recursion
      vector<vector<int> > ancestors(fam->members.size(), vector<int>(fam->members.size(), false));
      for (int f = 0;  f < fam->founder.size();  f++)
	{
	  int i = fam->founder[f]->individual_index;
	  ancestors[i][i] = true;
	}
      for (int x = 0;  x < top_down.size();  x++)
	{
	  int i = top_down[x]->individual_index;
	  ancestors[i][i] = true;
	  for (int y = 0;  y < top_down.size();  y++)
	    {
	      int j = top_down[y]->individual_index;
	      if (top_down[x]->mother_id != UNKNOWN)
		{
		  if (ancestors[top_down[x]->mother->individual_index][j] == true)
		    {
		      ancestors[i][j] = true;
		    }
		}
	      if (top_down[x]->father_id != UNKNOWN)
		{
		  if (ancestors[top_down[x]->father->individual_index][j] == true)
		    {
		      ancestors[i][j] = true;
		    }
		}
	    }
	}


            Rcout << "werwrerwer2" << std::endl;

      // initialize all the founder kinship coefficients
      for (int f = 0;  f < fam->founder.size();  f++)
	{
	  for (int g = 0;  g < fam->founder.size();  g++)
	    {
	      int f_index = fam->founder[f]->individual_index;
	      int g_index = fam->founder[g]->individual_index;
	      if (f == g)
		{
		  phi[f_index][f_index] = (1+fam->founder_kinship[f][f])/2;
		}
	      else
		{
		  phi[f_index][g_index] = fam->founder_kinship[f][g];
		}

	    }
	}



            Rcout << "werwrerwer3" << std::endl;

	    
      // Compute the kinship by top-down recursion
      for (int x = 0;  x < top_down.size(); x++)
	{
	  Individual* j = top_down[x];
	  
	  // for individual (in top-down order)
	  vector<OrderItem>::iterator top_down2;
	  for (int y = x;  y < top_down.size(); y++) // changed from y=0 to y=x by Sophie on mar.22
	    {
	      Individual* i = top_down[y];

	      if (i->father_id == UNKNOWN && i->mother_id == UNKNOWN && i == j)
		{ // i,j are same founder
		}
	      if (i->father_id == UNKNOWN && i->mother_id == UNKNOWN && i != j)
		{
		  if (j->father_id == UNKNOWN && j->mother_id == UNKNOWN)
		    {
		    }
		  else
		    {
		      // i is founder
		      if (j->father->father_id == UNKNOWN || j->father->mother_id == UNKNOWN || j->mother->father_id == UNKNOWN || j->mother->mother_id == UNKNOWN)
			{
			  // j is child of founder
			  phi[j->individual_index][i->individual_index] = (phi[j->mother->individual_index][i->individual_index] + phi[j->father->individual_index][i->individual_index])/2;
			  phi[i->individual_index][j->individual_index] = phi[j->individual_index][i->individual_index];
			}
		      else
			{
			  if (ancestors[j->individual_index][i->individual_index] == false)
			    {
			      phi[i->individual_index][j->individual_index] = 0;
			    }
			}
		    }
		}
	      if (j->father_id == UNKNOWN && j->mother_id == UNKNOWN && i != j)
		{ 
		  if (i->father_id == UNKNOWN && i->mother_id == UNKNOWN)
		    {
		    }
		  else
		    {
		      // j is founder
		      if (i->father->father_id == UNKNOWN || i->father->mother_id == UNKNOWN || i->mother->father_id == UNKNOWN || i->mother->mother_id == UNKNOWN)
			{ // i is child of founder
			  phi[i->individual_index][j->individual_index] = (phi[i->mother->individual_index][j->individual_index] + phi[i->father->individual_index][j->individual_index])/2;
			  phi[j->individual_index][i->individual_index] = phi[i->individual_index][j->individual_index];
			}
		      else
			{
			  if (ancestors[i->individual_index][j->individual_index] == false)
			    {
			      phi[j->individual_index][i->individual_index] = 0;
			    }
			}
		    }
		}

	      if (phi[i->individual_index][j->individual_index] == -1.0)
	      {
		if (i == j)
		  {
		    phi[i->individual_index][j->individual_index] = (1+ phi[i->mother->individual_index][i->father->individual_index])/2;
		  }
		if (ancestors[j->individual_index][i->individual_index] == false)
		  {
		    phi[i->individual_index][j->individual_index] = (phi[i->mother->individual_index][j->individual_index] + phi[i->father->individual_index][j->individual_index])/2;
		    phi[j->individual_index][i->individual_index] = phi[i->individual_index][j->individual_index];
		  } 
	      } // end if not set
		  
	    } // end for each j
	} // end for each individual i



      // Transformation to phi[i][i] to convert kinship to 
      // inbreeding coefficient
      for (int i=0;  i < phi.size(); i++)
      {
        phi[i][i] = 2*phi[i][i]-1;
      }


      Rcout << "werwrerwer" << std::endl;

      // print the kinship information to the file
      for (int i=0; i < fam->indiv_interest.size(); i++)
	{
	  for (int j=i;  j < fam->indiv_interest.size(); j++)
	    {
	      int x = fam->indiv_interest[i]->individual_index;
	      int y = fam->indiv_interest[j]->individual_index;
	      //	      fprintf(fp, "%i %i %i %f\n", fam->family_id, fam->indiv_interest[i]->individual_id, fam->indiv_interest[j]->individual_id, phi[x][y]);
	      Rcout << fam->family_id << " " << fam->indiv_interest[i]->individual_id << " " << fam->indiv_interest[j]->individual_id << " " << phi[x][y] << std::endl;
	    }
	}

    } // end for each family



  //  fclose(fp); // close out file 

      return(0);
}




//
// dynamic programming recursion to 
// sample identity states and estimate kinship
//
// founder CC membership is adjusted randomly according to 
// founder kinship coefficents
//
int sample_kinship(const char* out_kinship_file, const int S)
{

  // print out all genotypes (ordered by parental source)
  FILE* fp = fopen (out_kinship_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", out_kinship_file);
      exit(-1);
    }


  // For each family
  int hap_index = 0;
  int k = 0; // founder counter to read inbreeding coeffients

  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      vector<vector<double> > phi(fam->members.size(), vector<double>(fam->members.size(), 0.0));

      // For the number of samples S
      for (int s = 0;  s < S; s++)
	{

	  // choose an inheritance path by setting segregation indicators
	  // initialize the cc_membership as zeros
	  for(int i = 0;  i < fam->members.size();  i++)
	    {
	      while (fam->members[i]->segregation_indicator.size() < 2)
		{
		  fam->members[i]->segregation_indicator.push_back(false);
		}
	      while (fam->members[i]->cc_membership.size() < 2)
		{
		  fam->members[i]->cc_membership.push_back(0);
		}
	      for (int c = 0; c < 2; c++)
		{
		  double u = drand48();
		  if (u < 0.5)
		    {
		      fam->members[i]->segregation_indicator[c]=true;
		    }
		  else
		    {
		      fam->members[i]->segregation_indicator[c]=false;
		    }
		  fam->members[i]->cc_membership[c] = 0;
		}
	    }

	  // set cc_membership
	  int counter = 1;
	  vector<OrderItem>::iterator bottom_up;
	  for (bottom_up = fam->order.begin();  bottom_up != fam->order.end();  bottom_up++)
	    {
	      Individual* indiv = bottom_up->individual;

	      // leaf individuals
	      if (indiv->children.size() == 0)
		{
		  for (int c = 0; c < 2; c++)
		    { 
		      indiv->cc_membership[c] = counter;
		      counter++;
		    }
		}
	      // copy the cc_membership up the inheritance path
	      for (int c = 0; c < 2; c++)
		{ 
		  Individual* parent = NULL;
		  if (c == 0 && indiv->father_id != UNKNOWN)
		    {
		      parent = indiv->father;
		    }
		  if (c == 1 && indiv->mother_id != UNKNOWN)
		    {
		      parent = indiv->mother;
		    }

		  if (parent != NULL)
		    {
		      int s = indiv->segregation_indicator[c];
		      if (parent->cc_membership[s] < indiv->cc_membership[c])
			{
			  parent->cc_membership[s] = indiv->cc_membership[c];
			}	
		    }	      
		}
	    } // end of bottom-up order: cc_membership partially set



	  // draw the founder cc_membership that is inbred
	  for (int f = 0;  f < fam->founder.size();  f++)
	    {
	      // single founder inbreeding
	      double u = drand48();
	      if (u < fam->founder_kinship[f][f])
		{
		  //printf("inbred founder from prob %f\n", fam->founder_kinship[f][f]);
		  fam->founder[f]->cc_membership[0] = fam->founder[f]->cc_membership[1];
		}

	      for (int g = f+1;  g < fam->founder.size();  g++)
		{
		  // v == choice of which outbred pair of edges to consider
		  // (a1,b1)(a2,b2)  or  (a1,b2)(a2,b1)
		  // after choosing which pair use kinship prob to select
		  // whether the edges exist
		  double v = drand48();
		  for (int c = 0; c < 2; c++)
		    {
		      if (f != g)
			{
			  double u = drand48();
			  if (u < 2*fam->founder_kinship[f][g])
			    {
			      //printf("inbred founder pair (%i,%i) from prob %f\n", fam->founder[f]->individual_id, fam->founder[g]->individual_id, fam->founder_kinship[f][g]);
			      // copy this allele from a random founder allele
			      if (v < 0.5)
				{
				  fam->founder[f]->cc_membership[c] = fam->founder[g]->cc_membership[1-c];
				}
			      else
				{
				  fam->founder[f]->cc_membership[c] = fam->founder[g]->cc_membership[c];
				}
			    }
			}
		    }
		}
	    }


	  // for non-founder individual (in top-down order)
	  // propagate the cc_membership down the inheritance path
	  vector<OrderItem>::iterator top_down;
	  for (top_down = fam->order.end()-1;  top_down != fam->order.begin()-1;  top_down--)
	    {
	      Individual* indiv = top_down->individual;

	      // check if not founder
	      if (indiv->father_id != UNKNOWN && indiv->mother_id != UNKNOWN)
		{
		  for (int c = 0; c < 2; c++)
		    { 
		      Individual* parent = NULL;
		      if (c == 0 && indiv->father_id != UNKNOWN)
			{
			  parent = indiv->father;
			}
		      if (c == 1 && indiv->mother_id != UNKNOWN)
			{
			  parent = indiv->mother;
			}
		  
		      if (parent != NULL)
			{
			  int s = indiv->segregation_indicator[c];
			  indiv->cc_membership[c] = parent->cc_membership[s];
			}

		    }  
		}
	    }



	  // Create the identtity state graph, count the edges,
	  // and make the updates to the kinship, phi
	  for (int i=0; i < fam->indiv_interest.size(); i++)
	    {
	      Individual* indiv_i = fam->indiv_interest[i];
	      int i_idx = indiv_i->individual_index;


	      if (indiv_i->cc_membership[0] == indiv_i->cc_membership[1])
		{
		  //printf("inbreeding indiv %i cc=(%i,%i)\n", indiv_i->individual_id, indiv_i->cc_membership[0],indiv_i->cc_membership[1]);
		  phi[i_idx][i_idx] += 1.0/(double)S;
		}

	      for (int j=i+1;  j < fam->indiv_interest.size(); j++)
		{
		  Individual* indiv_j = fam->indiv_interest[j];
		  int j_idx = indiv_j->individual_index;

		  // count the number of edges
		  int edge_count = 0;
		  for (int c = 0;  c < 2;  c++)
		    {
		    for (int d = 0;  d < 2; d++)
		      {
			if (indiv_i->cc_membership[c] == indiv_j->cc_membership[d])
			  {
			    edge_count++;
			  }
		      }
		    }

		  // update phi
		  phi[i_idx][j_idx] += (double)edge_count/(4.0*(double)S);
		  //if (edge_count > 0){printf("IBD %i, phi %f, term %f\n",edge_count,phi[i_idx][j_idx], (double)edge_count/(4.0*(double)S));}
		  phi[j_idx][i_idx] += (double)edge_count/(4.0*(double)S);
		} // end for indiv intrest j
	    }// end for indiv interest i


	} // end for each sample





      // print the kinship information to the file
      for (int i=0; i < fam->indiv_interest.size(); i++)
	{
	  for (int j=i;  j < fam->indiv_interest.size(); j++)
	    {
	      int x = fam->indiv_interest[i]->individual_index;
	      int y = fam->indiv_interest[j]->individual_index;
	      fprintf(fp, "%i %i %i %f\n", fam->family_id, fam->indiv_interest[i]->individual_id, fam->indiv_interest[j]->individual_id, phi[x][y]);
	    }
	}


    } // end for each family

}





void print_genotypes(const char* hap_file)
{
  // print out all genotypes (ordered by parental source)
  FILE* fp = fopen (hap_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", hap_file);
      exit(-1);
    }


  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  fprintf(fp, "%i %i ", person->family_id, person->individual_id);
	  fprintf(fp, "%i %i ", person->father_id, person->mother_id);
	  fprintf(fp, "%i %i ", person->sex, person->affection);


	  for (unsigned int snp = 0; snp < person->genotypes.size(); snp++)
	    {
	      long allele1 = person->genotypes[snp].allele1;
	      long allele2 = person->genotypes[snp].allele2;

	      fprintf(fp," %i %i ", allele1, allele2);	      
	    }
	  fprintf(fp, "\n");
	}
    }

  fclose(fp);
}


void print_ibd(const char* hap_file)
{
  // print out all genotypes (ordered by parental source)
  FILE* fp = fopen (hap_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", hap_file);
      exit(-1);
    }


  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  fprintf(fp, "%i %i ", person->family_id, person->individual_id);
	  fprintf(fp, "%i %i ", person->father_id, person->mother_id);
	  fprintf(fp, "%i %i ", person->sex, person->affection);


	  for (unsigned int snp = 0; snp < person->genotypes.size(); snp++)
	    {
	      long allele1 = person->identity[snp].allele1;
	      long allele2 = person->identity[snp].allele2;

	      fprintf(fp," %i %i ", allele1, allele2);	      
	    }
	  fprintf(fp, "\n");
	}
    }

  fclose(fp);
}







void deleteped()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      // For each person in the family
      vector<Individual*>::iterator person;
      for (person = fam->members.begin();  person != fam->members.end();  person++)
	{
	  delete (*person);
	}
    }  
}




//////////////////////////////
//
//
int yegoodeoldmain (int argc, char *argv[])
{
  // read arguements
  int num_args = 5; // for exact
  if (argc <= num_args)  // argc is # argvs [1..n] with prog first
    {
      printf("\n");
      printf("  Program: PedKin version %s \n", version);
      printf("   Author: Bonnie Kirkpatrick \n");
      printf("  Company: Intrepid Net Computing \n");
      printf("     Date: Feb. 12, 2016 \n");
      printf("  License: GNU GPL v2.0 \n");
      printf(" \n");
      printf(" \n");
      printf("  Given a pedigree of individuals, founder kinship coefficients,  \n");
      printf("  and the individuals of interest compute the kinship on all \n");
      printf("  pairs of individuals individuals of interest in the pedigree. \n");
      printf("  Uses the founder kinship coefficients to initialize the kinship. \n");
      printf(" \n");
      printf("  Outputs kinship matrix with the diagonal transformed into the  \n");
      printf("  inbreeding coefficients. \n");
      printf(" \n");
      printf("  Algorithm properties: \n");
      printf("    Exact Algorithm \n");
      printf("      O(n^2) time for n individuals in a pedigree \n");
      printf("      O(n^2) space \n");
      printf("    Approximate Algorithm \n");
      printf("      O(n) time for n individuals and sqrt(n) individuals of interest \n");
      printf("      O(n) space \n");
      printf("\n\n");


      printf ("\nUSAGE: %s [algorithm] (iterations) [ped file] [founder kinship] [indiv interest] [kinship output]\n", argv[0]);
      printf ("[algorithm] choose [e|a] for exact or approximate\n");
      printf ("(iterations) if algorithm a, then give integer for iterations\n");
      printf ("[ped file]  ped file containing pedigree structure\n");
      printf ("[founder kinship]  list of (fam, indiv, indiv, kinship) separated by \'\\n\' \n");
      printf ("[indiv interest]  list of FamID and IndivID separated by \'\\n\' \n");
      printf ("[kinship output]  list of (fam, indiv, indiv, kinship) values separarated by \'\\n\'\n\n");
      printf ("\n");

      exit(-1);
    }

  printf("\n=======================\n");
  printf("  PedKin version %s\n", version);
  printf("=======================\n");

  const char* algorithm = argv[1];

  int param = 2;
  int iterations = 0;
  if (algorithm[0] == 'e')
    {
      printf("\nUsing EXACT algorithm\n\n");
    }
  else if (algorithm[0] == 'a')
    {
      const char* it = argv[param];
      param++;
      iterations = atoi(it);
      printf("\nUsing APPROXIMATE algorithm with %i iterations\n\n", iterations);
      if (argc <= num_args+1) // num args for approximate
	{
	  printf("FATAL ERROR: approximate algorithm needs number of iterations in parameter list\n");
	  printf("             run \n");
	  printf("               > ./%s \n", argv[0]);
	  printf("             for USAGE\n");
	  exit(-5345);
	}
    }
  else
    {
      printf("\nFATAL ERROR: invalid algorithm selection, choose [e|a]\n");
      exit(-254367);
    }
  printf("\n");


  const char* ped_file = argv[param];
  param++;
  const char* founder_file = argv[param];
  param++;
  const char* indiv_interest_file = argv[param];
  param++;
  const char* out_kinship_file = argv[param];



  printf("PARAMS: algorithm = %s, iterations = %i \\\n", algorithm, iterations);
  printf("        ped_file = %s, founder_file = %s, \\\n", ped_file, founder_file);
  printf("        indiv_interest_file = %s, out_kinship_file = %s\n", indiv_interest_file, out_kinship_file);
  printf("\n");





  srand48(time(NULL));

  printf("READING input pedigree...\n");

  // FIX HERE
  read_pedigree_file(ped_file);
  set_pedigree_pointers();

  set_traversal_orders();



  printf("READING the founder kinship...\n");
  read_founder_file(founder_file);



  printf("READING the individuals of interest...\n");
  read_indiv_interest_file(indiv_interest_file);


  if (algorithm[0] == 'e')  {
    //    compute_kinship(out_kinship_file);
  }
  else if (algorithm[0] == 'a') {
    sample_kinship(out_kinship_file, iterations);
  }
  


  deleteped();
}



// chrtype
// approximation
//' @export 
// [[Rcpp::export]]
void create_kinship(const IntegerVector& famid, const IntegerVector& id, const IntegerVector& fid, const IntegerVector& mid, const IntegerVector& sex) {

  // Sanity checks

  // First we create the pedigrees based on the information
  long N=famid.size();

  for (int i=0; i<N; i++) {
    // Find the family or make it
    Family newfam(famid[i]);
    vector<Family>::iterator fam;
    fam = find(pedigrees.begin(), pedigrees.end(), newfam);
    if (fam == pedigrees.end()){
      // make family
      pedigrees.push_back(newfam);
      fam = pedigrees.end() - 1;
    }

    // insert individual into pedigree
    // aff
    int isFounder = (fid[i]==0 && mid[i]==0);
    
    Individual* newperson = new Individual(famid[i], id[i], fid[i], mid[i], sex[i], sex[i], 1, i);
    vector<Individual*>::iterator person;
    (fam->members).push_back(newperson);
    person = fam->members.end() - 1;
    if (isFounder){
      (fam->founder).push_back(newperson);
    } else {
      (fam->nonfounder).push_back(newperson);
    }
    
    // Read all the alleles
    int isTyped = FALSE;

    // add individual to typed and untyped arrays
    if (isTyped){
      (fam->typed).push_back( *person );
    } else {
      (fam->untyped).push_back( *person );
    }
    (*person)->isTyped = isTyped;
    
  }

  
  set_pedigree_pointers();
  
  set_traversal_orders();


  set_founders();

  Rcout << "aloha" << std::endl;

  // Now set the founder information.
  // For now assume that founders are non-inbred

  //

  // Iterate over all families then all founders

  /*
  for (vector<Family>::iterator mainfam = pedigrees.begin();  
       mainfam != pedigrees.end(); 
       mainfam++) {

    int family_id = mainfam->family_id;

    // Scan families for family pointer
    vector<Family>::iterator family = pedigrees.end();
    
    for (vector<Family>::iterator fam = pedigrees.begin();  
	 fam != pedigrees.end(); 
	 fam++) {
      if (fam->family_id == family_id) {
	family = fam;
      }
    }
    
    if (family == pedigrees.end()) {
      stop("coding error (111)");
      //	printf("FATAL ERROR: (kinship file) family %i does not exist\n", family_id);
      //	exit(-643);
    }
    
    
    // create founder kinship table in family
    while (family->founder_kinship.size() < family->founder.size()) {
      family->founder_kinship.push_back(vector<double>(family->founder.size()));
    }
    

    // Now run through the founders and code them with NO inbreeding
    for (int ii=0;  ii < mainfam->founder_kinship.size();  ii++)  {

      int indiv_id = mainfam->founder[ii]->individual_id;
      int founder_index = -1;
      
      // scan the list of family founders for the individual pointer
      bool found = false;
      for (int i = 0;  i < family->founder.size();  i++) {
	  if (family->founder[i]->individual_id == indiv_id) {
	    founder_index = i;
	    found = true;
	  }
      }
      if (!found) {
	stop("Found an error here (222)");
	  //	printf("FATAL ERROR: (kinship file) founder %i does not exist in family %i\n", indiv_id, family_id);
	  //	exit(-644);
      }
      
      
      for (int j=i;  j < fam->founder_kinship[ii].size();  j++) {
	printf("    Founders %i,%i, %f\n", fam->founder[i]->individual_id, fam->founder[j]->individual_id, fam->founder_kinship[i][j]);
      }

    
    








    

    
    printf("  Family %i:\n", mainfam->family_id);
  }
  

  //printf("Founder Kinship Coefficients:  ");



  
       // read the first individual id
       if (fscanf(fp, "%s ", column) != EOF)
	 {
	   int indiv_id = atoi(column);
	   int founder_index = -1;

	   // scan the list of family founders for the individual pointer
	   bool found = false;
	   for (int i = 0;  i < family->founder.size();  i++)
	     {
	       if (family->founder[i]->individual_id == indiv_id)
		 {
		   founder_index = i;
		   found = true;
		 }
	     }
	   if (!found)
	     {
	       printf("FATAL ERROR: (kinship file) founder %i does not exist in family %i\n", indiv_id, family_id);
	       exit(-644);
	     }

	   
	   // read the second individual id
	   if (fscanf(fp, "%s ", column) != EOF)
	     {
	       int indiv_id2 = atoi(column);
	       int founder_index2 = -1;

	       // scan the list of family founders for the individual pointer
	       bool found = false;
	       for (int i = 0;  i < family->founder.size();  i++)
		 {
		   if (family->founder[i]->individual_id == indiv_id2)
		     {
		       founder_index2 = i;
		       found = true;
		     }
		 }
	       if (!found)
		 {
		   printf("FATAL ERROR: (kinship file) founder %i does not exist in family %i\n", indiv_id, family_id);
		   exit(-644);
		 }

	       
	       // read the kinship coefficient
	       if (fscanf(fp, "%s ", column) != EOF)
		 {
		   double kinship = atof(column);
	       
		   // Have two founder indices
		   family->founder_kinship[founder_index][founder_index2] = kinship;
		   family->founder_kinship[founder_index2][founder_index] = kinship;
		 }
	     }
	   else
	     {
	       printf("FATAL ERROR: kinship file must have (indiv id, indiv id, kinship) following every fam id.\n");
	       fflush(stdout);
	       exit(-45345);
	     }
	 }
     }
  fclose(fp);

  */

  // Will maybe need to set founders values with themselves

  























  compute_kinship();
  

  
  /*
  if (algorithm[0] == 'e')  {
    compute_kinship(out_kinship_file);
  }
  else if (algorithm[0] == 'a') {
    sample_kinship(out_kinship_file, iterations);
  }
  */



  // Garbage collection
  deleteped();
  //  return(0);
    
}


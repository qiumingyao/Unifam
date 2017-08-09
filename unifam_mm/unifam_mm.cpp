#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <string>

#include "utils.hpp"

#define WORKTAG    1
#define DIETAG     2

using namespace std;
// slave process, everybody is slave...
void MasterProcess(vector<string>& hmm_protein_pairs, const string & meta_info_file_hmm,const string & meta_info_file_protein);
void SlaveProcess(const string & file_dir_hmm,const string & file_dir_proteins, const string & out_directory, const string & numcpu);

vector<string> explode(const string& str, const char& ch) {
    string next;
    vector<string> result;

    // For each character in the string
    for (string::const_iterator it = str.begin(); it != str.end(); it++) {
        // If we've hit the terminal character
        if (*it == ch) {
            // If we have some characters accumulated
            if (!next.empty()) {
                // Add them to the result vector
                result.push_back(next);
                next.clear();
            }
        } else {
            // Accumulate the next character into the sequence
            next += *it;
        }
    }
    if (!next.empty())
         result.push_back(next);
    return result;
}

/* print usage */
void usage()
{
    cout << " [Usage]" << endl
    << "  unifam_mm [options] -cpu <cores> -ih <file_dir_hmm> -ip <file_dir_proteins> -fh <meta_info_file_hmm> -fp <meta_info_file_protein> -od <out_directory>" << endl

    << endl
	<< " This is MPI version for unifam to do the hmmsearch distributedly, "<<endl
	<< "  hmm models and proteins need to be both splitted accordingly, and have size information stored in two meta files respectively" << endl
    << " [Inputs]" << endl
	<< " Please don't put '\\' or '/' after the name of the directory" << endl
	<< " -ih <file_dir_hmm> file directory for hidden markov models used for protein families" << endl
	<< " -ip <file_dir_proteins> file directory for query protein chunks" << endl
	<< " -fh <meta_info_file_hmm> file with meta information of hidden markov model chunks, file_name (tab) max_len (tab) min_len" << endl
	<< " -fh <meta_info_file_hmm> file with meta information of protein chunk files, file_name (tab) max_len (tab) min_len" << endl
    << " -cpu <cores> mutlithreading numbers for hmmsearch"<<endl
	<< endl
    << " [Outputs]" << endl
    << " -od <out_directory> output direcotry for hmmer search results" << endl
    << endl
    << " [Options]" << endl
    << "  -h/--help Display this help message" << endl
    << endl;
}

bool initializeArguments(int argc, char **argv,
                         string & file_dir_hmm,
						 string & file_dir_proteins,
						 string & meta_info_file_hmm,
						 string & meta_info_file_protein,
						 string & out_directory,
						 string & numcpu
)
{
    vector<string> Arguments;
    file_dir_hmm = "";
    file_dir_proteins = "";
    meta_info_file_hmm = "";
    meta_info_file_protein = "";
    out_directory = "";
    numcpu="1";
    int paranum = 0;



    while(argc--)
        Arguments.push_back(*argv++);

    for(int i = 1; i <= (int)Arguments.size()-1; i++)
    {

        if (Arguments[i] == "-ih")
        {
        	file_dir_hmm = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-ip")
        {
        	file_dir_proteins = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-fh")
        {
        	meta_info_file_hmm = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-fp")
        {
        	meta_info_file_protein = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-od")
        {
        	out_directory = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-cpu")
        {
        	numcpu = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-h")
        {
            usage();
            return false;
        }
        else
        {
        	cout << "Unknown option " << Arguments[i] << endl << endl;
            usage();
            return false;
        }
    }
    // check required arguments
    if (paranum <6)
    {
    	cout << "miss necessary parameter(s)"<<endl;
    	cout << "use -h/--help for help" << endl;
        return false;
    }

    return true;

}



// main function

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  main function
 * =====================================================================================
 */
int main(int argc, char ** argv){

    /********************** variable declaration *********************/
	bool successflag=true; //whole flag control for error or exceptions

	string numcpu,file_dir_hmm,file_dir_proteins,meta_info_file_hmm,meta_info_file_protein,out_directory;


    vector<string> hmm_protein_pairs;

    double start_time, finish_time;

    /********************** initialize MPI  *********************/
    int myid, num_of_nodes;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_nodes);

    /********************** display work start and time record *********************/
    if (myid == 0) {
        start_time = MPI_Wtime();

        cout << endl
        << "============================================================================"
        << endl << Utils::currentDateTime() << endl
        << " Initialize commandline arguments -> " << endl;
    }

    /************ initialize arguments, for input, output, and options ***********/

    successflag = initializeArguments(argc, argv, file_dir_hmm,file_dir_proteins,meta_info_file_hmm,meta_info_file_protein,out_directory,numcpu);
    if(successflag == false)
    	{
    	cout << "[ERROR]::Initializing commands improperly!" << endl;
    	return 0;
    	}

    Utils::mkdirIfNonExist(out_directory);
//    string tempdir = out_directory+Utils::getPathSeparator()+"temp";
//    Utils::mkdirIfNonExist(tempdir);

    /* echo the arguments */
    if (myid == 0) {

        cout<< " The arguments are taken!" << endl << endl;
    }

    /********************** generate data for each node, and find chunk information *********************/

    if(num_of_nodes > 1){
	    if (myid == 0) {
		    MasterProcess(hmm_protein_pairs, meta_info_file_hmm,meta_info_file_protein);
		    cout << "  Running with total " << num_of_nodes << " processes" << endl;
	    }

	    else
		    SlaveProcess(file_dir_hmm,file_dir_proteins, out_directory, numcpu);
    }
    else{ /* If there is only one node, running in multi-thread mode */
	    cout << "[ERROR]::Only 1 process, reject to run. "<<endl;
	    return 0;
    }

    /********************** display work end and time record *********************/

    if( myid == 0 )
    {
        cout << " Done!" << endl;
        finish_time = MPI_Wtime();

        // display work end and time record
        cout << Utils::currentDateTime() << " all nodes done "<< endl
        << "Total Elapsed Time =  "
        << double(finish_time - start_time) << " [seconds]" << endl
        << "============================================================================"
        << std::endl << std::endl;
    }

     /*************************** Finalize MPI  *******************************/
    MPI_Finalize();
    return 1;
}

void MasterProcess(vector<string>& hmm_protein_pairs, const string & meta_info_file_hmm,const string & meta_info_file_protein)
{
	//create the hmm and protein file pairs
	ifstream myHMMFile, myProtFile;
	vector<string> hmmfiles;
	myHMMFile.open (meta_info_file_hmm.c_str());
	while(!myHMMFile.eof())
	{
		string text;
		getline(myHMMFile,text);
		if(text!=""&&text.find("file_name")!=0)
		 hmmfiles.push_back(text);
	}
	myHMMFile.close();

	myProtFile.open (meta_info_file_protein.c_str());
	while(!myProtFile.eof())
	{
		string text2;
		getline(myProtFile,text2);
		if(text2!=""&&text2.find("file_name")!=0)
		{
			for(unsigned int i = 0; i< hmmfiles.size(); i++)
			{
				string newpair = hmmfiles.at(i)+"\t"+text2;
				hmm_protein_pairs.push_back(newpair);
			}
		}
	}
	myProtFile.close();
	hmmfiles.clear();



	std::random_shuffle ( hmm_protein_pairs.begin(), hmm_protein_pairs.end());
	int num_of_jobs = hmm_protein_pairs.size();
	cout<< " There are total " <<num_of_jobs<<"  jobs generated!"<< endl << endl;



	//distribute the jobs to the slave processors
    int result, total_proc;
    int num_MPI, num_of_slaves, currentJobID;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &total_proc);

    num_of_slaves = total_proc - 1;
    num_MPI = (num_of_jobs <= num_of_slaves) ? num_of_jobs : num_of_slaves; /* number of slave MPI processes to start */


    for (int i = 1; i <= num_MPI; i++) {
        currentJobID = i - 1;
//        MPI_Send(&currentJobID, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
        string bla = hmm_protein_pairs.at(currentJobID);
        MPI_Send(&bla[0], bla.length(), MPI_CHAR, i, WORKTAG, MPI_COMM_WORLD);

    }

    if (num_of_jobs > num_of_slaves) {
        currentJobID = num_of_slaves;

        while (currentJobID < num_of_jobs) {
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            string bla = hmm_protein_pairs.at(currentJobID);
            MPI_Send(&bla[0], bla.length(), MPI_CHAR, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);

            currentJobID++;
        }
    }
    // all zeros for dying command
    int die = 0;
    for (int i = 1 ; i <= num_of_slaves; i++) {
        MPI_Send(&die, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
    }

}


void SlaveProcess(const string & file_dir_hmm,const string & file_dir_proteins, const string & out_directory, const string & numcpu)
{
    MPI_Status status;
    int myid;
    int number_amount;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    while (true) {

    	MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG) {
            cout << "   " << myid << ": I'm going to die " << endl;
            break;
        }

        MPI_Get_count(&status, MPI_CHAR, &number_amount);
        /* do the pipeline process for all the groups in this chunk */


        char *buf = new char[number_amount];
        MPI_Recv(buf, number_amount, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        string pairstring(buf, number_amount);
        delete [] buf;

        vector<string> paras = explode(pairstring, '\t');

        string hmmfile = paras[0];
        int maxhmmsize = std::stoi(paras[1]);
        int minhmmsize = std::stoi(paras[2]);
        string proteinfile = paras[3];
        int maxproteinsize = std::stoi(paras[4]);
        int minproteinsize = std::stoi(paras[5]);



        int maxrange = maxhmmsize*(1+0.2);
        int minrange = minhmmsize*(1-0.2);
        int hmmrange = maxrange-minrange+1;
        int proteinrange = maxproteinsize - minproteinsize +1;
        int overallmin = (minrange<minproteinsize)?minrange:minproteinsize;
        int overallmax = (maxrange>maxproteinsize)?maxrange:maxproteinsize;
        int overall = overallmax-overallmin+1;
        bool doflag = false;
        if(hmmrange+proteinrange>overall)//check if two ranges having overlap
        	doflag = true;



        if(doflag)
        {

//        string hmmsearchPath="/home/cjg/bin/hmmsearch"; //viper
        string hmmsearchPath="/lustre/atlas1/bip108/proj-shared/Unifam/bin/hmmsearch"; //titan
//        string hmmsearchPath="/ccs/home/chaij1/hmmsearch";//titan cannot see this node
//        string hmmsearchPath="/lustre/pfs1/cades-bsd/world-shared/qy2/Unifam/bin/hmmsearch"; //condo
        string Eval= "0.0001";
        string cpu = numcpu;
        string database = file_dir_hmm+Utils::getPathSeparator()+hmmfile;
        string protein = file_dir_proteins+Utils::getPathSeparator()+proteinfile;
        string databaseprefix = Utils::getFilename2(hmmfile);
        string proteinprefix = Utils::getFilename2(proteinfile);
        string thisoutdir = out_directory+Utils::getPathSeparator()+proteinprefix;

        Utils::mkdirIfNonExist(thisoutdir);

        string output_domtbfile = thisoutdir+Utils::getPathSeparator()+proteinprefix+"_"+databaseprefix+".domtab";

			if(Utils::isFileExist(output_domtbfile)==false||Utils::isFileEmpty(output_domtbfile)==true)
			{

			string hmmsearch_cmd = hmmsearchPath + " -Z 100000 -E " + Eval +" --noali --cpu " + cpu + " -o /dev/null --domtblout " + output_domtbfile + " " + database + " " + protein + "> /dev/null";

			cout<<"process# "<<myid<<" is currently doing "<<hmmfile<<" and "<<proteinfile<<endl;

			const int res = system(hmmsearch_cmd.c_str());
			}
			else
			{
				cout<<"process# "<<myid<<" does NOTNEED to do "<<hmmfile<<" and "<<proteinfile<<endl;
			}
        }
        else
        {
        	 cout<<"process# "<<myid<<" is NOTNOT doing "<<hmmfile<<" and "<<proteinfile<<endl;
        }

        MPI_Send(0,0,MPI_INT,0,0,MPI_COMM_WORLD);
    }
}

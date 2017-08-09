#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include <dirent.h>
#include <iostream>
#include <string>

#include "utils.hpp"

#define WORKTAG    1
#define DIETAG     2

using namespace std;
// slave process, everybody is slave...
void MasterProcess(vector<string>& folders, int num_of_jobs);
void SlaveProcess(const string & folder_for_protein_hmmresult, const string & proteinfasta_dir, const string & annotationdatabase_dir, const string & databasetype, const string & finaloutputfolder);
/* print usage */
void usage()
{
    cout << " [Usage]" << endl
    << "  unifam_mma [options] -hd <folder_for_protein_hmmresult> -fd <proteinfasta_dir> -data <annotation_database> -dtype <database_type> -od <final_output_folder>" << endl

    << endl
	<< "This is a MPI version of unifam to generate group file and annotations for each protein chunk after hmmsearch result"<<endl
	<< "usually run after 'unifam_mm'"<<endl
    << " [Inputs]" << endl
	<< " Please don't put '\\' or '/' after the string of the directory" << endl
	<< " -hd <folder_for_protein_hmmresult> this is folder of folders of protein chunks containing hmmsearch result *.domtab files" << endl
	<< " -fd <proteinfasta_dir> protein sequence fasta files folder" << endl
	<< " -data <annotation_database> this is the folder for annotation database downloaded from Uniprot" << endl
	<< " -dtype <database_type> annotation type, eg. 'prok',''euk', 'all'" << endl

    << endl
    << " Outputs of this program are:" << endl
    << "   annotation file, group file, protein fasta file with annotated header" << endl
	<< "  They will be copied to the folder defined by [output]:" << endl
    << "   -od <final_output_folder> directory for final annotation and other related files" << endl
    << endl
    << " [Options]" << endl
    << "   -h/--help Display this help message" << endl
    << endl;
}

bool initializeArguments(int argc, char **argv,
                         string & folder_for_protein_hmmresult,
						 string & proteinfasta_dir,
						 string & annotation_database,
						 string & database_type,
						 string & final_output_folder
)
{
    vector<string> Arguments;
    folder_for_protein_hmmresult = "";
    proteinfasta_dir = "";
    annotation_database = "";
    database_type = "";
    final_output_folder = "";
    int paranum = 0;



    while(argc--)
        Arguments.push_back(*argv++);

    for(int i = 1; i <= (int)Arguments.size()-1; i++)
    {

        if (Arguments[i] == "-hd")
        {
        	folder_for_protein_hmmresult = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-fd")
        {
        	proteinfasta_dir = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-data")
        {
        	annotation_database = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-dtype")
        {
        	database_type = Arguments[++i];
        	paranum++;
        }
        else if (Arguments[i] == "-od")
        {
        	final_output_folder = Arguments[++i];
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
    if (paranum <5)
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

    string folder_for_protein_hmmresult, proteinfasta_dir,annotation_database,database_type,final_output_folder;

    int num_of_folders=0;
    vector<string> folders;

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

    successflag = initializeArguments(argc, argv,folder_for_protein_hmmresult, proteinfasta_dir,annotation_database,database_type,final_output_folder);
    if(successflag == false)
    	{
    	cout << "[STOP Running here]" << endl;
    	return 0;
    	}

    if(!Utils::isDirectory(folder_for_protein_hmmresult))
    {
    	cout<< "[ERROR]::-hd is not set correctly!" << endl;
    	return 0;
    }
    else if(!Utils::isDirectory(proteinfasta_dir))
    {
    	cout<< "[ERROR]::-fd is not set correctly!" << endl;
    	return 0;
    }
    else if(!Utils::isDirectory(annotation_database))
    {
    	cout<< "[ERROR]::-data is not set correctly!" << endl;
    	return 0;
    }

    Utils::mkdirIfNonExist(final_output_folder);

    /*get the folder list from the directory*/
	folders.clear();

	DIR *d;
	struct dirent *dir;
	d = opendir(folder_for_protein_hmmresult.c_str());
	if (d)
	{
	  while ((dir = readdir(d)) != NULL)
	  {
		  if (dir->d_type == DT_DIR)
		  {
			  string foldname( dir->d_name );
			  std::size_t found=foldname.find('.');
			    if (found==std::string::npos)
			    	folders.push_back( foldname );
		  }
	  }

	  closedir(d);
	}

	num_of_folders = folders.size();




    /* echo the arguments */
    if (myid == 0) {

        cout<< " The arguments are taken!" << endl << endl;
    }

    /********************** generate data for each node, and find chunk information *********************/

    if(num_of_nodes > 1){
	    if (myid == 0) {
		    MasterProcess(folders, num_of_folders);
		    cout << "  Running with total " << num_of_nodes << " processes" << endl;
	    }

	    else
		    SlaveProcess(folder_for_protein_hmmresult, proteinfasta_dir, annotation_database,database_type,final_output_folder);
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

void MasterProcess(vector<string>& folders, int num_of_jobs)
{
	//distribute the jobs to the slave processors
    int result, total_proc;
    int num_MPI, num_of_slaves, currentJobID;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &total_proc);

    num_of_slaves = total_proc - 1;
    num_MPI = (num_of_jobs <= num_of_slaves) ? num_of_jobs : num_of_slaves; /* number of slave MPI processes to start */


    for (int i = 1; i <= num_MPI; i++) {
        currentJobID = i - 1;
        string bla = folders.at(currentJobID);
        MPI_Send(&bla[0], bla.length(), MPI_CHAR, i, WORKTAG, MPI_COMM_WORLD);

    }

    if (num_of_jobs > num_of_slaves) {
        currentJobID = num_of_slaves;

        while (currentJobID < num_of_jobs) {
            MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            string bla = folders.at(currentJobID);
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



void SlaveProcess(const string & folder_for_protein_hmmresult, const string & proteinfasta_dir, const string & annotationdatabase_dir, const string & databasetype, const string & finaloutputfolder)
{
    MPI_Status status;
    int myid;
    int number_amount;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    while (true)
    {
    	MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG)
        {
            cout << "   " << myid << ": I'm going to die " << endl;
            break;
        }

        MPI_Get_count(&status, MPI_CHAR, &number_amount);
        /* do the pipeline process for all the groups in this chunk */


        char *buf = new char[number_amount];
        MPI_Recv(buf, number_amount, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        string pairstring(buf, number_amount);
        delete [] buf;



//        string unifampath = "/chongle/qiuming/unifam_test/src/UniFam.py"; //Viper
//        string pythonpath = "/home/cjg/bin/python"; //viper
        string unifampath = "/lustre/atlas1/bif107/proj-shared/Unifam/unifam/src/UniFam.py"; //Titan
        string pythonpath = "python"; //Titan


		string proteinchunk_name=pairstring;
		string proteinchunk_dir = folder_for_protein_hmmresult+Utils::getPathSeparator()+proteinchunk_name;



        string config_str="\n\
[hmmsearch]\n\
output = /dev/null\n\
hmmsearchPath = hmmsearch\n\
cpu = 1\n\
eval = 0.0001\n\
\n\
[UniFam]\n\
seqCoverage = 0.5\n\
hmmCoverage = 0.5\n\
dataDir= "+annotationdatabase_dir+Utils::getPathSeparator()+"\n\
doParse = True\n\
dohmmsearch = False\n\
database = "+databasetype+"\n\
inputFormat = proteins\n\
doProdigal = False\n\
doPathway = False\n\
doRNAmmer = False\n\
dotRNAscan = False\n\
\n\
name = "+proteinchunk_name
        +"\nworkDir = "+proteinchunk_dir
        +"\ntmpDir = /dev/null\n";


        string configfile = proteinchunk_dir+ Utils::getPathSeparator() + proteinchunk_name+".cfg";
    	ofstream filePointer;
    	filePointer.open(configfile.c_str());
    	filePointer<<config_str;
    	filePointer.close();

        string proteinfile = proteinfasta_dir+ Utils::getPathSeparator() + proteinchunk_name+".fasta";



        string combine_file = "cat "+proteinchunk_dir+ Utils::getPathSeparator() +"*.domtab > "+proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+".domtab";
        system(combine_file.c_str());

        string unifam_cmd = pythonpath+ " "+unifampath+" -i "+proteinfile+" -c "+ configfile;
        const int res = system(unifam_cmd.c_str());
        if (res != 0 )
            cout<<"*** Error: Failed command: " + unifam_cmd<<endl;
        else
        {
        	string move_command1 = "mv "+proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+".domtab"+" "+finaloutputfolder+Utils::getPathSeparator();
        	system(move_command1.c_str());
        	string move_command2 = "mv "+proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+".annot"+" "+finaloutputfolder+Utils::getPathSeparator();
        	system(move_command2.c_str());
        	string move_command3 = "mv "+proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+".group"+" "+finaloutputfolder+Utils::getPathSeparator();
        	system(move_command3.c_str());
        	string move_command4 = "mv "+proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+"_annot.faa"+" "+finaloutputfolder+Utils::getPathSeparator();
        	system(move_command4.c_str());
        }


        /*cleanout*/
        string file1 = proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+".domtab";
        Utils::ifFileExistRemove(file1);
        string file2 = proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+".annot";
        Utils::ifFileExistRemove(file2);
        string file3 = proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+".group";
        Utils::ifFileExistRemove(file3);
        string file4 = proteinchunk_dir+ Utils::getPathSeparator()+proteinchunk_name+"_annot.faa";
        Utils::ifFileExistRemove(file4);

        MPI_Send(0,0,MPI_INT,0,0,MPI_COMM_WORLD);
    }
}

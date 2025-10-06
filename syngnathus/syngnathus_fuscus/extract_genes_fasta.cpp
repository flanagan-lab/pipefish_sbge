//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com
//Start Date: 30 April 2016
//Last Updated: 30 April 2016
//Keep only some fasta entries based on a list...but uses less memory than subset_fasta_file


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#define INT_MAX 2147483647

using namespace std;


bool FileTest(ifstream& file, string filename)
{
	cout << filename;
	if (file.is_open())
		cout << " open\n";
	else
	{
		while (!file.is_open())
		{
			cout << " not open. Please re-enter filename: ";
			getline(cin, filename, '\n');
			file.open(filename);
		}
	}
	return true;
}

istream& universal_getline(istream& is, string& t)
{
	//this code is adapted from a post on stackoverflow:
	// http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	//written by user763305
	t.clear();
	istream::sentry se(is, true);
	streambuf* sb = is.rdbuf();//sets pointer to stream buffer object

	for (;;)
	{
		int c = sb->sbumpc();//get current character and advance to the next position
		switch (c)//tests for equality against a list of variables (like multiple if statements)
		{
		case '\n'://if the next character is '\n', return the line
			return is;
		case '\r'://if the character is '\n', see if the next one is '\n' or return the line
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:// Also handle the case when the last line has no line ending
			if (t.empty())//if there's nothing there, set it to be the end of file
				is.setstate(ios::eofbit);//set it to be the end of the file and return it
			return is;
		default://if none of the above, continue on.
			t += (char)c;
		}
	}

}

class fasta_record
{
public:
	string seq_id;
	string sequence;

	fasta_record()
	{
		seq_id = "";
		sequence = "";
	}
};


int main(int argc, char* argv[])
{
	int end, i, j, num_ids, num_records;
	string fasta_in_name, fasta_out_name, list_name, name, seq,line, out_name, prefix;
	ifstream fasta_in, list_in;
	ofstream fasta_out;
	bool interactivemode, combine, found, simplify_name;
	string query, tempstring1, tempstring2;
	vector <string> scaffold_names;
	combine = false;
	simplify_name = false;
	prefix = "";

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nextract_genes_fasta:\n";
			cout << "Remove specific genes from a fasta file based on a list\n";
			cout << "-f:\tfasta file input\n";
			cout << "-l:\tList of Fasta Record IDs\n";
			cout << "-o:\tLabel to be added on end of input file.\n";
			cout << "-c:\tIf included, it means that all records in list should be Combined into one file. Otherwise, one file per fasta ID.\n";
			cout << "-s:\tIf included, the header names will be simplified to create more reasonable file names. Will not have an effect in combination with -c.\n";
			cout << "-h:\tPrints this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
		if (query == "I" || query == "i")
			interactivemode = true;
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nextract_genes_fasta:\n";
			cout << "Remove specific genes from a fasta file based on a list\n";
			cout << "-f:\tfasta file input\n";
			cout << "-l:\tList of Fasta Record IDs\n";
			cout << "-o:\tLabel to be added on the start of input file.\n";
			cout << "-c:\tIf included, it means that all records in list should be Combined into one file. Otherwise, one file per fasta ID.\n";
			cout << "-s:\tIf included, the header names will be simplified to create more reasonable file names. Will not have an effect in combination with -c.\n";
			cout << "-h:\tPrints this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc; i++)
	{
		interactivemode = false;
		tempstring1 = argv[i];
		//tempstring2 = argv[i + 1];
		if (tempstring1 == "-f")
			fasta_in_name = argv[i + 1];
		if (tempstring1 == "-l")
			list_name = argv[i + 1];
		if (tempstring1 == "-c")
			combine = true;
		if (tempstring1 == "-o")
			prefix = argv[i + 1];
		if(tempstring1 == "-s")
			simplify_name = true;
	}

	if (interactivemode)
	{
		cout << "Please provide input fasta file name\n";
		cin >> fasta_in_name;
		cout << "Please provide name of file with list of fasta record IDs\n";
		cin >> list_name;
		cout << "Please provide a string to add to the start of the input file name.\n";
		cin >> prefix;
		//cout << "Provide output path\n";
		//cin >> fasta_out_name;
		cout << "Do you want all of the records in the list to be combined into one file? Default is one file per fasta ID. Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "Y" || tempstring1 == "y")
		{
			combine = true;
		}
		if(combine != true){
			cout << "Do you want to simplify headers for more reliable filenames? Y or N\n";
			cin >> tempstring1;
			if(tempstring1 == "Y" || tempstring1=="y" || tempstring1=="yes" || tempstring1 == "Yes")
			{
				simplify_name = true;
			}
		}
	}

	if (combine)
		fasta_out_name = fasta_in_name.substr(0, fasta_in_name.size() - 6) + prefix + ".fasta";
	cout << "\nInput file:\t" << fasta_in_name;
	cout << "\nList file:\t" << list_name;
	if (combine)
		cout << "\nOne file will be output, called " << fasta_out_name << ".\n";
	else
		cout << "\nMultiple files will be created, one per fasta ID.\n";
	if (simplify_name)
		cout << "\nHeaders will be simplified for creating the file names.\n";
	if (interactivemode)
	{
		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;

		if (query != "y" && query != "Y")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> i;
			return 0;
		}
	}
	else
		cout << "\n\nProceeding...\n";

	list_in.open(list_name);
	FileTest(list_in, list_name);
	while (!list_in.eof())
	{
		if (!list_in.eof())
		{
			universal_getline(list_in, line);
			cout << line << '\t';
			scaffold_names.push_back(line);
		}
	}
	list_in.close();
	num_ids = scaffold_names.size();
	cout << "\nFound " << num_ids << " fasta IDs.\n";

	

	fasta_in.open(fasta_in_name);
	FileTest(fasta_in, fasta_in_name);
	int count = 0;
	num_records = 0;
	if (combine)
	{
		cout << "\nWriting to new file, " << fasta_out_name << '\n';
		fasta_out.open(fasta_out_name);
	}
	while (!fasta_in.eof())
	{
		if (!fasta_in.eof())
		{
			
			while(universal_getline(fasta_in, line))
			{
				if(line.empty())
					continue;
				if(line[0]=='>')
				{
					if(found == true)
					{
						fasta_out << seq;
						//cout << seq.length() << " bases have been written to file:\n" << seq << '\n';
						if(!combine){
							fasta_out.close();
							//cout << fasta_out_name << " has been closed.\n";
						}
						name.clear();	
					}
					found = false;
					name = line.substr(1, line.length());
					seq.clear();
					for (i = 0; i < scaffold_names.size(); i++)
					{
						if (scaffold_names[i] == name)
						{
							cout << "\nMatched " << name << " to " << scaffold_names[i] << ", and i=" << i <<'\n';
						
							found = true;
				
							if (combine && count == 0)
								fasta_out << ">" << name << '\n' << seq;
							if (combine && count > 0)
								fasta_out << "\n>" << name << '\n' << seq;
							if (!combine)
							{
								if(simplify_name)
								{
									string short_name = scaffold_names[i].substr(0, scaffold_names[i].find(' '));
									fasta_out_name = prefix + short_name + ".fasta";
								}
								else 
								{
									fasta_out_name = prefix + scaffold_names[i] + ".fasta";
								}
								cout << "Writing to new file, " << fasta_out_name << '\n';
								fasta_out.open(fasta_out_name);
								fasta_out << ">" << name << '\n';// << seq;
							}
							count++;
						}
					}
					num_records++;
				}
				else
				{	
					if(found == true)
					{
						seq += line;
					}
				}
			}
		}
	}
	fasta_in.close();
	if (combine)
		fasta_out.close();
	cout << '\n'<< fasta_in_name << " had " << num_records << " sequences.\n";
	cout << count << " records were matched and saved to file.\n";
	
	if (interactivemode)
	{
		cout << "\nDone! Input integer to quit: ";
		cin >> end;
	}
	return 0;
}


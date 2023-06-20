#include <iostream>
#include <string>

#ifndef READ_PARSESCRIPT_CONFIG_H
#define READ_PARSESCRIPT_CONFIG_H

class Configfile
{
	//Defining variables read in by the configuration file
	TCut m_globalcut="";
	int m_pass{0};
	int m_kine{4};
	int m_sbsfieldscale{30};
	TString m_target="LD2";
	TString m_output_dir="";
	TString m_outfilename_preint="";
	TString m_input_dir="";
	
public:

	// The following function reads the config file and returns 0 if the all the necessary paramers were specified and returns -1 if not.
	int readin_parsescript_configfile( const char* configfilename)
	{
		ifstream configfile(configfilename);
	  	TString currentline;

	  	std::cout <<'\n'<<"--- Reading configuration file: " << configfilename << " --- \n";

	  	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") )
	    {
	    	if( !currentline.BeginsWith("#") )
	    	{
	      		m_globalcut += currentline;
	      		std::cout << "Global cuts requested: " << m_globalcut << '\n';
	    	}    
	  	}

	  	std::string globalcut_string = m_globalcut.GetTitle();

	  	if ( globalcut_string.compare("") == 0 )
	  	{
	  		std::cerr << "Error: The input gobal cut is empty!\n";
	  		return -1;
	  	}

	  	while( currentline.ReadLine( configfile ) )
	  	{
	  		if ( !currentline.BeginsWith("#") )
	  		{
		  		TObjArray *tokens = currentline.Tokenize(" ");
		    	int ntokens = tokens->GetEntries();
		    	if( ntokens>1 )
		    	{	
		      		TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
		      		
		      		if( skey == "pass" )
		      		{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_pass = sval.Atoi();
						std::cout << "Pass: " << m_pass << endl;
		      		}

		      		else if( skey == "kine" )
		      		{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_kine = sval.Atoi();
						std::cout << "SBS kine number: " << m_kine << endl;
		      		}

		      		else if( skey == "sbsfieldscale" )
		      		{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_sbsfieldscale = sval.Atoi();
						std::cout << "SBS fieldscale: " << m_sbsfieldscale << endl;
		      		}

		      		else if( skey == "target" )
		      		{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_target = sval;
						std::cout << "Target: " << m_target << endl;
		      		}

		      		else if( skey == "output_dir" )
		      		{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_output_dir = sval;
						std::cout << "Output directory: " << m_output_dir << endl;
		      		}

		      		else if( skey == "outfilename_preint" )
		      		{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_outfilename_preint = sval;
						std::cout << "Output filename preint: " << m_outfilename_preint << endl;
		      		}

		      		else if( skey == "input_dir" )
		      		{
						TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
						m_input_dir = sval;
						std::cout << "Input mass replayed ROOT file directory: " << m_input_dir << endl;
		      		}

		      		delete tokens;
		    	}
		    	else if ( ntokens==1 )
		    	{
		    		TString skey = ( (TObjString*)(*tokens)[0] )->GetString();

		    		if( skey == "pass" )
		      		{
						std::cerr << "Error: Pass number not provied! \n ";
						return -1;
		      		}

		      		else if( skey == "kine" )
		      		{
						std::cerr << "Error: SBS kinematic number not provided! \n ";
						return -1;
		      		}

		      		else if( skey == "sbsfieldscale" )
		      		{
						std::cerr << "Error: SBS magnet fieldscale not provided! Program stoppig.\n";
						return -1;
		      		}

		      		else if( skey == "target" )
		      		{
						std::cerr << "Error : Target not provided! \n";
						return -1;
		      		}

		      		else if( skey == "output_dir" )
		      		{
						std::cerr << "Erro: Output directory not provied! \n";
						return -1;
		      		}

		      		else if( skey == "input_dir" )
		      		{
		      			m_input_dir = Form("/work/halla/sbs/sbs-gmn/pass%i/SBS%i/%s/rootfiles", m_pass, m_kine, m_target.Data());
						// If the above does not work, try volatile disk.
						//m_input_dir = Form("/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", m_pass, m_kine, m_target.Data());
						std::cout << "Warning: Input mass replayed ROOT file dir not provided. Using defaul dir: " << m_input_dir <<'\n';
					}
		    	}
	  		}
	  		
	    }

	    return 0;
	}

	TCut return_globalcut()
	{
		return m_globalcut;
	}

	int return_pass_num()
	{
		return m_pass;
	}

	int return_kin_num()
	{
		return m_kine;
	}

	int return_sbsfieldscale()
	{
		return m_sbsfieldscale;
	}

	TString return_target()
	{
		return m_target;
	}

	TString return_outputdir()
	{
		return m_output_dir;
	}

	TString return_outfilename_preint()
	{
		return m_outfilename_preint;
	}

	TString return_inputdir()
	{
		return m_input_dir;
	}
};

#endif
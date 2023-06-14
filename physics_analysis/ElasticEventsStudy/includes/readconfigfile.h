#ifndef READCONFIGFILE_H
#define READCONFIGFILE_H

//Defining variables read in by the configuration file
//TChain *C = new TChain("T"); //Initialize the TChain to chain the root files for analysis.
TCut globalcut{""};
double W2_min{0.};
double W2_max{0.};
double vzcut{0.};
double pse_min{0.};
double eoverp_leftcut{0.};
double eoverp_rightcut{0.};
double hcal_energycut{0.};
double bbcal_energycut{0.};
double trigdiff_leftcut{0.};
double trigdiff_rightcut{0.};
double hcal_clusblk_ADCtime_leftcut{0.};
double hcal_clusblk_ADCtime_rightcut{0.};
double mc_Neventsgen{0.};
double mc_Ibeam_uA{0.};

void readin_anaconfigfile( const char* configfilename, TChain* C )
{
	ifstream configfile(configfilename);
  	TString currentline;

  	std::cout <<'\n'<<"--- Reading configuration file " << configfilename << " --- \n";

  	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") )
  	{
    	if( !currentline.BeginsWith("#") )
    	{
      	C->Add(currentline);
      	cout << "Loaded root file: " << currentline << '\n';
      	}
    }

  	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") )
    {
    	if( !currentline.BeginsWith("#") )
    	{
      		globalcut += currentline;
      		cout << "Global cuts requested: " << globalcut << '\n';
    	}    
  	}

  	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") )
  	{
  		TObjArray *tokens = currentline.Tokenize(" ");
    	int ntokens = tokens->GetEntries();
    	if( ntokens>1 )
    	{	
      		TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      		/*if( skey == "rootfile_path")
      		{
      			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      			TString runnum_string{std::to_string(runnum)};
      			C->Add(sval+"e1209019_fullreplay_"+runnum_string+"_stream0_seg*_*.root");
      			cout << "All the root files of run number "<< runnum << " were loaded for the analysis" <<'\n';
      		}*/
      		if( skey == "W2_min" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			W2_min = sval.Atof();
			cout << "W2 min (GeV): " << W2_min << endl;
      		}

      		if( skey == "W2_max" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			W2_max = sval.Atof();
			cout << "W2 max (GeV): " << W2_max << endl;
      		}

      		if( skey == "vzcut" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			vzcut = sval.Atof();
			cout << "vzcut (m): " << vzcut << endl;
      		}

      		if( skey == "pse_min" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			pse_min = sval.Atof();
			cout << "pse_min(GeV): " << pse_min << endl;
      		}

      		if( skey == "eoverp_leftcut" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			eoverp_leftcut = sval.Atof();
			cout << "eoverp_leftcut: " << eoverp_leftcut << endl;
      		}

      		if( skey == "eoverp_rightcut" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			eoverp_rightcut = sval.Atof();
			cout << "eoverp_rightcut: " << eoverp_rightcut << endl;
      		}

      		if( skey == "hcal_energycut" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			hcal_energycut = sval.Atof();
			cout << "hcal_energycut (m): " << hcal_energycut << endl;
      		}

      		if( skey == "bbcal_energycut" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			bbcal_energycut = sval.Atof();
			cout << "bbcal_energycut (GeV): " << bbcal_energycut << endl;
      		}

      		if( skey == "trigdiff_leftcut" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			trigdiff_leftcut = sval.Atof();
			cout << "trigdiff_leftcut (ns): " << trigdiff_leftcut << endl;
      		}

      		if( skey == "trigdiff_rightcut" )
      		{
			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
			trigdiff_rightcut = sval.Atof();
			cout << "trigdiff_rightcut (ns): " << trigdiff_rightcut << endl;
      		}

      		if( skey == "hcal_clusblk_ADCtime_leftcut" )
      		{
      			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
				hcal_clusblk_ADCtime_leftcut = sval.Atof();
				cout << "hcal_clusblk_ADCtime_leftcut (ns): " << hcal_clusblk_ADCtime_leftcut << endl;
      		}

      		if( skey == "hcal_clusblk_ADCtime_rightcut" )
      		{
      			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
				hcal_clusblk_ADCtime_rightcut = sval.Atof();
				cout << "hcal_clusblk_ADCtime_rightcut (ns): " << hcal_clusblk_ADCtime_rightcut << endl;
      		}

      		if( skey == "mc_Neventsgen")
      		{
      			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
				mc_Neventsgen = sval.Atof();
				cout << "Number of events simulated: " << mc_Neventsgen << endl;
      		}

      		if( skey == "mc_Ibeam_uA")
      		{
      			TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
				mc_Ibeam_uA = sval.Atof();
				cout << "Beam current used in the simulation: " << mc_Ibeam_uA << endl;
      		}

      		delete tokens;
    	}
    }

}

#endif